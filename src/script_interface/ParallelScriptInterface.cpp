/*
 * Copyright (C) 2015-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ParallelScriptInterface.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/utility.hpp>

#include "RemoteObjectHandle.hpp"

#include <cassert>

namespace {
Communication::MpiCallbacks *m_cb = nullptr;

void make_remote_handle() { new ScriptInterface::RemoteObjectHandle(m_cb); }
} // namespace

REGISTER_CALLBACK(make_remote_handle)

namespace ScriptInterface {
ParallelScriptInterface::ParallelScriptInterface(std::string const &name)
    : m_callback_id(m_cb, [](CallbackAction) {}) {
  assert(m_cb && "Not initialized!");

  /* Create the slaves */
  m_cb->call(make_remote_handle);

  /* Create local object */
  m_p = ObjectHandle::make_shared(name, ObjectHandle::CreationPolicy::LOCAL);
}

ParallelScriptInterface::~ParallelScriptInterface() {
  /* Delete the instances on the other nodes */
  try {
    call(CallbackAction::DELETE);
  } catch (...) {
    /* We have to do MPI calls in the destructor, which
       can cause exceptions. To avoid sanitizer warnings
       about this, and to clarify we directly explicitly
       call terminate, which would happen anyway.
    */
    std::terminate();
  }
}

void ParallelScriptInterface::initialize(Communication::MpiCallbacks &cb) {
  m_cb = &cb;
}

void ParallelScriptInterface::do_construct(VariantMap const &params) {
  /* Bcast class name and global id to the slaves */
  std::pair<ObjectId, std::string> what = std::make_pair(m_p->id(), name());
  boost::mpi::broadcast(m_cb->comm(), what, 0);

  call(CallbackAction::CONSTRUCT);

  auto p = unwrap_variant_map(params);
  boost::mpi::broadcast(m_cb->comm(), p, 0);

  m_p->construct(p);
}

void ParallelScriptInterface::do_set_parameter(const std::string &name,
                                               const Variant &value) {
  std::pair<std::string, Variant> d(name, Variant());

  if (is_type<ObjectId>(value)) {
    d.second = map_parallel_to_local_id(value);
  } else {
    d.second = value;
  }

  call(CallbackAction::SET_PARAMETER);

  boost::mpi::broadcast(m_cb->comm(), d, 0);

  m_p->set_parameter(d.first, d.second);

  collect_garbage();
}

Variant ParallelScriptInterface::do_call_method(const std::string &name,
                                                const VariantMap &parameters) {
  call(CallbackAction::CALL_METHOD);
  VariantMap p = unwrap_variant_map(parameters);

  auto d = std::make_pair(name, p);
  /* Broadcast method name and parameters */
  boost::mpi::broadcast(m_cb->comm(), d, 0);

  auto ret = map_local_to_parallel_id(m_p->call_method(name, p));

  collect_garbage();

  return ret;
}

Variant ParallelScriptInterface::get_parameter(std::string const &name) const {
  auto p = m_p->get_parameter(name);

  return map_local_to_parallel_id(p);
}

VariantMap ParallelScriptInterface::unwrap_variant_map(VariantMap const &map) {
  /* Copy parameters into a non-const buffer, needed by boost::mpi */
  auto p = map;

  /* Unwrap the object ids */
  for (auto &it : p) {
    if (is_type<ObjectId>(it.second)) {
      it.second = map_parallel_to_local_id(it.second);
    }
  }

  return p;
}

Variant
ParallelScriptInterface::map_local_to_parallel_id(Variant const &value) const {
  if (is_type<ObjectId>(value)) {
    /** Check if the objectid is the empty object (ObjectId()),
     * if so it does not need translation, the empty object
     * has the same id everywhere.
     */
    auto oid = get_value<ObjectId>(value);

    if (oid != ObjectId()) {
      return obj_map.at(oid)->id();
    }
    return oid;
  }
  if (is_type<std::vector<Variant>>(value)) {
    auto const &in_vec = get_value<std::vector<Variant>>(value);
    std::vector<Variant> out_vec;
    out_vec.reserve(in_vec.size());

    for (auto const &e : in_vec) {
      out_vec.emplace_back(map_local_to_parallel_id(e));
    }

    return out_vec;
  }
  return value;
}

Variant
ParallelScriptInterface::map_parallel_to_local_id(Variant const &value) {
  const auto outer_id = get_value<ObjectId>(value);

  auto so_ptr = get_instance(outer_id).lock();

  auto po_ptr = std::dynamic_pointer_cast<ParallelScriptInterface>(so_ptr);

  if (po_ptr != nullptr) {
    auto inner_id = po_ptr->get_underlying_object()->id();

    /* Store a pointer to the object */
    obj_map[inner_id] = po_ptr;

    /* and return the id of the underlying object */
    return inner_id;
  }
  if (so_ptr == nullptr) {
    /* Release the object */
    obj_map.erase(outer_id);

    /* Return None */
    return ObjectId();
  }
  throw std::runtime_error(
      "Parameters passed to Parallel entities must also be parallel.");
}

void ParallelScriptInterface::collect_garbage() {
  /* Removal condition, the instance is removed iff
     its payload is not used anywhere. In this case
     the reference count is one, because the host object
     still holds a pointer.
  */
  auto pred = [](map_t::value_type const &e) -> bool {
    return e.second->get_underlying_object().use_count() == 1;
  };

  for (auto it = obj_map.begin(); it != obj_map.end();) {
    if (pred(*it)) {
      obj_map.erase(it);
    }
    ++it;
  }
}
} /* namespace ScriptInterface */
