/*
 * Copyright (C) 2016-2022 The ESPResSo project
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

#ifndef OBSERVABLES_PIDOBSERVABLE_HPP
#define OBSERVABLES_PIDOBSERVABLE_HPP

#include <particle_observables/observable.hpp>

#include "Observable.hpp"
#include "Particle.hpp"
#include "ParticleTraits.hpp"
#include "communication.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/flatten.hpp>

#include <boost/range/algorithm/copy.hpp>

#include <cstddef>
#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace Observables {

using ParticleReferenceRange =
    std::vector<std::reference_wrapper<const Particle>>;

/** %Particle-based observable.
 *
 *  Base class for observables extracting raw data from particle subsets and
 *  returning either the data or a statistic derived from it.
 */
class PidObservable : virtual public Observable {
  /** Identifiers of particles measured by this observable */
  std::vector<int> m_ids;

  virtual std::vector<double>
  evaluate(ParticleReferenceRange const &particles,
           const ParticleObservables::traits<Particle> &traits) const = 0;

public:
  explicit PidObservable(std::vector<int> ids) : m_ids(std::move(ids)) {}
  std::vector<double> operator()() const final;
  std::vector<int> const &ids() const { return m_ids; }
};

namespace detail {
/**
 * Recursive implementation for finding the shape of a given `std::vector` of
 * types. A vector of extents is constructed starting at
 * the template specialization for `std::vector<T>`.
 */
template <class T> struct shape_impl;

template <> struct shape_impl<double> {
  static std::vector<std::size_t> eval(std::size_t /* n_part */) { return {1}; }
};
template <class _, std::size_t N> struct shape_impl<Utils::Vector<_, N>> {
  static std::vector<std::size_t> eval(std::size_t /* n_part */) { return {N}; }
};
template <class T> struct shape_impl<std::vector<T>> {
  static std::vector<std::size_t> eval(std::size_t n_part) {
    std::vector<std::size_t> ret{n_part};
    boost::copy(shape_impl<T>::eval(n_part), std::back_inserter(ret));

    return ret;
  }
};
template <class T, class U> struct shape_impl<std::pair<T, U>> {
  static std::vector<std::size_t> eval(std::size_t n_part) {
    return shape_impl<T>::eval(n_part);
  }
};
} // namespace detail

/**
 * This class implements an interface to the `particle_observables` library that
 * implements necessary algorithms needed for observables that are based on
 * single particle properties.
 * @tparam ObsType An observables composed of an algorithm from
 * src/particle_observables/include/particle_observables/algorithms.hpp and two
 * particle properties.
 *
 *  Example usage:
 *  @code{.cpp}
 *  using namespace ParticleObservables;
 *  using CenterOfMass = ParticleObservable<WeightedAverage<Position, Mass>>;
 *  @endcode
 */
template <class ObsType> class ParticleObservable : public PidObservable {
public:
  using PidObservable::PidObservable;
  std::vector<std::size_t> shape() const override {
    using std::declval;

    return detail::shape_impl<decltype(declval<ObsType>()(
        declval<ParticleReferenceRange const &>()))>::eval(ids().size());
  }

  template <typename T> struct is_map : std::false_type {};
  template <typename T>
  struct is_map<ParticleObservables::Map<T>> : std::true_type {};

  std::vector<double>
  evaluate(ParticleReferenceRange const &particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    if constexpr (is_map<ObsType>::value) {
      std::vector<double> local_traits;
      Utils::flatten(ObsType{}(particles), std::back_inserter(local_traits));
      std::vector<int> local_pids;
      Utils::flatten(ParticleObservables::Identities{}(particles),
                     std::back_inserter(local_pids));
      auto const pid_begin = std::begin(local_pids);
      auto const pid_end = std::end(local_pids);

      auto const n_dims = local_traits.size() / local_pids.size();
      std::vector<double> output;
      output.reserve(local_traits.size());
      for (auto const pid : ids()) {
        auto const pid_pos = std::find(pid_begin, pid_end, pid);
        auto const i =
            static_cast<std::size_t>(std::distance(pid_begin, pid_pos));
        for (std::size_t j = 0; j < n_dims; ++j) {
          output.emplace_back(local_traits[i * n_dims + j]);
        }
      }
      return output;
    } else {
      auto const local_result = ObsType{}(particles);
      std::remove_const_t<decltype(local_result)> result;

      boost::mpi::reduce(
          comm_cart, local_result, result,
          ObsType::template reduction<decltype(std::declval<ObsType>()(
              std::declval<ParticleReferenceRange const &>()))>,
          0);

      return result.first;
    }
  }
};

} // namespace Observables
#endif
