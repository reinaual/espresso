/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifndef FETCH_PARTICLES_HPP
#define FETCH_PARTICLES_HPP

#include "grid.hpp"
#include "particle_node.hpp"
#include "cells.hpp"

#include "PidObservable.hpp"
#include <utils/Span.hpp>

#include <algorithm>
#include <cstddef>
#include <vector>
#include <set>
#include <functional>

/** Fetch a group of particles.
 *
 *  @param ids particle identifiers
 *  @return array of particle copies, with positions in the current box.
 */
inline auto fetch_particles(std::vector<int> const &ids) {
  auto const ids_set = std::set<int>{ids.begin(), ids.end()};
  auto const local_particles = ::cell_structure.local_particles();
  Observables::ParticleReferenceRange local_particle_refs;
  std::copy_if(local_particles.begin(), local_particles.end(),
	       std::back_inserter(local_particle_refs),
	       [&ids_set](Particle &p) { return ids_set.count(p.id()) != 0; });
  return local_particle_refs;
}

/** Fetch a group of particles.
 *
 *  @param ids particle identifiers
 *  @return array of particle copies, with positions in the current box.
 */
inline std::vector<Particle> old_fetch_particles(std::vector<int> const &ids) {
  std::vector<Particle> particles;
  particles.reserve(ids.size());

  auto const chunk_size = fetch_cache_max_size();
  for (std::size_t offset = 0; offset < ids.size();) {
    auto const this_size = std::clamp(chunk_size, std::size_t{0},
                                      std::size_t{ids.size() - offset});
    auto const chunk_ids =
        Utils::make_const_span(ids.data() + offset, this_size);

    prefetch_particle_data(chunk_ids);

    for (auto id : chunk_ids) {
      particles.push_back(get_particle_data(id));

      auto &p = particles.back();
      p.pos() += image_shift(p.image_box(), box_geo.length());
      p.image_box() = {};
    }

    offset += this_size;
  }

  return particles;
}
#endif
