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
#ifndef OBSERVABLES_FLUXDENSITYPROFILE_HPP
#define OBSERVABLES_FLUXDENSITYPROFILE_HPP

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "PidProfileObservable.hpp"
#include "grid.hpp"

#include "communication.hpp"

#include <utils/Histogram.hpp>
#include <utils/Span.hpp>

#include <cstddef>
#include <vector>

namespace Observables {
class FluxDensityProfile : public PidProfileObservable {
public:
  using PidProfileObservable::PidProfileObservable;
  std::vector<std::size_t> shape() const override {
    auto const b = n_bins();
    return {b[0], b[1], b[2], 3};
  }

  std::vector<double>
  evaluate(ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    using pos_type = decltype(traits.position(std::declval<Particle>()));

    std::vector<pos_type> local_folded_positions;
    local_folded_positions.reserve(local_particles.size());

    for (auto const &p : local_particles) {
      local_folded_positions.emplace_back(folded_position(traits.position(p), box_geo));
    }

    std::vector<std::vector<pos_type>> global_folded_positions;
    boost::mpi::gather(comm_cart, local_folded_positions, global_folded_positions, 0);

    if (comm_cart.rank() != 0) {
      return {};
    }

    Utils::Histogram<double, 3> histogram(n_bins(), limits());

    for (auto const &vec : global_folded_positions) {
      for (auto const &p : vec) {
        histogram.update(p);
      }
    }

    histogram.normalize();
    return histogram.get_histogram();
  }
};

} // Namespace Observables

#endif
