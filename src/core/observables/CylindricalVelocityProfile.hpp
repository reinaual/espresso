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
#ifndef OBSERVABLES_CYLINDRICALVELOCITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALVELOCITYPROFILE_HPP

#include "BoxGeometry.hpp"
#include "CylindricalPidProfileObservable.hpp"
#include "grid.hpp"

#include <utils/Histogram.hpp>
#include <utils/Span.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <boost/range/combine.hpp>

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace Observables {
class CylindricalVelocityProfile : public CylindricalPidProfileObservable {
public:
  using CylindricalPidProfileObservable::CylindricalPidProfileObservable;

  std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    using pos_type = Utils::Vector3d;
    using vel_type = Utils::Vector3d;

    std::vector<pos_type> local_folded_positions;
    local_folded_positions.reserve(local_particles.size());
    std::vector<vel_type> local_velocities;
    local_velocities.reserve(local_particles.size());

    for (auto const &p : local_particles) {
      auto const pos = folded_position(traits.position(p), box_geo) -
                       transform_params->center();
      local_folded_positions.emplace_back(
          Utils::transform_coordinate_cartesian_to_cylinder(
              pos, transform_params->axis(), transform_params->orientation()));
      local_velocities.emplace_back(
          Utils::transform_vector_cartesian_to_cylinder(
              traits.velocity(p), transform_params->axis(), pos));
    }

    std::vector<std::vector<pos_type>> global_folded_positions;
    std::vector<std::vector<vel_type>> global_velocities;
    boost::mpi::gather(comm, local_folded_positions, global_folded_positions,
                       0);
    boost::mpi::gather(comm, local_velocities, global_velocities, 0);

    if (comm.rank() != 0) {
      return {};
    }

    Utils::CylindricalHistogram<double, 3> histogram(n_bins(), limits());

    for (auto const &[pos_vec, vel_vec] :
         boost::combine(global_folded_positions, global_velocities)) {
      for (auto const &[pos, vel] : boost::combine(pos_vec, vel_vec)) {
        histogram.update(pos, vel);
      }
    }

    auto hist_tmp = histogram.get_histogram();
    auto tot_count = histogram.get_tot_count();
    for (std::size_t ind = 0; ind < hist_tmp.size(); ++ind) {
      if (tot_count[ind] > 0) {
        hist_tmp[ind] /= static_cast<double>(tot_count[ind]);
      }
    }
    return hist_tmp;
  }

  std::vector<std::size_t> shape() const override {
    auto const b = n_bins();
    return {b[0], b[1], b[2], 3};
  }
};

} // Namespace Observables

#endif
