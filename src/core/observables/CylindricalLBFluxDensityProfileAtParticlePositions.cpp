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
#include "CylindricalLBFluxDensityProfileAtParticlePositions.hpp"

#include "BoxGeometry.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"

#include "communication.hpp"

#include <utils/Histogram.hpp>
#include <utils/Span.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <boost/mpi/collectives/gather.hpp>

#include <vector>

namespace Observables {
std::vector<double>
CylindricalLBFluxDensityProfileAtParticlePositions::evaluate(
    ParticleReferenceRange const &local_particles,
    const ParticleObservables::traits<Particle> &traits) const {
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

  Utils::CylindricalHistogram<double, 3> histogram(n_bins(), limits());
  // First collect all positions (since we want to call the LB function to
  // get the fluid velocities only once).

  for (auto const &pos_vec : global_folded_positions) {
    for (auto const &pos : pos_vec) {
      auto const v = LB::get_interpolated_velocity(pos) * LB::get_lattice_speed();
      auto const flux_dens = LB::get_interpolated_density(pos) * v;

      histogram.update(Utils::transform_coordinate_cartesian_to_cylinder(
                          pos - transform_params->center(),
                          transform_params->axis(),
                          transform_params->orientation()),
                      Utils::transform_vector_cartesian_to_cylinder(
                          flux_dens, transform_params->axis(),
                          pos - transform_params->center()));
    }
  }

  // normalize by number of hits per bin
  auto hist_tmp = histogram.get_histogram();
  auto tot_count = histogram.get_tot_count();
  std::transform(hist_tmp.begin(), hist_tmp.end(), tot_count.begin(),
                 hist_tmp.begin(), [](auto hi, auto ci) {
                   return ci > 0 ? hi / static_cast<double>(ci) : 0.;
                 });
  return hist_tmp;
}
} // namespace Observables
