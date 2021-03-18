/*
 * Copyright (C) 2021 The ESPResSo project
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

#define BOOST_TEST_MODULE elc test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <utils/Vector.hpp>

#include "electrostatics_magnetostatics/elc.hpp"

BOOST_AUTO_TEST_CASE(mirror_z_position) {
  const Utils::Vector3d position{1, 2, 3};
  auto const lower_mirror = detail::mirror_z_position(position, 0);
  const Utils::Vector3d expected_lower{position[0], position[1], -position[2]};
  BOOST_TEST(lower_mirror == expected_lower, boost::test_tools::per_element());

  auto const mirror_equal = detail::mirror_z_position(position, position[2]);
  BOOST_TEST(mirror_equal == position, boost::test_tools::per_element());

  auto const upper_mirror_position = 10.;
  auto const upper_mirror =
      detail::mirror_z_position(position, upper_mirror_position);
  const Utils::Vector3d expected_upper_mirror = {
      position[0], position[1], 2 * upper_mirror_position - position[2]};
  BOOST_TEST(upper_mirror == expected_upper_mirror,
             boost::test_tools::per_element());
}