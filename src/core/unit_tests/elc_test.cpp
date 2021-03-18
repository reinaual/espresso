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

BOOST_AUTO_TEST_CASE(mirror_position) {
  const double position = 4;
  auto const lower_mirror = detail::mirror_position(position, 0);
  const double expected_lower = -position;
  BOOST_CHECK_EQUAL(lower_mirror, expected_lower);

  auto const mirror_equal = detail::mirror_position(position, position);
  BOOST_CHECK_EQUAL(mirror_equal, position);

  auto const upper_mirror_position = 10.;
  auto const upper_mirror =
      detail::mirror_position(position, upper_mirror_position);
  const double expected_upper_mirror = 2 * upper_mirror_position - position[2];
  BOOST_CHECK_EQUAL(upper_mirror, expected_upper_mirror);
}