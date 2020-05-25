#
# Copyright (C) 2017-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import unittest as ut
import espressomd
import espressomd.observables
import numpy as np


class RdfTest(ut.TestCase):
    s = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.s.box_l = 3 * [10]
        self.s.part.clear()

    def bin_volumes(self, midpoints):
        bin_size = midpoints[1] - midpoints[0]
        r = midpoints - 0.5 * bin_size
        r2 = np.power(r, 2)

        # Volumes of the bins
        return 4. * np.pi * (r2 * bin_size + r *
                             bin_size**2 + bin_size**3 / 3.)

    def test_single_type(self):
        s = self.s

        n_part = 99
        dx = self.s.box_l[0] / float(n_part + 1)

        for i in range(n_part):
            s.part.add(
                id=i, pos=[i * dx, 0.5 * s.box_l[1], 0.5 * s.box_l[2]], type=0)

        r_bins = 50
        r_min = 0.5 * dx
        r_max = r_bins * dx
        obs = espressomd.observables.RDF(ids1=s.part[:].id, min_r=r_min,
                                         max_r=r_max, n_r_bins=r_bins)
        rdf = obs.calculate()
        r = obs.bin_centers()
        rv = self.bin_volumes(r)
        rho = n_part / (s.box_l[0]**3)

        parts_in_bin = rdf * rv * rho

        # All but the last bin should contain 2 particles
        np.testing.assert_allclose(parts_in_bin[:-1], 2.0, rtol=1e-1)

    def test_mixed(self):
        s = self.s

        n_part = 99
        dx = self.s.box_l[0] / float(n_part + 1)

        for i in range(n_part):
            s.part.add(
                id=i, pos=[i * dx, 0.5 * s.box_l[1], 0.5 * s.box_l[2]], type=(i % 2))

        r_bins = 50
        r_min = 0.5 * dx
        r_max = r_bins * dx
        obs = espressomd.observables.RDF(ids1=s.part[:].id[0::2],
                                         ids2=s.part[:].id[1::2],
                                         min_r=r_min, max_r=r_max,
                                         n_r_bins=r_bins)
        rdf01 = obs.calculate()

        r = obs.bin_centers()
        rv = self.bin_volumes(r)
        rho = 0.5 * n_part / (s.box_l[0]**3)
        parts_in_bin = rdf01 * rv * rho

        # Every even bin should contain two parts
        np.testing.assert_allclose(parts_in_bin[0::2], 2.0, rtol=1e-1)
        # Every odd bin should contain zero parts
        np.testing.assert_allclose(parts_in_bin[1::2], 0.0)

        # Check symmetry
        obs = espressomd.observables.RDF(ids1=s.part[:].id[1::2],
                                         ids2=s.part[:].id[0::2],
                                         min_r=r_min, max_r=r_max,
                                         n_r_bins=r_bins)
        rdf10 = obs.calculate()

        np.testing.assert_allclose(rdf10, rdf01)


if __name__ == "__main__":
    ut.main()
