# Copyright (C) 2021 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np
import espressomd.electrostatics
import espressomd.electrostatic_extensions

def symmetry_testing(energies, forces, accuracy):
    # test symmetry of results
    lower_half = len(energies) // 2

    np.testing.assert_equal(np.argmax(energies), lower_half)
    np.testing.assert_almost_equal(forces[lower_half], 0, accuracy)

    # energy has to be symmetric to the central plane between the walls
    np.testing.assert_almost_equal(energies[:lower_half], energies[lower_half + 1:][::-1], accuracy)

    # force has to be symmetric to the central plane between the walls
    np.testing.assert_almost_equal(forces[:lower_half], -forces[lower_half + 1:][::-1], accuracy)

    np.testing.assert_equal(np.argmin(energies[:lower_half]), 0)
    np.testing.assert_equal(np.argmin(forces[:lower_half]), 0)


@utx.skipIfMissingFeatures(["P3M"])
class TEST2D_symmetry(ut.TestCase):
    # Handle to espresso system
    box_l = 20
    system = espressomd.System(box_l=[box_l, box_l, 2 * box_l])
    system.periodicity = [1, 1, 1]
    p3m_common_params = {"check_neutrality": False,
                         "accuracy": 1e-5,
                         "verbose": False}
    check_accuracy = 1e-12
    system.time_step = 0.01
    distance = 1

    delta_mid_bot = -1.
    delta_mid_top = -1.

    # number_samples has to be odd to fulfill symmetry condition
    number_samples = 21
    z_positions = np.linspace(distance, box_l - distance, number_samples)
    q = (-5.0, 5.0)

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def test_test2d_coulomb_prefactor(self):
        self.tearDown()
        p1 = self.system.part.add(pos=(self.box_l / 2, self.box_l / 2, self.box_l / 2), q=self.q[0])

        p3m = espressomd.electrostatics.P3M(prefactor=2.,
                                            **self.p3m_common_params)
        self.system.actors.add(p3m)

        test2d = espressomd.electrostatic_extensions.TEST2D(potential_difference=0)
        self.system.actors.add(test2d)

        energies_2, forces_2 = self.scan(p1)

        p3m.set_params(prefactor=1.)

        energies_1, forces_1 = self.scan(p1)

        np.testing.assert_almost_equal(energies_2, 2 * energies_1, self.check_accuracy)
        np.testing.assert_almost_equal(forces_2, 2 * forces_1, self.check_accuracy)

    def test_test2d_potential_difference(self):
        self.tearDown()
        p1 = self.system.part.add(pos=(self.box_l / 2, self.box_l / 2, self.box_l / 2), q=self.q[0])

        p3m = espressomd.electrostatics.P3M(prefactor=1.,
                                            **self.p3m_common_params)
        self.system.actors.add(p3m)

        test2d = espressomd.electrostatic_extensions.TEST2D(potential_difference=0)
        self.system.actors.add(test2d)

        energies_no_difference, forces_no_difference = self.scan(p1)

        potential_difference = 2.0

        test2d.set_params(potential_difference=potential_difference)

        energies_with_difference, forces_with_difference = self.scan(p1)

        np.testing.assert_almost_equal(energies_with_difference - energies_no_difference, self.q[0] * self.z_positions * potential_difference / self.box_l)
        np.testing.assert_almost_equal(forces_with_difference - forces_no_difference, -self.q[0] * np.full(self.z_positions.shape, potential_difference / self.box_l))

    def test_test2d_symmetry(self):
        self.tearDown()

        p1 = self.system.part.add(pos=(self.box_l / 2, self.box_l / 2, self.box_l / 2), q=self.q[0])

        p3m = espressomd.electrostatics.P3M(prefactor=2.,
                                            **self.p3m_common_params)
        self.system.actors.add(p3m)

        test2d = espressomd.electrostatic_extensions.TEST2D(potential_difference=0)
        self.system.actors.add(test2d)

        for q in self.q:
            p1.q = q

            energies, forces = self.scan(p1)
            symmetry_testing(energies, forces, self.check_accuracy)

    def scan(self, particle):
        energy_array = np.empty(len(self.z_positions))
        force_array = np.empty_like(energy_array)
        for i, z in enumerate(self.z_positions):
            pos = particle.pos
            particle.pos = [pos[0], pos[1], z]

            self.system.integrator.run(0)
            energy_array[i] = self.system.analysis.energy()["total"]
            force_array[i] = particle.f[2]
        return energy_array, force_array

if __name__ == "__main__":
    ut.main()
