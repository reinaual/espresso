# Copyright (C) 2019 The ESPResSo project
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
from espressomd import electrostatic_extensions
from scipy.special import digamma

import matplotlib.pyplot as plt


@utx.skipIfMissingFeatures(["P3M"])
class ELC_vs_analytic(ut.TestCase):
    # Handle to espresso system
    box_l = 20.
    system = espressomd.System(box_l=[box_l, box_l, box_l])
    accuracy = 1e-7
    check_accuracy = 1e-4
    elc_gap = box_l
    system.time_step = 0.01
    delta_mid_top = -1.
    delta_mid_bot = -1.
    pot_diff = 0

    offset = 1
    number_samples = 25
    zPos = np.linspace(offset, box_l - offset, number_samples)
    q = np.array([-5])  # np.arange(-5.0, 5, 1)

    def test_elc_nonneutral(self):
        """
        Testing ELC against the analytic solution for an infinite system with dielectric contrast on top and bottom which is taken from:        
        Barrera, R. G., O. Guzman, and B. Balaguer. "Point charge in a three‐dielectric medium with planar interfaces." American Journal of Physics 46.11 (1978): 1172-1179 (http://www.fisica.unam.mx/personales/rbarrera/pdf/pub/int/8-AJP-46-1172-78.pdf)
        """
        self.system.box_l = [self.box_l, self.box_l, self.box_l + self.elc_gap]
        self.system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        self.system.periodicity = [1, 1, 1]

        self.system.part.add(id=1, pos=self.system.box_l / 2., q=self.q[0])

        p3m = espressomd.electrostatics.P3M(prefactor=1.,
                                            check_neutrality=False,
                                            accuracy=self.accuracy)
        self.system.actors.add(p3m)

        self.system.part.add(id=2, pos=self.system.box_l / 2., q=-self.q[0])
        self.move_charge_image(1, 2, self.box_l)

        p3m_results = self.scan(True)

        self.system.part[2].remove()

        elc = electrostatic_extensions.ELC(gap_size=self.elc_gap,
                                           maxPWerror=self.accuracy,
                                           delta_mid_bot=self.delta_mid_bot,
                                           delta_mid_top=self.delta_mid_top,
                                           check_neutrality=False,
                                           const_pot=True,
                                           pot_diff=self.pot_diff)

        self.system.actors.add(elc)

        elc_results = self.scan(False)

        # ANALYTIC SOLUTION
        charge_reshaped = np.square(self.q.reshape(-1, 1))

        # this is not periodic in x and y dimension
        zeta = (self.zPos - self.box_l / 2) / self.box_l
        analytic_energy = charge_reshaped / (4 * self.box_l) * (
                digamma(0.5 - zeta) + digamma(0.5 + zeta) - 2 * digamma(0.5))

        plt.figure()
        for i in range(len(self.q)):
            plt.plot(self.zPos, analytic_energy[i], label=f'analytic: q={self.q[i]}')
            plt.plot(self.zPos, elc_results[i, :, 1], label=f'ELC: q={self.q[i]}')
            plt.plot(self.zPos, p3m_results[i, :, 1], label=f'P3M: q={self.q[i]}')
        #          plt.plot(self.zPos, analytic_energy[i] - elc_results[i, :, 1], label=f'ELC: q={self.q[i]}')

        print(analytic_energy[i] - elc_results[i, :, 1])

        #        plt.plot(self.zPos, (elc_results[0, :, 1] + elc_results[-1, :, 1]) / 2, color='black', dashes=[2, 2])
        plt.legend()
        plt.show()

        analytic_force = charge_reshaped * self.zPos  # * (1 / self.distance ** 2 + self.delta_mid_bot * (
        #            1 / np.square(2 * self.zPos) - 1 / np.square(2 * self.zPos + self.distance)))

        analytic_results = np.dstack((analytic_force, analytic_energy))

        # np.testing.assert_allclose(
        #     elc_results[..., 1], analytic_results[..., 1], rtol=0, atol=self.check_accuracy)

    def scan(self, move: bool = False):
        result_array = np.empty((len(self.q), len(self.zPos), 2))
        for chargeIndex, charge in enumerate(self.q):
            self.system.part[1].q = charge
            for i, z in enumerate(self.zPos):
                pos = self.system.part[1].pos
                self.system.part[1].pos = [pos[0], pos[1], z]
                if move:
                    self.move_charge_image(1, 2, self.box_l)

                self.system.integrator.run(0)
                result_array[chargeIndex, i, 0] = self.system.part[1].f[2]
                result_array[chargeIndex, i, 1] = self.system.analysis.energy()["total"]  # / charge**2
        return result_array

    def move_charge_image(self, id1, id2, box_l):
        q = self.system.part[id1].q

        self.system.part[id2].q = -q

        pos = self.system.part[id1].pos
        if pos[2] < 0.5 * box_l:
            self.system.part[id2].pos = [pos[0], pos[1], -pos[2]]
        else:
            self.system.part[id2].pos = [pos[0], pos[1], 2 * box_l - pos[2]]


if __name__ == "__main__":
    ut.main()
