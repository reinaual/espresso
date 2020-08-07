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


@utx.skipIfMissingFeatures(["P3M"])
class TEST2D_vs_analytic(ut.TestCase):
    # Handle to espresso system
    box_l = 200.
    system = espressomd.System(box_l=[box_l, box_l, 2 * box_l])
    accuracy = 1e-7
    check_accuracy = 1e-4
    system.time_step = 0.01
    distance = 0.5

    number_samples = 25
    zPos = np.linspace(distance, box_l - distance, number_samples)
    q = np.arange(-5.0, 5.1, 2.5)

    def test_test2d(self):
        self.system.part.add(id=1, pos=(self.box_l/2, self.box_l/2, self.box_l/2), q=self.q[0])
        self.system.part.add(id=2, pos=(self.box_l/2, self.box_l/2, 3*self.box_l/2), q=-self.q[0])

        self.system.box_l = [self.box_l, self.box_l, self.box_l]
        self.system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        self.system.periodicity = [1, 1, 1]
        p3m = espressomd.electrostatics.P3M(prefactor=1.,
                                            accuracy=self.accuracy)
        self.system.actors.add(p3m)
        
        test2d = espressomd.electrostatic_extensions.TEST2D()
        self.system.actors.add(test2d)
        
        self.system.part[2].remove()
        
        

#        elc = electrostatic_extensions.ELC(gap_size=self.elc_gap,
#                                           maxPWerror=self.accuracy,
#                                           delta_mid_bot=self.delta_mid_bot,
#                                           delta_mid_top=self.delta_mid_top)
#        self.system.actors.add(elc)

        elc_results = self.scan()
        
        import matplotlib.pyplot as plt
        plt.figure()
        plt.subplot(121)
        for i, charge in enumerate(self.q):
          plt.plot(self.zPos, elc_results[i, :, 0], label=f'{charge}')
        plt.title('force')
        
        plt.subplot(122)
        for i, charge in enumerate(self.q):
          plt.plot(self.zPos, elc_results[i, :, 1], label=f'{charge}')
        plt.title('energy')
        plt.legend()
        
        plt.show()
        
        

    def scan(self):
        result_array = np.empty((len(self.q), len(self.zPos), 2))
        for chargeIndex, charge in enumerate(self.q):
            self.system.part[1].q = charge
            for i, z in enumerate(self.zPos):
                pos = self.system.part[1].pos
                self.system.part[1].pos = [pos[0], pos[1], z]

                self.system.integrator.run(0)
                result_array[chargeIndex, i, 0] = self.system.part[1].f[2]
                result_array[chargeIndex, i, 1] = self.system.analysis.energy()[
                    "total"]
        return result_array


if __name__ == "__main__":
    ut.main()
