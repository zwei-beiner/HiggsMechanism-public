import math
from unittest import TestCase

import numpy as np
import pandas

from main_package.simulation.InitialConditioners import RSFInitialConditioner_PacketType, InitialConditioner, \
    RSFInitialConditioner_HiggsTypeMatchedToRSFAABBInteraction
from main_package.simulation.Integrator import LeapfrogIntegrator
from main_package.simulation.Interaction import RSFAABBInteraction
from main_package.simulation.WaveThing import RSFProperties, SpaceProperties, RSFWaveThing


class TestRSFProperties(TestCase):
    def test_default_init(self):
        result = RSFProperties()

        self.assertEqual(result.mSq, 0)
        self.assertEqual(result.quarticTerm, 0)

    def test_init(self):
        mSq = 3.
        quarticTerm = 4.

        result = RSFProperties(mSq, quarticTerm)

        self.assertEqual(result.mSq, mSq)
        self.assertEqual(result.quarticTerm, quarticTerm)

    def test_is_sometimes_height_of_pot_min(self):
        def V(msq, q):
            x = math.sqrt(-msq / q)
            return 0.5 * msq * (x ** 2) + 0.25 * q * (x ** 4)

        test_cases = [[-3., 4., V(-3., 4.)], [0., 4., 0.], [3., 4., 0.],
                      [-3., 0., 0.], [0., 0., 0.], [3., 0., 0.],
                      [-3., -4., 0.], [0., -4., 0.], [3., -4., 0.]]

        for case in test_cases:
            result = RSFProperties(case[0], case[1]).isSometimesHeightOfPotMin()
            self.assertAlmostEqual(result, case[2], places=7)

    def test_get_quadratic_mass_about_vev_or_zero(self):
        test_cases = [[-3., 4., math.sqrt(-2. * (-3))], [0., 4., 0.], [3., 4., math.sqrt(3)],
                      [-3., 0., 0.], [0., 0., 0.], [3., 0., math.sqrt(3)],
                      [-3., -4., 0.], [0., -4., 0.], [3., -4., 0.]]

        for case in test_cases:
            result = RSFProperties(case[0], case[1]).getQuadraticMassAboutVevOrZero()
            self.assertAlmostEqual(result, case[2], places=7)

    def test_non_neg_vev_or_zero(self):
        test_cases = [[-3., 4., math.sqrt(3. / 4.)], [0., 4., 0.], [3., 4., 0.],
                      [-3., 0., 0.], [0., 0., 0.], [3., 0., 0.],
                      [-3., -4., 0.], [0., -4., 0.], [3., -4., 0.]]

        for case in test_cases:
            result = RSFProperties(case[0], case[1]).nonNegVevOrZero()
            self.assertAlmostEqual(result, case[2], places=7)

    ##########################################

    def test_non_neg_vev_returns0_1(self):
        mSq = 3.
        quarticTerm = 4.

        result = RSFProperties(mSq, quarticTerm).nonNegVev()

        self.assertEqual(result, 0)

    def test_non_neg_vev_returns0_2(self):
        mSq = 3.
        quarticTerm = 0

        result = RSFProperties(mSq, quarticTerm).nonNegVev()

        self.assertEqual(result, 0)

    def test_non_neg_vev_returns0_3(self):
        mSq = 0
        quarticTerm = 4.

        result = RSFProperties(mSq, quarticTerm).nonNegVev()

        self.assertEqual(result, 0)

    def test_non_neg_vev_UFB_1(self):
        mSq = 3.
        quarticTerm = -4.

        result = RSFProperties(mSq, quarticTerm)
        with self.assertRaises(Exception) as context:
            result.nonNegVev()
        self.assertEqual('UFB', str(context.exception))

    def test_non_neg_vev_UFB_2(self):
        mSq = .0
        quarticTerm = -4.

        result = RSFProperties(mSq, quarticTerm)
        with self.assertRaises(Exception) as context:
            result.nonNegVev()
        self.assertEqual('UFB', str(context.exception))

    def test_non_neg_vev_UFB_3(self):
        mSq = -3.
        quarticTerm = 0

        result = RSFProperties(mSq, quarticTerm)
        with self.assertRaises(Exception) as context:
            result.nonNegVev()
        self.assertEqual('UFB', str(context.exception))

    def test_non_neg_vev_UFB_4(self):
        mSq = -3.
        quarticTerm = -4.

        result = RSFProperties(mSq, quarticTerm)
        with self.assertRaises(Exception) as context:
            result.nonNegVev()
        self.assertEqual('UFB', str(context.exception))

    def test_non_neg_vev_FLAT(self):
        mSq = .0
        quarticTerm = .0

        result = RSFProperties(mSq, quarticTerm)
        with self.assertRaises(Exception) as context:
            result.nonNegVev()
        self.assertEqual('FLAT', str(context.exception))

    def test_non_neg_vev_returns(self):
        mSq = -3.
        quarticTerm = 4.

        result = RSFProperties(mSq, quarticTerm).nonNegVev()

        self.assertEqual(result, math.sqrt(- mSq / quarticTerm))


class TestRSFWaveThing(TestCase):

    @staticmethod
    def createData() -> tuple[np.ndarray, np.ndarray]:
        initialConditioners: list[InitialConditioner] = []
        integrator: LeapfrogIntegrator = LeapfrogIntegrator()

        rsfWaveThing: dict[str, RSFWaveThing] = {}

        ##################################################

        n = 1001
        L = 60

        propSpace = SpaceProperties(n, L)

        pCentral = 4
        pSpread = 1
        centralX = 5
        normScale = 6

        ##################################################

        rsfWaveThing['A'] = RSFWaveThing(propSpace, RSFProperties(mSq=0, quarticTerm=0))

        initialConditioners.append(RSFInitialConditioner_PacketType(
            rsfWaveThing['A'],
            pCentral,
            pSpread,
            centralX,
            normScale))
        integrator.add(rsfWaveThing['A'])

        ##################################################

        rsfWaveThing['B'] = RSFWaveThing(propSpace, RSFProperties(mSq=2.2 * 2.2, quarticTerm=0))

        initialConditioners.append(RSFInitialConditioner_PacketType(
            rsfWaveThing['B'],
            pCentral,
            pSpread,
            centralX,
            normScale))
        integrator.add(rsfWaveThing['B'])

        ##################################################

        rsfWaveThing['C'] = RSFWaveThing(propSpace, RSFProperties(mSq=-2 * 2, quarticTerm=0.02 * 2 * 2))
        rsfAABBInteraction_AC = RSFAABBInteraction(rsfWaveThing['A'], rsfWaveThing['C'], 0.1)

        initialConditioners.append(
            RSFInitialConditioner_HiggsTypeMatchedToRSFAABBInteraction(rsfWaveThing['C'], rsfWaveThing['A'],
                                                                       rsfAABBInteraction_AC))

        integrator.add(rsfWaveThing['C'])
        integrator.add(rsfAABBInteraction_AC)

        ##################################################

        rsfWaveThing['D'] = RSFWaveThing(propSpace, RSFProperties(mSq=0, quarticTerm=0))

        initialConditioners.append(
            RSFInitialConditioner_PacketType(rsfWaveThing['D'], pCentral, pSpread, centralX, normScale))
        integrator.add(rsfWaveThing['D'])

        ##################################################

        # Set initial conditions
        for initialConditioner in initialConditioners:
            initialConditioner.setInitialCondition()

        ##################################################

        shape = (500, 4, 1001)
        x: np.ndarray = np.zeros(shape)
        v: np.ndarray = np.zeros(shape)

        timeSteps: int = 500
        t: float = 0
        dt: float = 0.0008
        sorted_list_of_dic_elements = [v for k, v in sorted(rsfWaveThing.items(), key=lambda item: item[0])]
        for i_t in range(timeSteps):
            for wave, i_wave in zip(sorted_list_of_dic_elements, range(len(sorted_list_of_dic_elements))):
                x[i_t, i_wave, :] = np.copy(wave.data.x)
                v[i_t, i_wave, :] = np.copy(wave.data.v)

            integrator.advance(dt)
            t += dt

        return x, v

    def test_advance(self):
        ########################################
        # Create data to be tested
        x, v = self.createData()

        ########################################
        # Load reference data
        data = pandas.read_csv('data/runSimulationForFiniteTime.csv', comment='#',
                               dtype={'n': int, 'x0': float, 'v0': float, 'x1': float, 'v1': float, 'x2': float,
                                      'v2': float, 'x3': float, 'v3': float})
        x_test: np.ndarray = np.zeros(shape=(500, 4, 1001))
        v_test: np.ndarray = np.zeros(shape=(500, 4, 1001))

        for i_t in range(500):
            for i_wave in range(4):
                lower = 0 + i_t * 1001
                higher = 1001 + i_t * 1001
                x_test[i_t, i_wave, :] = data.iloc[lower:higher, [1 + i_wave * 2]].to_numpy().reshape((1001,))
                v_test[i_t, i_wave, :] = data.iloc[lower:higher, [2 + i_wave * 2]].to_numpy().reshape((1001,))

        ########################################
        # Compare
        diff_x = np.abs(x - x_test)
        diff_v = np.abs(v - v_test)

        print("Max value for x:", np.max(diff_x))
        print("Max index for x where (time step, wave number, coordinate index):", np.unravel_index(diff_x.argmax(), diff_x.shape))

        print("Max value for v:", np.max(diff_v))
        print("Max index for v where (time step, wave number, coordinate index):", np.unravel_index(diff_v.argmax(), diff_v.shape))

        self.assertTrue(np.all(diff_x < 1e-02))
        self.assertTrue(np.all(diff_v < 1e-01))

    def test_zero_the_wave(self):
        self.fail()

    def test_energy_at_index(self):
        self.fail()
