from unittest import TestCase
import math
import time

import numpy as np
import pandas
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

from main_package.simulation.WaveThing import RSFWaveThing, SpaceProperties, RSFProperties
from main_package.simulation.InitialConditioners import RSFInitialConditioner_PacketType, \
    RSFInitialConditioner_HiggsTypeMatchedToRSFAABBInteraction, \
    RSFInitialConditioner_HiggsTypeMatchedToFullyGeneralInteraction
from main_package.simulation.Interaction import RSFAABBInteraction


class TestRSFInitialConditioner_PacketType(TestCase):
    mSq = 3
    quarticTerm = 4

    @staticmethod
    def naive_putPacket(wave: RSFWaveThing, propRSF: RSFProperties, propSpace: SpaceProperties, centralP: float,
                        pSpread: float, centralX: float, norm: float):
        """Algorithm used in the original Java version. Used for testing purposes."""

        def addModeWithMomIndexKn(kn: int, xnCen: int, norm: float):
            piKnOnN = propSpace.piOnN * kn
            twoPiKnOnN = piKnOnN * 2.
            erm = math.sin(piKnOnN) * 2. / propSpace.delta
            en = math.sqrt(erm * erm + propRSF.mSq)

            for xn in range(propSpace.n):
                arg = twoPiKnOnN * (xn - xnCen)
                psiReal = math.cos(arg) * norm
                psiImag = -math.sin(arg) * norm

                psiDotReal = -en * psiImag

                wave.data.x[xn] += psiReal
                wave.data.v[xn] += psiDotReal

        bottomIp = int(centralP * propSpace.oneOverDp - 0.5 * propSpace.n)
        topIp = bottomIp + propSpace.n

        xnCen = int(centralX // propSpace.delta)

        for ip in range(bottomIp, topIp):
            pThatMultX = propSpace.dp * ip
            pDistSigs = (pThatMultX - centralP) / pSpread
            pDistSigsSq = pDistSigs ** 2
            weight = math.exp(-pDistSigsSq * 0.5) / (math.sqrt(2. * math.pi) * pSpread) * propSpace.dp
            addModeWithMomIndexKn(ip, xnCen, weight * norm)

    def test_put_packet_X_data_comparison_with_data_generated_by_Python(self):
        ###########################################################
        # Create wave packet to be tested (using FFT)

        n = 1001
        L = 60

        propSpace = SpaceProperties(n, L)

        pCentral = 4
        pSpread = 1
        centralX = 5
        normScale = 6

        propRSF_A = RSFProperties(mSq=self.mSq, quarticTerm=self.quarticTerm)
        rsfWaveThing_A = RSFWaveThing(propSpace, propRSF_A)

        initialConditioner = RSFInitialConditioner_PacketType(
            rsfWaveThing_A,
            pCentral,
            pSpread,
            centralX,
            normScale)

        tic = time.perf_counter()
        # initialConditioner.putPacket(rsfWaveThing_A, pCentral, pSpread, centralX, normScale)
        initialConditioner.setInitialCondition()
        toc = time.perf_counter()
        print('Execution time: ', toc - tic, ' seconds')

        x = rsfWaveThing_A.data.x

        ###########################################################
        # Create reference wave packet (using the original algorithm)

        propRSF_test = RSFProperties(mSq=self.mSq, quarticTerm=self.quarticTerm)
        propSpace_test = SpaceProperties(n, L)
        rsfWaveThing_test = RSFWaveThing(propSpace_test, propRSF_test)
        self.naive_putPacket(rsfWaveThing_test, propRSF_test, propSpace_test, pCentral, pSpread, centralX, normScale)

        x_test = rsfWaveThing_test.data.x

        ###########################################################
        # Test if the wave packets agree within floating point error
        diff = np.abs(x - x_test)
        self.assertTrue(np.all(diff < 1e-04))

        # fig1, ax1 = plt.subplots(figsize=(8, 6), dpi=300)
        # ax1.plot(x, ':', label='Python', linewidth=0.5, )
        # ax1.plot(x_test, label='Java', linewidth=0.5)
        # ax1.legend()
        # fig1.show()

        # fig2, ax2 = plt.subplots(figsize=(8, 6), dpi=300)
        # ax2.plot(x-x_test)
        # # ax2.set_ylim([-1,1])
        # ax2.legend()
        # fig2.show()

        # ax1.plot(x - x_test, linewidth=0.5)
        # # ax2.set_ylim([-1,1])
        # ax1.legend()
        # fig1.show()

    def test_put_packet_X_data_comparison_with_data_generated_by_Java(self):
        ###########################################################
        # Create wave packet to be tested (using FFT)

        n = 1001
        L = 60

        propSpace = SpaceProperties(n, L)

        pCentral = 4
        pSpread = 1
        centralX = 5
        normScale = 6

        propRSF_A = RSFProperties(mSq=self.mSq, quarticTerm=self.quarticTerm)
        rsfWaveThing_A = RSFWaveThing(propSpace, propRSF_A)

        initialConditioner = RSFInitialConditioner_PacketType(
            rsfWaveThing_A,
            pCentral,
            pSpread,
            centralX,
            normScale)

        # initialConditioner.putPacket(rsfWaveThing_A, pCentral, pSpread, centralX, normScale)
        initialConditioner.setInitialCondition()

        x = rsfWaveThing_A.data.x

        ###########################################################
        # Get reference wave packet data (data was generated using the original Java version of the simulation)

        data = pandas.read_csv('data/initialisation.csv')
        xData: np.ndarray = data['x0'].values

        ###########################################################
        # Test if the wave packets agree within floating point error
        diff = np.abs(x - xData)
        self.assertTrue(np.all(diff < 1e-04))

        # print(np.argmax(x - xData))
        # # self.assertTrue(np.array_equal(x, xData))
        #
        # print(len(x))
        # fig1, ax1 = plt.subplots(figsize=(8, 6), dpi=300)
        # ax1.plot(rsfWaveThing_A.data.x, label='Python')
        # ax1.plot(data['x0'], label='Java')
        # ax1.legend()
        # fig1.show()
        #
        # for i in range(-10, 10):
        #     fig2, ax2 = plt.subplots(figsize=(8, 6), dpi=300)
        #     ax2.plot(np.roll(rsfWaveThing_A.data.x, i) - data['x0'], label=str(i))
        #     ax2.set_ylim([-1, 1])
        #     ax2.legend()
        #     fig2.show()

        # fig2, ax2 = plt.subplots(figsize=(8, 6), dpi=300)
        # ax2.plot(data['x0'])
        # ax2.set_title('Java')
        # fig2.show()

        # plt.plot(data['x0'])
        # plt.title('Java')
        # print(len(data['x0']))
        # plt.show()
        # plt.figure(figsize=(8, 6), dpi=300)

    def test_put_packet_V_data_comparison_with_data_generated_by_Python(self):
        ###########################################################
        # Create wave packet to be tested (using FFT)

        n = 1001
        L = 60

        propSpace = SpaceProperties(n, L)

        pCentral = 4
        pSpread = 1
        centralX = 5
        normScale = 6

        propRSF_A = RSFProperties(mSq=self.mSq, quarticTerm=self.quarticTerm)
        rsfWaveThing_A = RSFWaveThing(propSpace, propRSF_A)

        initialConditioner = RSFInitialConditioner_PacketType(
            rsfWaveThing_A,
            pCentral,
            pSpread,
            centralX,
            normScale)

        tic = time.perf_counter()
        # initialConditioner.putPacket(rsfWaveThing_A, pCentral, pSpread, centralX, normScale)
        initialConditioner.setInitialCondition()
        toc = time.perf_counter()
        print('Execution time: ', toc - tic, ' seconds')

        v = rsfWaveThing_A.data.v

        ###########################################################
        # Create reference wave packet (using the original algorithm)

        propRSF_test = RSFProperties(mSq=self.mSq, quarticTerm=self.quarticTerm)
        propSpace_test = SpaceProperties(n, L)
        rsfWaveThing_test = RSFWaveThing(propSpace_test, propRSF_test)
        self.naive_putPacket(rsfWaveThing_test, propRSF_test, propSpace_test, pCentral, pSpread, centralX, normScale)

        v_test = rsfWaveThing_test.data.v

        ###########################################################
        # Test if the wave packets agree within floating point error
        diff = np.abs(v - v_test)
        self.assertTrue(np.all(diff < 1e-04))

        # fig1, ax1 = plt.subplots(figsize=(8, 6), dpi=300)
        # ax1.plot(v, label='Python')
        # ax1.plot(v_test, label='Java')
        # ax1.plot(diff, label='Difference')
        # ax1.legend()
        # fig1.show()
        # print(np.max(diff))

    def test_put_packet_V_data_comparison_with_data_generated_by_Java(self):
        ###########################################################
        # Create wave packet to be tested (using FFT)

        n = 1001
        L = 60

        propSpace = SpaceProperties(n, L)

        pCentral = 4
        pSpread = 1
        centralX = 5
        normScale = 6

        propRSF_A = RSFProperties(mSq=self.mSq, quarticTerm=self.quarticTerm)
        rsfWaveThing_A = RSFWaveThing(propSpace, propRSF_A)

        initialConditioner = RSFInitialConditioner_PacketType(
            rsfWaveThing_A,
            pCentral,
            pSpread,
            centralX,
            normScale)

        # initialConditioner.putPacket(rsfWaveThing_A, pCentral, pSpread, centralX, normScale)
        initialConditioner.setInitialCondition()

        v = rsfWaveThing_A.data.v

        ###########################################################
        # Get reference wave packet data (data was generated using the original Java version of the simulation)

        data = pandas.read_csv('data/initialisation.csv')
        vData: np.ndarray = data['v0'].values

        ###########################################################
        # Test if the wave packets agree within floating point error
        diff = np.abs(v - vData)
        self.assertTrue(np.all(diff < 1e-04))

    def test_RSFAABBHiggsInteraction(self):
        ###########################################################
        # Create wave packets to be tested
        n = 1001
        L = 60

        propSpace = SpaceProperties(n, L)

        pCentral = 4
        pSpread = 1
        centralX = 5
        normScale = 6

        initialConditioners = []

        # A
        propRSF_A = RSFProperties(mSq=self.mSq, quarticTerm=self.quarticTerm)
        rsfWaveThing_A = RSFWaveThing(propSpace, propRSF_A)

        initialConditioners.append(RSFInitialConditioner_PacketType(
            rsfWaveThing_A,
            pCentral,
            pSpread,
            centralX,
            normScale))

        # C
        propRSF_C = RSFProperties(-2 * 2, 0.02 * 2 * 2)
        rsfWaveThing_C = RSFWaveThing(propSpace, propRSF_C)

        # Interaction
        rsfAABBInteraction_AC = RSFAABBInteraction(rsfWaveThing_A, rsfWaveThing_C, 0.1)

        initialConditioners.append(
            RSFInitialConditioner_HiggsTypeMatchedToRSFAABBInteraction(rsfWaveThing_C, rsfWaveThing_A,
                                                                       rsfAABBInteraction_AC))

        # Set initial conditions
        for i in initialConditioners:
            i.setInitialCondition()

        x = rsfWaveThing_C.getX()
        v = rsfWaveThing_C.getV()

        ###########################################################
        # Data from simulation
        data = pandas.read_csv('data/initialisation.csv')
        x_data: np.ndarray = data['x2'].values
        v_data: np.ndarray = data['v2'].values

        ###########################################################
        # Compare
        diff = np.abs(x - x_data)
        # fig1, ax1 = plt.subplots(figsize=(8, 6), dpi=300)
        # ax1.plot(x, label='Python')
        # ax1.plot(x_data, label='Java')
        # ax1.plot(diff, label='diff')
        # print(np.max(np.abs(diff)))
        # ax1.legend()
        # fig1.show()

        self.assertTrue(np.all(diff < 1e-10))


class TestRSFInitialConditioner_HiggsTypeMatchedToFullyGeneralInteraction(TestCase):
    def test_find_global_min_of_polynomial(self):
        p: Polynomial = Polynomial([0, 0, 1, 0, -5, 0, 5])
        domain: list[float, float] = [-10, 12]
        x_min_true: float = 0.737666

        x_min = RSFInitialConditioner_HiggsTypeMatchedToFullyGeneralInteraction.find_global_min_of_polynomial(p, domain)

        self.assertTrue(abs(x_min - x_min_true) < 1e-5)
