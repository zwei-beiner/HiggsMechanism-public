import math
import sys

import numpy as np
from main_package.simulation.EnergyDensity import EnergyDensity


# Abstract class
from main_package.simulation.PlotOptions import PlotOptions


class LeapfroggableWaveThing:
    def zeroTheWave(self):
        raise NotImplementedError()

    def advanceVBasedOnX(self, dt: float):
        raise NotImplementedError()

    def advanceXBasedOnV(self, dt: float):
        raise NotImplementedError()


class RealData:
    def __init__(self, n):
        self.x = np.zeros(n, dtype=np.float64)
        self.v = np.zeros(n, dtype=np.float64)

    def setToZero(self):
        self.x.fill(0)
        self.v.fill(0)


class RSFProperties():
    """
    Stores the "diagonal" (i.e. self-interaction) energy terms: V(x)=1/2 m^2 x^2 + 1/4 q x^4
    Handles calculations related to the vacuum expectation value (vev) and the minimum of the potential term.

    Class uses Exception handling for control flow (if no valid value can be returned, a string is returned,
    which is done by throwing an exception containing the string).
    """

    def __init__(self, mSq=0, quarticTerm=0):
        self.mSq = mSq
        self.quarticTerm = quarticTerm

    def get_mSq(self) -> float:
        return self.mSq

    def set_mSq(self, mSq: float) -> None:
        self.mSq = mSq

    def get_quarticTerm(self) -> float:
        return self.quarticTerm

    def set_quarticTerm(self, quarticTerm: float) -> None:
        self.quarticTerm = quarticTerm

    def isSometimesHeightOfPotMin(self) -> float:
        """
        Returns height of potential minimum V(x_min)
        """
        try:
            x = self.nonNegVev()
            return 0.5 * self.mSq * x ** 2 + 0.25 * self.quarticTerm * x ** 4
        except Exception:
            return 0

    def getQuadraticMassAboutVevOrZero(self) -> float:
        msq = self.mSq
        q = self.quarticTerm
        if msq > 0:
            return math.sqrt(msq) if q >= 0 else 0
        elif msq == 0:
            return 0
        else:
            # msq < 0
            return math.sqrt(-2. * msq) if q > 0 else 0

    def nonNegVevOrZero(self) -> float:
        """
        Called when a numerical value of the potential minimum is needed for calculations.
        The method calls nonNegVev but catches any exception and returns 0 instead.
        """
        try:
            return self.nonNegVev()
        except Exception:
            return 0

    def nonNegVev(self):
        """
        Returns x position of potential minimum V(x)=1/2 m^2 x^2 + 1/4 q x^4, if it exists uniquely.
        Otherwise, the error messages UFB (unbounded from below) or FLAT are raised.
        """
        msq = self.mSq
        q = self.quarticTerm
        if msq > 0:
            if q >= 0:
                return 0
            else:
                raise Exception('UFB')
        elif msq == 0:
            if q > 0:
                return 0
            elif q == 0:
                raise Exception('FLAT')
            else:
                raise Exception('UFB')
        else:
            # msq < 0
            if q > 0:
                return math.sqrt(-msq / q)
            else:
                raise Exception('UFB')


# TODO: Possibly make this a static class?
class SpaceProperties:
    def __init__(self, n: int, L: float):
        self.n = n
        self.L = L
        self.piOnN = math.pi / n
        self.delta = L / n
        self.maxMom = math.pi * 2. * n / L
        self.dp = math.pi * 2. / L
        self.oneOverDp = 1. / self.dp


class RSFWaveThing(LeapfroggableWaveThing, EnergyDensity):
    """
    Handles a single wave.
    """
    def __init__(self, propSpace: SpaceProperties, propRSF: RSFProperties, plot_options: PlotOptions = None):
        if plot_options is None:
            self._plot_options = PlotOptions()
        else:
            self._plot_options = plot_options

        self.propSpace = propSpace
        self.propRSF = propRSF
        self.data = RealData(propSpace.n)

        self.dTwoPhiByDXSquared = np.zeros(self.propSpace.n, dtype=np.float64)
        self.vDot = np.zeros(self.propSpace.n, dtype=np.float64)

    def get_plot_options(self) -> PlotOptions:
        return self._plot_options

    def get_propRSF(self) -> RSFProperties:
        return self.propRSF

    def get_N(self) -> int:
        return self.propSpace.n

    def getX(self) -> np.ndarray:
        return self.data.x

    def setX(self, new_x: np.ndarray) -> None:
        self.data.x = new_x

    def incrementX(self, dx: np.ndarray) -> None:
        self.data.x += dx

    def getV(self) -> np.ndarray:
        return self.data.v

    def setV(self, new_v: np.ndarray) -> None:
        self.data.v = new_v

    def incrementV(self, dv: np.ndarray) -> None:
        self.data.v += dv


    def advanceVBasedOnX(self, dt: float) -> None:
        # Alias for readability
        x = self.getX()

        # Calculate Phi'' = ( Phi_(i+1) + Phi_(i-1) - 2 Phi_i ) / delta^2
        # self.dTwoPhiByDXSquared[1:-1] = x[2:] + x[:-2] - 2 * x[1:-1]
        # self.dTwoPhiByDXSquared[0] = x[1] + x[-1] - 2 * x[0]
        # self.dTwoPhiByDXSquared[-1] = x[0] + x[-2] - 2 * x[-1]
        # self.dTwoPhiByDXSquared /= self.propSpace.delta ** 2
        self.dTwoPhiByDXSquared[:] = np.diff(x, 2, append=x[0], prepend=x[-1]) / self.propSpace.delta ** 2

        # Calculate PhiDotDot = Phi'' - m^2 Phi - q Phi^3
        self.vDot[:] = self.dTwoPhiByDXSquared - self.propRSF.mSq * x - self.propRSF.quarticTerm * (x ** 3)

        self.incrementV(dt * self.vDot)


    def advanceXBasedOnV(self, dt: float) -> None:
        self.incrementX(dt * self.getV())

    def zeroTheWave(self) -> None:
        self.data.setToZero()

    def energy(self) -> np.ndarray:
        # dPhiByDT = self.data.v
        #
        # dPhiByDX = np.zeros(self.propSpace.n)
        # dPhiByDX[:-1] = dPhiByDX[1:] - dPhiByDX[:-1]
        # dPhiByDX[-1] = dPhiByDX[0] - dPhiByDX[-1]
        # dPhiByDX /= self.propSpace.delta
        #
        # phi = self.data.x
        #
        # return (0.5 * dPhiByDT ** 2 + 0.5 * dPhiByDX ** 2 + 0.5 * self.propRSF.mSq * phi ** 2
        #         + 0.25 * self.propRSF.quarticTerm * phi ** 4
        #         - self.propRSF.isSometimesHeightOfPotMin()) * self.propSpace.delta

        # NOTE: Need to create copy as otherwise the arrays are updated by the integrators while doing calculations
        dPhiByDT = np.copy(self.data.v)
        phi = np.copy(self.data.x)

        dPhiByDX = np.diff(phi, 1, append=phi[0]) / self.propSpace.delta

        res = (0.5 * dPhiByDT ** 2 + 0.5 * dPhiByDX ** 2 + 0.5 * self.propRSF.mSq * phi ** 2
                + 0.25 * self.propRSF.quarticTerm * phi ** 4
                - self.propRSF.isSometimesHeightOfPotMin()) * self.propSpace.delta

        # term1 = 0.5 * dPhiByDT ** 2
        # term2 = 0.5 * dPhiByDX ** 2
        # term3 = 0.5 * self.propRSF.mSq * phi ** 2
        # term4 = 0.25 * self.propRSF.quarticTerm * phi ** 4
        # term5 = self.propRSF.isSometimesHeightOfPotMin()


        # print(f'max: {np.max(res)}')
        # print(f'min: {np.min(res)}')
        # if self.propRSF.mSq == -400 and np.min(res) < -0.0001:
        #     print(f'min: {np.min(res)}, max: {np.max(res)}, mSq: {self.propRSF.mSq}, q: {self.propRSF.quarticTerm}, heightOfPotMin: {self.propRSF.isSometimesHeightOfPotMin()}, delta: {self.propSpace.delta}')
            # with open("python.log", "a") as log_file:
            #     sys.stdout = log_file
            #     np.set_printoptions(threshold=sys.maxsize)
            #     print(phi)
            #     print(dPhiByDX)
            #     print(dPhiByDT)
            #     print(res)
            #
            #     print(term1)
            #     print(term2)
            #     print(term3)
            #     print(term4)
            #     print(term5)


        return res
