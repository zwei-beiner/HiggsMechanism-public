import sys
from abc import ABC, abstractmethod

import numpy as np
from PyQt5 import QtCore

from main_package.simulation.WaveThing import RSFWaveThing


# Abstract class
class LeapfroggableInteraction(ABC):
    @abstractmethod
    def get_parameter(self) -> float:
        pass

    @abstractmethod
    def set_parameter(self, lam: float) -> float:
        pass

    @abstractmethod
    def advanceVBasedOnX(self, dt: float) -> None:
        pass


class RSFAABBInteraction(LeapfroggableInteraction):
    '''
    Adds term of type
    1/2 * lambdaAABB * PhiA^2 * PhiB^2
    to the Lagrangian
    '''

    def __init__(self, A: RSFWaveThing, B: RSFWaveThing, lambdaAABB: float):
        self.A = A
        self.B = B
        self.lambdaAABB = lambdaAABB

    def get_parameter(self) -> float:
        return self.lambdaAABB

    def set_parameter(self, lambdaAABB: float) -> None:
        self.lambdaAABB = lambdaAABB

    def advanceVBasedOnX(self, dt: float) -> None:
        # print(self.lambdaAABB)
        self.A.incrementV(-dt * self.lambdaAABB * self.A.getX() * (self.B.getX() ** 2))
        self.B.incrementV(-dt * self.lambdaAABB * self.B.getX() * (self.A.getX() ** 2))


class RSFABBInteraction(LeapfroggableInteraction):
    '''
    Adds term of type
    lambdaABB * PhiA * PhiB^2
    to the Lagrangian
    '''

    def __init__(self, A: RSFWaveThing, B: RSFWaveThing, lambdaABB: float):
        self.A = A
        self.B = B
        self.lambdaABB = lambdaABB

    def get_parameter(self) -> float:
        return self.lambdaABB

    def set_parameter(self, lam: float) -> None:
        self.lambdaABB = lam

    def advanceVBasedOnX(self, dt: float) -> None:
        self.A.data.v -= dt * self.lambdaABB * self.B.data.x ** 2
        self.B.data.v -= dt * 2 * self.lambdaABB * self.A.data.x * self.B.data.x


class RSFAnBnInteraction(LeapfroggableInteraction):
    '''
    Adds term of type
    lam * PsiA^n_A * PsiB^n_B
    to the Lagrangian

    NOTE: Class doesn't use a factor in front of lam. Hence, to replace the class RSFAABBInteraction, a different parameter lam is required.
    '''

    def __init__(self, A: RSFWaveThing, n_A: int, B: RSFWaveThing, n_B: int, lam: float):
        '''
        Takes two waves and and their exponents in the interaction term and the interaction parameter.
        :param A: RSFWaveThing of wave A
        :param n_A: Exponent of wave A
        :param B: RSFWaveThing of wave B
        :param n_B: Exponent of wave B
        :param lam: Interaction parameter
        '''

        if n_A < 1 or n_B < 1:
            raise Exception('Interaction term not valid')

        self.A = A
        self.n_A = n_A

        self.B = B
        self.n_B = n_B

        self.lam = lam

    def get_parameter(self) -> float:
        return self.lam

    def set_parameter(self, lam: float) -> None:
        self.lam = lam

    def advanceVBasedOnX(self, dt: float) -> None:
        self.A.incrementV(-dt * self.lam * self.n_A * (self.A.getX() ** (self.n_A - 1)) * (self.B.getX() ** self.n_B))
        self.B.incrementV(-dt * self.lam * (self.A.getX() ** self.n_A) * self.n_B * (self.B.getX() ** (self.n_B - 1)))


class RSFInteraction_FullGenerality(LeapfroggableInteraction):
    def __init__(self, lam: float, **waves: tuple[int, RSFWaveThing]):
        self.lam = lam
        self.waves: dict[str, tuple[int, RSFWaveThing]] = waves

    def get_parameter(self) -> float:
        return self.lam

    def set_parameter(self, lam: float) -> None:
        self.lam = lam

    # TODO: Optimise this?
    def advanceVBasedOnX(self, dt: float) -> None:
        name: str
        n: int
        wave: RSFWaveThing

        name2: str
        n2: int
        wave2: RSFWaveThing

        for name, (n, wave) in self.waves.items():
            res: np.ndarray = self.lam * n * wave.getX() ** (n - 1)

            for name2, (n2, wave2) in self.waves.items():
                if name != name2:
                    res *= wave2.getX() ** n2

            wave.incrementV(-dt * res)

# class RSFInteraction_FullGenerality(LeapfroggableInteraction):
#     def __init__(self, lam: float, **waves: tuple[int, RSFWaveThing]):
#         self.lam: float = lam
#         self.waves: dict[str, tuple[int, RSFWaveThing]] = waves
#         # Number of waves
#         self.m: int = len(waves)
#         # Number of x values per wave
#         self.n = list(waves.values())[0][1].propSpace.n
#
#         # Keep list of RSFWaveThings
#         self.waveThings: list[RSFWaveThing] = [w[1] for w in waves.values()]
#         # Keep references to x arrays, 2D array
#         self.xs: np.ndarray = np.array([w[1].getX() for w in waves.values()])
#         # Prepare boolean mask
#         b = np.diag(np.ones(self.m)).astype(bool)
#         b = np.array([b] * self.n)
#         self.b = np.transpose(b)
#
#     def get_parameter(self) -> float:
#         return self.lam
#
#     def set_parameter(self, lam: float) -> None:
#         self.lam = lam
#
#     def advanceVBasedOnX(self, dt: float) -> None:
#         xs_stacked = np.array([self.xs] * self.m)
#         xs_stacked[self.b] = 1
#         res = np.prod(xs_stacked, axis=1)
#
#         # 1D array
#         e: np.ndarray = self.lam * np.fromiter((phi.getX() ** (n - 1) for n, phi in self.waves.values()), dtype=np.float)
#         # 2D array
#         E: np.ndarray = e * res
#
#         self.xs += -dt * E
