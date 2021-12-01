from abc import ABC, abstractmethod

import numpy as np


# Abstract class
class EnergyDensity(ABC):
    @abstractmethod
    def energy(self) -> np.ndarray:
        raise NotImplementedError()