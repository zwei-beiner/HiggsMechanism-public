import numpy as np

from main_package.simulation.DemoSetter import FrozenMultiset
from main_package.simulation.EnergyDensity import EnergyDensity
from main_package.simulation.Interaction import LeapfroggableInteraction
from main_package.simulation.WaveThing import RSFWaveThing


class InteractionEnergyDensity(EnergyDensity):
    def __init__(self, waves: dict[str, RSFWaveThing], interactions: dict[FrozenMultiset, LeapfroggableInteraction]):
        self._waves = waves
        self._interactions = interactions

        # Get n from the first waveThing in the dict
        self._n: int = self._waves[list(self._waves.keys())[0]].get_N()

    def get_N(self) -> int:
        return self._n

    def energy(self) -> np.ndarray:
        # Create local copy of the wave data
        wave_xs: dict[str, np.ndarray] = {name: np.copy(wave.getX()) for name, wave in self._waves.items()}

        multiset: FrozenMultiset
        interaction: LeapfroggableInteraction
        energyDensity: np.ndarray = np.zeros(self._n)
        # Loop through each interaction term, e.g. lambda_A_A_B_B
        for multiset, interaction in self._interactions.items():
            items: list[str] = multiset.get_items()
            term: np.ndarray = np.ones(self._n) * interaction.get_parameter()
            # In each interaction term, multiply the waves together which appear in that term
            for name in items:
                term *= wave_xs[name]

            # Add term to the energy density
            energyDensity += term

        return energyDensity