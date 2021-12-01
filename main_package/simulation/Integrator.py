from functools import singledispatchmethod

from main_package.simulation.WaveThing import LeapfroggableWaveThing
from main_package.simulation.Interaction import LeapfroggableInteraction


class LeapfrogIntegrator:
    def __init__(self):
        self.waveThings: set[LeapfroggableWaveThing] = set()
        self.interactions: set[LeapfroggableInteraction] = set()

    def __str__(self):
        return "waveThings:" + str(self.waveThings) + "\n interactions:" + str(self.interactions)

    def _advanceVBasedOnX(self, dt: float) -> None:
        for waveThing in self.waveThings:
            waveThing.advanceVBasedOnX(dt)
        for interaction in self.interactions:
            interaction.advanceVBasedOnX(dt)

    def _advanceXBasedOnV(self, dt: float) -> None:
        for waveThing in self.waveThings:
            waveThing.advanceXBasedOnV(dt)

    def advance(self, dt: float) -> None:
        if dt > 0:
            self._advanceVBasedOnX(dt * 0.5)
            self._advanceXBasedOnV(dt * 0.5)
        elif dt < 0:
            self._advanceXBasedOnV(dt * 0.5)
            self._advanceVBasedOnX(dt * 0.5)

    def clear(self):
        self.waveThings.clear()
        self.interactions.clear()

    # functools.singledispatchmethod is Python's way of overloading methods based on argument type
    @singledispatchmethod
    def add(self, thing) -> None:
        raise TypeError("This type isn't supported: {}".format(type(thing)))

    @add.register
    def _(self, thing: LeapfroggableWaveThing) -> None:
        self.waveThings.add(thing)

    @add.register
    def _(self, thing: LeapfroggableInteraction) -> None:
        self.interactions.add(thing)