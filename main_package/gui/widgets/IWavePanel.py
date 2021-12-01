from enum import auto, IntEnum

import numpy as np
from PyQt5 import QtWidgets

from main_package.gui.Setting import Setting
from main_package.simulation.DemoSetter import FrozenMultiset
from main_package.simulation.Interaction import LeapfroggableInteraction
from main_package.simulation.WaveThing import RSFWaveThing


# Abstract base class (but can't make it inherit from ABC because QWidget doesn't allow multiple inheritance with its classes)
class IInfoBox(QtWidgets.QWidget):
    MAXIMUM_WIDTH: int = 240
    MAXIMUM_HEIGHT: int = 122


# Abstract base class (but can't make it inherit from ABC because QWidget doesn't allow multiple inheritance with its classes)
class IWavePanel(QtWidgets.QWidget):
    ENERGY_COLOUR = '#547fff'

    class Mode(IntEnum):
        PLOT_NORMAL = auto()
        PLOT_3D = auto()

    def get_name(self) -> str:
        raise NotImplementedError()

    def updatePlot(self, plot_3D: bool, setting: Setting) -> None:
        raise NotImplementedError()

    def set_waveThings(self, reset_parameters: bool, wave: RSFWaveThing = None,
                       waves: dict[str, RSFWaveThing] = None,
                       interactions: dict[FrozenMultiset, LeapfroggableInteraction] = None) -> None:
        raise NotImplementedError()

    def _get_3D_data(self, data: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Maps data into a circle for 3D display.
        :param data: The y-coordinates of the data
        :param n: The number of points on the x-axis. I.e. the x-axis consists of n points
        :return: x- and y- coordinates of the data on a circle
        """

        # width = self._viewBox.width()
        # height = self._viewBox.height()
        # x_min, x_max = self._plot.dataBounds(0)
        # width = x_max - x_min
        # y_min, y_max = self._plot.dataBounds(1)
        # height = y_max - y_min

        # Circle
        angles: np.ndarray = np.linspace(0, 2 * np.pi, n)
        R = 5.
        x_base: np.ndarray = R * np.cos(angles)
        y_base: np.ndarray = R * np.sin(angles)

        # Offset y-values of the circle by y-values of the plot
        y: np.ndarray = data + y_base

        return x_base, y
