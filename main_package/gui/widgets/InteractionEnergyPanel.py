from enum import Enum, auto

import numpy as np
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout
import pyqtgraph as pg
from pyqtgraph.parametertree import Parameter

from main_package.gui.Setting import Setting
from main_package.gui.widgets.IWavePanel import IWavePanel, IInfoBox
from main_package.simulation.DemoSetter import FrozenMultiset
from main_package.simulation.Interaction import LeapfroggableInteraction
from main_package.simulation.InteractionEnergyDensity import InteractionEnergyDensity
from main_package.simulation.WaveThing import RSFWaveThing

class InteractionInfoBox(IInfoBox):
    def __init__(self):
        super().__init__()

        self._name = 'Interaction Energy'

        self._parameter, self._parameterTree = self._make_parameterTree()

        self._layout = QVBoxLayout()
        self._layout.addWidget(self._parameterTree)

        self.setLayout(self._layout)
        self.setFixedWidth(self.MAXIMUM_WIDTH)
        self.setMaximumHeight(self.MAXIMUM_HEIGHT)


    def _make_parameterTree(self):
        params = [
            {'name': 'Name', 'type': 'str', 'value': self._name, 'readonly': True}
        ]

        p = Parameter.create(name='params', type='group', children=params)

        # Create the ParameterTree widget
        parametertree = pg.parametertree.ParameterTree()
        parametertree.setParameters(p, showTop=False)

        return p, parametertree


class InteractionEnergyPanel(IWavePanel):
    def __init__(self, waves: dict[str, RSFWaveThing], interactions: dict[FrozenMultiset, LeapfroggableInteraction], setting: Setting = None):
        super().__init__()

        self._waves = waves
        self._interactions = interactions
        self._interactionEnergyDensity = InteractionEnergyDensity(waves, interactions)
        self._current_setting = setting

        self._plotWidget: pg.PlotWidget = pg.PlotWidget()
        self._plotWidget.setBackground('w')
        # Lock the x-axis scaling + position
        self._plotWidget.setMouseEnabled(x=False)

        # Store PlotItem and ViewBox for convenience
        self._plotItem: pg.PlotItem = self._plotWidget.getPlotItem()
        self._plotItem.hideAxis('bottom')
        self._viewBox: pg.ViewBox = self._plotItem.getViewBox()

        # Disable auto scale button in the bottom left corner
        self._plotItem.autoBtn.disable()

        # Make left axis transparent so that it still takes the same space as in the other WavePanels
        transparent = (255,255,255,0)
        self._plotWidget.getPlotItem().getAxis('left').setLabel('Displacement', color=transparent)
        self._plotWidget.getPlotItem().getAxis('left').setTextPen(pg.mkPen(color=transparent))

        # Create second y axis
        self._energy_viewBox: pg.ViewBox = pg.ViewBox()
        self._energy_viewBox.setMouseEnabled(y=False)
        self._plotItem.showAxis('right')
        self._plotItem.scene().addItem(self._energy_viewBox)
        self._plotItem.getAxis('right').linkToView(self._energy_viewBox)
        self._energy_viewBox.setXLink(self._plotWidget)
        self._energy_viewBox.setYLink(self._plotWidget)
        self._plotItem.getAxis('right').setLabel('Energy', color=self.ENERGY_COLOUR)
        def update_views():
            """ Method to synchronise the left and right y axis (Boilerplate code) """
            self._energy_viewBox.setGeometry(self._viewBox.sceneBoundingRect())
            self._energy_viewBox.linkedViewChanged(self._viewBox, self._energy_viewBox.XAxis)
        self._viewBox.sigResized.connect(update_views)

        # Add plot elements in the order in which they should displayed. We want the wave plot to be on top of the energy plot, so we add the wave plot last.
        self._energy_plot: pg.PlotDataItem = pg.PlotDataItem()
        self._energy_plot.setPen(pg.mkPen(color=self.ENERGY_COLOUR))
        self._energy_viewBox.addItem(self._energy_plot)

        self._energy_reference: pg.PlotDataItem = pg.PlotDataItem()
        self._energy_reference.setPen(pg.mkPen(color=self.ENERGY_COLOUR))
        # Set initial data for PLOT_NORMAL mode because self._energy_fill will only be displayed if this is initialised before the creation of self._energy_fill
        self._energy_reference.setData(np.zeros(self._interactionEnergyDensity.get_N()))
        self._energy_viewBox.addItem(self._energy_reference)

        self._energy_fill: pg.FillBetweenItem = pg.FillBetweenItem(
            self._energy_reference,
            self._energy_plot,
            pg.mkBrush(color=self.ENERGY_COLOUR))
        self._energy_viewBox.addItem(self._energy_fill)


        layout = QHBoxLayout()
        layout.addWidget(InteractionInfoBox())
        layout.addWidget(self._plotWidget)
        self.setLayout(layout)

        self._energy_range_for_normal_mode: float = 0
        self._current_mode = None


    def get_name(self) -> str:
        return 'Interaction Energy'

    def set_waveThings(self, reset_parameters: bool, wave: RSFWaveThing = None,
                       waves: dict[str, RSFWaveThing] = None,
                       interactions: dict[FrozenMultiset, LeapfroggableInteraction] = None) -> None:
        assert (wave is None) and (waves is not None) and (interactions is not None)

        self._waves = waves
        self._interactions = interactions
        # Also update the energy density calculator
        self._interactionEnergyDensity = InteractionEnergyDensity(waves, interactions)

    def updatePlot(self, plot_3D: bool, setting: Setting) -> None:
        # Convert bool to variable of type Mode
        new_plot_mode = self.Mode.PLOT_3D if plot_3D else self.Mode.PLOT_NORMAL

        # Store this here because it's used multiple times below
        energy_data: np.ndarray = self._interactionEnergyDensity.energy()

        if new_plot_mode == self.Mode.PLOT_NORMAL:
            self._energy_plot.setData(energy_data)
            self._energy_reference.setData(np.zeros(self._interactionEnergyDensity.get_N()))

            self._viewBox.setMouseEnabled(y=True)
            self._energy_viewBox.setMouseEnabled(y=False)

            # If the previous mode was 3D mode...
            if (self._current_mode is None) or (new_plot_mode != self._current_mode):
                self._viewBox.clear()
                self._energy_viewBox.clear()

                self._energy_viewBox.addItem(self._energy_plot)
                self._energy_viewBox.addItem(self._energy_reference)
                self._energy_viewBox.addItem(self._energy_fill)

                self._energy_viewBox.setXRange(0, self._interactionEnergyDensity.get_N())
                self._energy_viewBox.setYRange(0, self._energy_range_for_normal_mode)

                self._plotItem.showAxis('left')
                self._plotItem.showAxis('right')

            # Store energy plot range from the first time the energy data is set in normal mode
            # self._current_setting is None occurs when wavepanels are created for the first time
            if (self._current_setting is None) or (setting != self._current_setting):
                self._energy_range_for_normal_mode = 2 * np.max(energy_data)
                self._energy_viewBox.setYRange(0, self._energy_range_for_normal_mode)
                self._current_setting = setting

        else: # 3D
            # If the previous mode was PLOT_NORMAL...
            if new_plot_mode != self._current_mode:
                self._viewBox.clear()
                self._energy_viewBox.clear()

                self._viewBox.addItem(self._energy_plot)
                self._viewBox.addItem(self._energy_reference)
                self._viewBox.addItem(self._energy_fill)

                self._plotItem.hideAxis('left')
                self._plotItem.hideAxis('right')

                self._viewBox.enableAutoRange()

            self._viewBox.setMouseEnabled(y=True)
            self._energy_viewBox.setMouseEnabled(y=True)

            n: int = self._interactionEnergyDensity.get_N()
            self._energy_plot.setData(*self._get_3D_data(energy_data, n))
            self._energy_reference.setData(*self._get_3D_data(np.zeros(n), n))

        # Update the current mode
        self._current_mode = new_plot_mode
