from enum import Enum, auto

import numpy as np
import pyqtgraph as pg
import pyqtgraph.parametertree
from PyQt5.QtWidgets import QCheckBox, QWidget, QVBoxLayout, QHBoxLayout
from pyqtgraph.parametertree import Parameter

from main_package.gui.Setting import Setting
from main_package.gui.widgets.IWavePanel import IWavePanel, IInfoBox
from main_package.simulation.DemoSetter import FrozenMultiset
from main_package.simulation.Interaction import LeapfroggableInteraction
from main_package.simulation.WaveThing import RSFWaveThing


class InfoBox(IInfoBox):
    """
    Implements a ParameterTree widget, which is basically a table with (Parameter, Value) pairs which can be modified.
    The ParameterTree widget is not resizable itself so it has to be wrapped inside another widget which is resizable.
    """

    def __init__(self, name: str, rsfWaveThing: RSFWaveThing):
        super().__init__()

        self._name = name
        self._rsfWaveThing = rsfWaveThing
        # Names of the fields in the table. Below, their names will be referenced by their position in this list.
        self._field_names = ['Name', 'Vev', 'mSq', 'quartic term']

        # Get the ParameterTree widget (which contains the Parameter object but we store it additionally as we will need it below)
        self._parameter, self._parameterTree = self._make_parameterTree()

        # Wrap the ParameterTree in a Layout
        self._layout = QVBoxLayout()
        self._layout.addWidget(self._parameterTree)

        # Wrap the Layout in a Widget
        self.setLayout(self._layout)
        self.setFixedWidth(self.MAXIMUM_WIDTH)
        self.setMaximumHeight(self.MAXIMUM_HEIGHT)

    def _make_parameterTree(self) -> tuple[Parameter, pg.parametertree.ParameterTree]:
        """Uses self._rsfWaveThing to make a ParameterTree widget."""

        # Get the Vacuum expectation value
        try:
            vev: str = np.format_float_positional(self._rsfWaveThing.get_propRSF().nonNegVev(), precision=2)
        except Exception as e:
            vev: str = str(e)

        # Set the fields in the ParameterTree
        params = [{'name': self._field_names[0], 'type': 'str', 'value': self._name, 'readonly': True},
                  {'name': self._field_names[1], 'type': 'str', 'value': vev, 'readonly': True},
                  {'name': self._field_names[2], 'type': 'float', 'value': self._rsfWaveThing.get_propRSF().get_mSq()},
                  {'name': self._field_names[3], 'type': 'float', 'value': self._rsfWaveThing.get_propRSF().get_quarticTerm()}
                 ]
        p = Parameter.create(name='params', type='group', children=params)

        # Connect the values that can be changed (mSq and quarticTerm) to changes in self._rsfWaveThing
        p.child(self._field_names[2]).sigValueChanged.connect(self._mSq_changed)
        p.child(self._field_names[3]).sigValueChanged.connect(self._quarticTerm_changed)

        # Create the ParameterTree widget
        parametertree = pg.parametertree.ParameterTree()
        parametertree.setParameters(p, showTop=False)

        # Return both the Parameter object and the ParameterTree. The Parameter object is needed because we need to
        # later access the individual fields in the table, which can be only done using the Parameter object.
        return p, parametertree

    def set_rsfWaveThing(self, reset_parameters: bool, wave: RSFWaveThing):
        """Called when the RSFWaveThing is replaced in the WavePanel. All fields in the table must be changed to the new
        values, and they must be re-connected to the new RSFWaveThing."""

        self._rsfWaveThing = wave

        # Get the new vev
        try:
            vev: str = np.format_float_positional(self._rsfWaveThing.get_propRSF().nonNegVev(), precision=2)
        except Exception as e:
            vev: str = str(e)

        # Set new default values (default values are the ones used when the 'undo' buttons in the ParameterTree are clicked)
        self._parameter.child(self._field_names[1]).setDefault(vev)
        self._parameter.child(self._field_names[2]).setDefault(self._rsfWaveThing.get_propRSF().get_mSq())
        self._parameter.child(self._field_names[3]).setDefault(self._rsfWaveThing.get_propRSF().get_quarticTerm())

        # Reconnect the modifiable fields to the new RSFWaveThing
        self._parameter.child(self._field_names[2]).sigValueChanged.connect(self._mSq_changed)
        self._parameter.child(self._field_names[3]).sigValueChanged.connect(self._quarticTerm_changed)

        if reset_parameters:
            # Change the values displayed in the ParameterTree.
            # This piece of code is only executed if reset_parameters is True, which occurs when the button 'Restart and reset parameters' is pressed.
            # If reset_parameters is False, the values that are currently displayed in the ParameterTree are kept (which
            # might not be the default value because the user has changed the value). The next piece of code then
            # synchronises the values displayed and the ones stored in the WaveThing.
            self._parameter.child(self._field_names[1]).setValue(vev)
            self._parameter.child(self._field_names[2]).setValue(self._rsfWaveThing.get_propRSF().get_mSq())
            self._parameter.child(self._field_names[3]).setValue(self._rsfWaveThing.get_propRSF().get_quarticTerm())

        # Synchronise the values in the ParameterTree and in the WaveThing.
        # This is done by emitting sigValueChanged (the signal emitted when the value is changed in the GUI by the user)
        # and passing to it the value that is CURRENTLY displayed in the ParameterTree.
        # This causes self._mSq_changed and self._quarticTerm_changed to be called, which ensure that the value displayed
        # in the ParameterTree is the same as the value stored in the WaveThing.
        for i in (2, 3): # Do this for i=2 (mSq) and i=3 (quarticTerm)
            parameter = self._parameter.child(self._field_names[i]) # Alias for readability
            parameter.sigValueChanged.emit(parameter, parameter.value())

    def _mSq_changed(self, param: pg.parametertree.parameterTypes.SimpleParameter, value: float) -> None:
        """
        Called when the value of mSq is changed in the ParameterTree.
        :param param: A Parameter object containing the value. Can be ignored for our purposes here, but we must accept
        it as an argument because the signal which calls this method requires this.
        :param value: The value entered by the user.
        """

        # Set mSq in the RSFWaveThing to the new value
        self._rsfWaveThing.get_propRSF().set_mSq(value)

        # Update the Vacuum expectation value in ParameterTree (which depends on mSq and quarticTerm)
        try:
            vev: str = np.format_float_positional(self._rsfWaveThing.get_propRSF().nonNegVev(), precision=2)
        except Exception as e:
            vev: str = str(e)
        self._parameter.child(self._field_names[1]).setValue(vev)

    def _quarticTerm_changed(self, param: pg.parametertree.parameterTypes.SimpleParameter, value: float) -> None:
        """
        Called when the value of quarticTerm is changed in the ParameterTree.
        :param param: A Parameter object containing the value. Can be ignored for our purposes here, but we must accept
        it as an argument because the signal which calls this method requires this.
        :param value: The value entered by the user.
        """

        # Set quarticTerm in the RSFWaveThing to the new value
        self._rsfWaveThing.get_propRSF().set_quarticTerm(value)

        # Update the Vacuum expectation value in ParameterTree (which depends on mSq and quarticTerm)
        try:
            vev: str = np.format_float_positional(self._rsfWaveThing.get_propRSF().nonNegVev(), precision=2)
        except Exception as e:
            vev: str = str(e)
        self._parameter.child(self._field_names[1]).setValue(vev)


class WavePanel(IWavePanel):

    def __init__(self, rsfWaveThing: RSFWaveThing, name: str, setting: Setting = None):
        super().__init__()

        self._rsfWaveThing = rsfWaveThing
        self.name = name
        self._current_setting = setting

        self._plotWidget: pg.PlotWidget = pg.PlotWidget()
        self._plotWidget.setBackground('w')
        # Lock the x-axis scaling + position
        self._plotWidget.setMouseEnabled(x=False)

        # Store PlotItem and ViewBox for convenience
        self._plotItem: pg.PlotItem = self._plotWidget.getPlotItem()
        self._viewBox: pg.ViewBox = self._plotItem.getViewBox()

        self._plot: pg.PlotDataItem = self._plotWidget.plot()
        self._plot.setPen(pg.mkPen(color=self._rsfWaveThing.get_plot_options().get_colour_rgba(), width=2))
        self._plotItem.getAxis('left').setLabel('Displacement', color=self._rsfWaveThing.get_plot_options().get_colour_rgba())

        # Hide x-axis
        self._plotItem.hideAxis('bottom')
        # Disconnect the autoRange button (in the bottowm left corner) from all slots it's connected to
        # This has the effect that the button doesn't do anything when clicked.
        # self._plotItem.autoBtn.disconnect()
        self._plotItem.autoBtn.disable()
        # self._plotItem.autoBtn.clicked.connect(lambda:print('clicked'))

        # The two Vev lines, placed at ±vev
        self._vev_lines: tuple[pg.PlotDataItem,...] = (pg.PlotDataItem(), pg.PlotDataItem())
        self._viewBox.addItem(self._vev_lines[0])
        self._viewBox.addItem(self._vev_lines[1])

        # Initialise energy plots
        self._energy_viewBox: pg.ViewBox
        self._energy_plot: pg.PlotDataItem
        self._energy_reference: pg.PlotDataItem
        self._energy_fill: pg.FillBetweenItem
        self._init_energy_plot()

        # Make a small button
        self._autoRange_checkBox: QCheckBox
        self._init_auto_checkbox()

        # Make an InfoBox
        self._infoBox = InfoBox(self.name, self._rsfWaveThing)

        # Arrange everything in a layout
        left_layout = QVBoxLayout()
        left_layout.addWidget(self._infoBox)
        left_layout.addWidget(self._autoRange_checkBox)
        self._layout = QHBoxLayout()
        self._layout.addLayout(left_layout)
        self._layout.addWidget(self._plotWidget)
        self.setLayout(self._layout)

        self._energy_set_first_time: bool = False
        self._energy_range_for_normal_mode: float = 0

        self._current_mode = None

    def _init_energy_plot(self):
        # Create second y axis
        self._energy_viewBox: pg.ViewBox = pg.ViewBox()
        # self._energy_viewBox.setMouseEnabled(y=False)
        self._plotItem.showAxis('right')
        self._plotItem.scene().addItem(self._energy_viewBox)
        self._plotItem.getAxis('right').linkToView(self._energy_viewBox)
        self._energy_viewBox.setXLink(self._plotWidget)
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
        # self._energy_plot.setBrush(pg.mkBrush(color=cyan))
        # self._energy_plot.setFillLevel(0)

        self._energy_reference: pg.PlotDataItem = pg.PlotDataItem()
        self._energy_reference.setPen(pg.mkPen(color=self.ENERGY_COLOUR))
        # Set initial data for PLOT_NORMAL mode because self._energy_fill will only be displayed if this is initialised before the creation of self._energy_fill
        self._energy_reference.setData(np.zeros(self._rsfWaveThing.propSpace.n))
        self._energy_viewBox.addItem(self._energy_reference)

        self._energy_fill: pg.FillBetweenItem = pg.FillBetweenItem(
            self._energy_reference,
            self._energy_plot,
            pg.mkBrush(color=self.ENERGY_COLOUR))
        self._energy_viewBox.addItem(self._energy_fill)

    def _init_auto_checkbox(self):
        self._autoRange_checkBox: QCheckBox = QCheckBox('Auto scale')
        self._autoRange_checkBox.setStyleSheet('QCheckBox::indicator { width: 10px; height: 10px;}')
        self._autoRange_checkBox.setChecked(False)

        # If the checkbox is unchecked, set the range to -11 to 11 (and let the user have manual control over the range after that)
        def conditional_set_range():
            if not self._autoRange_checkBox.isChecked():
                self.set_y_range(
                    self._rsfWaveThing.get_plot_options().get_plot_range()[0],
                    self._rsfWaveThing.get_plot_options().get_plot_range()[1]
                )
                self._plotWidget.setMouseEnabled(y=True)

        self._autoRange_checkBox.clicked.connect(conditional_set_range)

    def get_name(self) -> str:
        return self.name

    def set_waveThings(self, reset_parameters: bool, wave: RSFWaveThing = None,
                       waves: dict[str, RSFWaveThing] = None,
                       interactions: dict[FrozenMultiset, LeapfroggableInteraction] = None) -> None:
        assert (wave is not None) and (waves is None) and (interactions is None)

        self._rsfWaveThing = wave
        self._infoBox.set_rsfWaveThing(reset_parameters, wave)

        # Get new initial y range from PlotOptions
        self.set_y_range(
            self._rsfWaveThing.get_plot_options().get_plot_range()[0],
            self._rsfWaveThing.get_plot_options().get_plot_range()[1]
        )

        # Get new initial energy y range from PlotOptions
        self.set_energy_y_range(*wave.get_plot_options().get_energy_plot_range())

    def set_y_range(self, min: float, max: float) -> None:
        self._viewBox.setYRange(min, max)

    def set_x_range(self, min: float, max: float) -> None:
        self._viewBox.setXRange(min, max)

    def set_energy_y_range(self, min: float, max: float):
        self._energy_viewBox.setYRange(min, max)

    def set_energy_x_range(self, min: float, max: float):
        self._energy_viewBox.setXRange(min, max)

    def updatePlot(self, plot_3D: bool, setting: Setting) -> None:
        # Convert bool to variable of type Mode
        new_plot_mode = self.Mode.PLOT_3D if plot_3D else self.Mode.PLOT_NORMAL

        # Store this here because it's used multiple times below
        energy_data: np.ndarray = self._rsfWaveThing.energy()

        if new_plot_mode == self.Mode.PLOT_NORMAL:
            # Set data
            self._plot.setData(self._rsfWaveThing.getX())
            self._energy_plot.setData(energy_data)
            self._energy_reference.setData(np.zeros(self._rsfWaveThing.propSpace.n))
            # Set Data for Vev lines
            try:
                vev: float = self._rsfWaveThing.get_propRSF().nonNegVev()
                self._vev_lines[0].setData(np.ones(self._rsfWaveThing.get_N()) * vev)
                self._vev_lines[1].setData(np.ones(self._rsfWaveThing.get_N()) * (-vev))
                for vev_line in self._vev_lines:
                    vev_line.setPen(pg.mkPen(color=(252, 182, 3, 200), width=1.5))
            except Exception as e:
                if str(e) == 'FLAT':
                    for vev_line in self._vev_lines:
                        vev_line.setData(np.zeros(self._rsfWaveThing.get_N()))
                        vev_line.setPen(pg.mkPen(color=(200,200,200,200), width=1.5))

            # If the previous mode was 3D mode...
            if (self._current_mode is None) or (new_plot_mode != self._current_mode):
                self._viewBox.clear()
                self._energy_viewBox.clear()

                self._plotItem.showAxis('left')
                self._plotItem.showAxis('right')

                self._viewBox.addItem(self._plot)
                self._viewBox.addItem(self._vev_lines[0])
                self._viewBox.addItem(self._vev_lines[1])
                self._energy_viewBox.addItem(self._energy_plot)
                self._energy_viewBox.addItem(self._energy_reference)
                self._energy_viewBox.addItem(self._energy_fill)

                self.set_y_range(
                    self._rsfWaveThing.get_plot_options().get_plot_range()[0],
                    self._rsfWaveThing.get_plot_options().get_plot_range()[1]
                )
                self.set_energy_y_range(*self._rsfWaveThing.get_plot_options().get_energy_plot_range())

            ## Store energy plot range from the first time the energy data is set in normal mode
            # self._current_setting is None occurs when wavepanels are created for the first time
            if (self._current_setting is None) or (setting != self._current_setting):
                # # The initial maximum height of the energy plot will occupy half of the plotting area
                # self._energy_range_for_normal_mode = 2 * np.max(energy_data)
                # self.set_energy_y_range(0, self._energy_range_for_normal_mode)
                # self._current_setting = setting
                self.set_energy_y_range(*self._rsfWaveThing.get_plot_options().get_energy_plot_range())

            # Autoscale algorithm: Rescale axes if the data is out of bounds or too small to be seen
            if self._autoRange_checkBox.isChecked():
                self._plotWidget.setMouseEnabled(y=False)

                y_min, y_max = self._plot.dataBounds(1)
                axis_min, axis_max = self._plotItem.getAxis('left').range
                if (y_max - y_min) < 0.2 * (axis_max - axis_min):
                    self.fit_axes_to_data(new_plot_mode)
                elif (y_max - y_min) > 1.0 * (axis_max - axis_min) or (y_max > axis_max or y_min < axis_min):
                    self.fit_axes_to_data(new_plot_mode)

        else: # 3D
            # If the previous mode was PLOT_NORMAL...
            if new_plot_mode != self._current_mode:
                self._viewBox.clear()
                self._energy_viewBox.clear()

                self._plotItem.hideAxis('left')
                self._plotItem.hideAxis('right')

                self._viewBox.addItem(self._plot)
                self._viewBox.addItem(self._vev_lines[0])
                self._viewBox.addItem(self._vev_lines[1])
                self._viewBox.addItem(self._energy_plot)
                self._viewBox.addItem(self._energy_reference)
                self._viewBox.addItem(self._energy_fill)

                self._viewBox.enableAutoRange()

            # self._current_setting is None occurs when wavepanels are created for the first time
            if (self._current_setting is None) or (setting != self._current_setting):
                self._viewBox.enableAutoRange()

            n: int = self._rsfWaveThing.get_N()
            self._plot.setData(*self._get_3D_data(self._rsfWaveThing.getX(), n))
            # fudge factor of 16 to make the energy plots smaller... without that factor they appear too large
            self._energy_plot.setData(*self._get_3D_data(energy_data/16., n))
            self._energy_reference.setData(*self._get_3D_data(np.zeros(self._rsfWaveThing.propSpace.n), n))

            # Set Data for Vev lines
            try:
                vev: float = self._rsfWaveThing.get_propRSF().nonNegVev()
                self._vev_lines[0].setData(*self._get_3D_data(np.ones(self._rsfWaveThing.get_N()) * vev, n))
                self._vev_lines[1].setData(*self._get_3D_data(np.ones(self._rsfWaveThing.get_N()) * (-vev), n))
                for vev_line in self._vev_lines:
                    vev_line.setPen(pg.mkPen(color=(252, 182, 3, 200), width=1.5))
            except Exception as e:
                if str(e) == 'FLAT':
                    for vev_line in self._vev_lines:
                        vev_line.setData(*self._get_3D_data(np.zeros(n), n))
                        vev_line.setPen(pg.mkPen(color=(200, 200, 200, 200), width=1.5))

        # Update the current mode
        self._current_mode = new_plot_mode
        # Update the current setting
        self._current_setting = setting

    # Rescales plot axes to the current min and max values of the data
    def fit_axes_to_data(self, mode: IWavePanel.Mode):
        x_min, x_max = self._plot.dataBounds(0)
        y_min, y_max = self._plot.dataBounds(1)
        # Set range to (-11,11) if there is no data yet (then dataBounds returns None)
        y_min, y_max = (self._rsfWaveThing.get_plot_options().get_plot_range()[0], self._rsfWaveThing.get_plot_options().get_plot_range()[1]) if y_min is None or y_max is None else (y_min, y_max)

        if mode is self.Mode.PLOT_NORMAL:
            self.set_y_range(y_min, y_max)
            self.set_x_range(0, self._rsfWaveThing.propSpace.n)
        else: # 3D mode
            self.set_x_range(x_min - 2, x_max + 2)
            self.set_y_range(y_min - 10, y_max + 10)
            # self.set_y_range(-20, 20)
            # self.set_x_range(-9,9)



    # def updatePlot2(self, plot_3D: bool, setting: Setting) -> None:
    #     # Convert bool to variable of type Mode
    #     new_plot_mode = self.Mode.PLOT_3D if plot_3D else self.Mode.PLOT_NORMAL
    #
    #     # Store this here because it's used multiple times below
    #     energy_data: np.ndarray = self._rsfWaveThing.energy()
    #
    #     # Set the data
    #     if new_plot_mode is self.Mode.PLOT_NORMAL:
    #         # Curiously the GUI lags significantly when named arguments
    #         # are used instead: self._plot.setData(y=self._rsfWaveThing.data.x)
    #         self._plot.setData(self._rsfWaveThing.getX())
    #         self._energy_plot.setData(energy_data)
    #     else: # 3D mode
    #         self._plot.setData(*self._get_3D_data(self._rsfWaveThing.getX()))
    #         self._energy_plot.setData(*self._get_3D_data(energy_data))
    #
    #     # Store energy plot range from the first time the energy data is set in normal mode
    #     if setting is not self._current_setting:
    #         self._energy_range_for_normal_mode = 4 * np.max(energy_data)
    #         self._energy_set_first_time = True
    #         self.set_energy_y_range(0, self._energy_range_for_normal_mode)
    #         self._current_setting = setting
    #
    #
    #     # Autoscale algorithm: Rescale axes if the data is out of bounds or too small to be seen
    #     if self._autoRange_checkBox.isChecked():
    #         self._plotWidget.setMouseEnabled(y=False)
    #
    #         y_min, y_max = self._plot.dataBounds(1)
    #         axis_min, axis_max = self._plotItem.getAxis('left').range
    #         if (y_max - y_min) < 0.2 * (axis_max - axis_min):
    #             self.set_axes_scaling(new_plot_mode)
    #         elif (y_max - y_min) > 1.0 * (axis_max - axis_min) or (y_max > axis_max or y_min < axis_min):
    #             self.set_axes_scaling(new_plot_mode)
    #
    #     if new_plot_mode is not self._current_mode:
    #         # Only rescale the axes if we're switching modes
    #         self.set_axes_scaling(new_plot_mode)
    #         # Set energy reference fill level
    #         if new_plot_mode is self.Mode.PLOT_NORMAL:
    #             self._energy_reference.setData(np.zeros(self._rsfWaveThing.propSpace.n))
    #         else: # 3D mode
    #             self._energy_reference.setData(*self._get_3D_data(np.zeros(self._rsfWaveThing.propSpace.n)))
    #
    #     # Update the current mode
    #     self._current_mode = new_plot_mode

# class CustomDockArea(pg.dockarea.DockArea):
#     def __init__(self):
#         super().__init__()
#         self.state = None
#
#     def done(self):
#         self.state = super().saveState()
#
#     def restoreState(self, **kwargs):
#         super().restoreState(self.state)
