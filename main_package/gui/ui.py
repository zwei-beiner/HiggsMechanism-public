import functools

from PyQt5 import QtWidgets, QtCore, QtGui

from main_package.gui import GlobalStyle
from main_package.gui.Setting import Setting
from main_package.gui.Simulation import Simulation
from main_package.gui.widgets.Clock import Clock
from main_package.gui.widgets.ErrorDialog import ErrorDialog
from main_package.gui.widgets.InteractionEnergyPanel import InteractionEnergyPanel
from main_package.gui.widgets.WavePanel import WavePanel
from main_package.gui.widgets.CustomDock import CustomDock, CustomDockArea
from main_package.gui.widgets.InteractionCheckBox import InteractionCheckBox, InteractionCheckBoxLayout, CustomCheckBox
from main_package.simulation.DemoSetter import FrozenMultiset
from main_package.simulation.Interaction import LeapfroggableInteraction
from main_package.simulation.WaveThing import RSFWaveThing


class FormWidget(QtWidgets.QWidget):
    def __init__(self, layout: QtWidgets.QLayout):
        super().__init__()
        self.setLayout(layout)


class UI(QtWidgets.QMainWindow):
    """
    Handles user interface.
    Inherits from QWidget to be displayed as a window.
    Inherits from QObject to run the repainter loop.
    """

    def __init__(self, simulation: Simulation, fps: int = 30):
        super().__init__()

        self._fps = fps

        self._simulation = simulation
        self._simulation.build_finished.connect(self.remake_checkboxes_and_wavepanels)

        self.interactions: list[InteractionCheckBox] = []
        self.interactions_layout: InteractionCheckBoxLayout = InteractionCheckBoxLayout()

        self.configNames = [Setting.lookup(setting) for setting in Setting]
        self.currentConfigName = self.configNames[0]

        self._dockarea = CustomDockArea(Setting.lookup(self.currentConfigName))

        # Only sets the wave packets to their initial displacement at t=0. The parameters are not reset.
        self._restartButton = QtWidgets.QPushButton('Restart waves')
        # Sets all parameters to their initial values
        self._resetButton = QtWidgets.QPushButton('Restart and reset parameters')

        self._plot_3D: bool = False
        self._checkBox_3D: CustomCheckBox = CustomCheckBox('3D', add_stretch_at_the_end=False)

        self._box_dt = InteractionCheckBox.create_for_dt(True, simulation.get_dt(), b'\xce\x94'.decode() + 't')

        self._comboBox = QtWidgets.QComboBox()
        self._layout = QtWidgets.QHBoxLayout()

        self._clock: Clock = Clock(self._simulation)

        # Error Dialog is shown when overflow occurs because dt is too large
        self._error_dialog: ErrorDialog = ErrorDialog()

        self.makeGUI()

    def makeGUI(self):

        self._comboBox.addItems(self.configNames)
        # Pass in reset_parameters=True because we want to reset all simulation parameters when the setting is changed.
        self._comboBox.currentTextChanged.connect(functools.partial(self._simulation.restart, True))
        self._comboBox.currentTextChanged.connect(lambda s: self.__setattr__('currentConfigName', s))
        # self._comboBox.currentTextChanged.connect(lambda: self.__setattr__('_remake_checkboxes_and_wavepanels', True))
        self._comboBox.currentTextChanged.connect(
            lambda: self._simulation.set_setting(Setting.lookup(self.currentConfigName)))

        # self._restartButton.clicked.connect(lambda: self.__setattr__('_remake_checkboxes_and_wavepanels', True))

        self._checkBox_3D.setChecked(self._plot_3D)
        self._checkBox_3D.stateChanged.connect(lambda: self.__setattr__('_plot_3D', self._checkBox_3D.isChecked()))

        self._box_dt.signal_value_changed.connect(self._simulation.set_dt)

        # self._restartButton.setMaximumHeight(self._box_dt.height())
        self._restartButton.setStyleSheet(f'background-color: rgba{GlobalStyle.YELLOW};')
        # We don't want to set all simulation parameters to their default values when the restart button is clicked.
        self._restartButton.clicked.connect(functools.partial(self._simulation.restart, False))
        self._restartButton.setToolTip('Set waves to their initial displacement.')
        self._restartButton.setFixedHeight(30)

        self._resetButton.setStyleSheet(f'background-color: rgba{GlobalStyle.YELLOW};')
        # We want to set all simulation parameters to their default values when the restart button is clicked.
        self._resetButton.clicked.connect(functools.partial(self._simulation.restart, True))
        self._resetButton.setToolTip('Set waves to their initial displacement and set all parameters to their default values.')
        self._resetButton.setFixedHeight(30)

        # Pause the simulation and show an Error Dialog when a RuntimeWarning occurs in the simulation loop.
        def display_error_message(message: str):
            # Unclick the dt checkbox. This stops the simulation loop from running.
            self._box_dt._cb._checkBox.setChecked(False)
            # Show the Error Dialog.
            self._error_dialog.show_with_message(message)
        self._simulation.sig_RuntimeWarning_occurred.connect(display_error_message)

        upper_top_left_layout = QtWidgets.QHBoxLayout()
        upper_top_left_layout.addWidget(self._comboBox)
        upper_top_left_layout.addWidget(self._checkBox_3D)
        top_left_layout = QtWidgets.QVBoxLayout()
        top_left_layout.addWidget(self._dockarea.get_restore_button())
        top_left_layout.addLayout(upper_top_left_layout)

        top_layout = QtWidgets.QHBoxLayout()
        top_layout.addSpacing(45)
        top_layout.addLayout(top_left_layout)
        # Add a label as a spacing which is one quarter of the window size. self.resizeEvent() dynamicallt adjusts the width of the label.
        self._label_for_spacing = QtWidgets.QLabel()
        self._label_for_spacing.setFixedWidth(self.width() // 4)
        top_layout.addWidget(self._label_for_spacing)
        top_layout.addWidget(self._clock)
        top_layout.addStretch()

        # Build checkboxes and wavepanels to make them appear simultaneously with the other GUI elements (even though
        # they haven't been initialised yet by self._simulation)
        # As a result, self._simulation.getParameterCheckBoxes() and self._simulation.getWaveThings() is called twice
        # when the application is started. (The other method call is due to self.restart() being called at the beginning
        # of self._simulation.run())
        self.remake_checkboxes_and_wavepanels(False)

        left_layout = QtWidgets.QVBoxLayout()
        left_layout.addLayout(top_layout)
        left_layout.addWidget(self._dockarea)

        right_layout = QtWidgets.QVBoxLayout()
        right_layout.addSpacing(10)
        right_layout.addWidget(self._restartButton)
        right_layout.addWidget(self._resetButton)
        right_layout.addWidget(self._box_dt)
        right_layout.addLayout(self.interactions_layout)

        # bottom_layout = QtWidgets.QHBoxLayout()
        # bottom_layout.addWidget(self._dockarea)
        # bottom_layout.addLayout(self.interactions_layout)

        self._layout.addLayout(left_layout)
        self._layout.addLayout(right_layout)
        self.setCentralWidget(FormWidget(self._layout))

    def make_wave_panels(self, reset_parameters: bool) -> None:
        self._dockarea.make_wave_panels(reset_parameters, Setting.lookup(self.currentConfigName), *self._simulation.get_waveThings_and_parameters())

    def re_make_checkboxes(self, reset_parameters: bool) -> None:
        """
        Creates a list of new interaction checkboxes (on the right of the GUI) and adds them to the layout. The previously
        existing interaction checkboxes are deleted.
        :param reset_parameters: True if the parameter (lambda) for each checkbox should be reset to its default value.
        """

        # Keep list of existing checkboxes
        existing_names: dict[str, InteractionCheckBox] = {interaction.get_identifier(): interaction for interaction in
                                                          self.interactions}

        # Create list of new checkboxes
        new_interactions: dict[FrozenMultiset, LeapfroggableInteraction] = self._simulation.getParameterCheckBoxes()

        cbs: list[InteractionCheckBox] = []
        for frozenmultiset, interaction in new_interactions.items():
            # Alias for readability
            name = frozenmultiset.get_name()

            # Copy existing values, if present.
            # We look if the interaction checkbox already exists. If it does, the checked state and the currently
            # displayed parameter are adopted. This has the effect that, if the user changes the value displayed, these
            # changed will be carried over to the newly created checkboxes.
            cb = InteractionCheckBox.create_for_interaction_term(
                frozenmultiset.get_on() if not (name in existing_names) else existing_names[name].getInitialOnState(),
                interaction.get_parameter() if not (name in existing_names) else existing_names[name].get_current_value(),
                frozenmultiset,
                interaction)

            # If reset_parameters is True, the checkbox will display the default value.
            if reset_parameters:
                cb.reset_to_default_value()

            # When the value is changed in the checkbox, the lambda stored in LeapfroggableInteraction must be updated, too.
            # NOTE: Do NOT use cb.signal_value_changed.connect(lambda val: interaction.set_parameter(val))
            cb.signal_value_changed.connect(interaction.set_parameter)

            # Add the newly created checkbox to the list.
            cbs.append(cb)

        for cb in cbs:
            # Synchronise Interaction object and checkbox.
            # This is necessary because, sometimes after creation, the interaction parameter (lambda) in the LeapfroggableInteraction object
            # does not match the value currently displayed. Hence, to ensure that they are synchronised, we synchronise them
            # here manually.
            cb.value_changed()

        # Delete all previous checkboxes
        for interaction in self.interactions:
            interaction.deleteLater()
        self.interactions.clear()

        # Use new checkboxes
        self.interactions = cbs
        for cb in cbs:
            self.interactions_layout.addWidget(cb)

    def remake_checkboxes_and_wavepanels(self, reset_parameters: bool) -> None:
        # These methods are always called together so for convenience group them in one wrapper method.
        self.make_wave_panels(reset_parameters)
        self.re_make_checkboxes(reset_parameters)
        if reset_parameters:
            self._box_dt.reset_to_default_value()

        self._simulation.unlock()

    def startSimulation(self):
        self._thread = QtCore.QThread()
        self._simulation.moveToThread(self._thread)
        self._thread.started.connect(self._simulation.run)
        self._thread.start()

    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        '''
        Overrides QtCore.MainWindow.closeEvent()

        IMPORTANT: This method is automatically called when the MainWindow is closed.
        It quits the thread running in the background.
        Without this method, the program crashes with a segfault.
        All three commands stop(), quit() and wait() are necessary.
        '''

        # Stop the infinite loop in Simulation.run()
        self._simulation.stop()
        # Quit the thread
        self._thread.quit()
        # Wait until the thread quits
        self._thread.wait()

    def startRepainterThread(self):
        # # Repainter and thread object have to be created as class members.
        # # Otherwise (i.e. leaving out the "self.") they go out of scope and the thread is terminated.
        # # Note: self._repainter and self._thread are members which are dynamically added in this method and
        # # are not declared in the constructor.
        # self._repainter = Repainter(self.waves)
        # # self._repainter.updateWavePanels.connect(self.updateWavePanels)
        #
        # self._thread = QtCore.QThread()
        # self._repainter.moveToThread(self._thread)
        # self._thread.started.connect(self._repainter.repainterLoop)
        #
        # self._thread.start()

        self._timer = QtCore.QTimer()
        self._timer.timeout.connect(self._updateWavePanels)
        def update_clock():
            self._clock.update_time(self._box_dt.get_current_value())
        self._timer.timeout.connect(update_clock)
        self._timer.start(1000 // self._fps)

    def _updateWavePanels(self):
        # if self._remake_checkboxes_and_wavepanels:
        #     # for wave in self.waves:
        #     #     wave.deleteLater()
        #     # for interaction in self.interactions:
        #     #     interaction.deleteLater()
        #     #     self.interactions.remove(interaction)
        #
        #     # Make only the wave panels
        #     self.make_wave_panels()
        #     self._remake_checkboxes_and_wavepanels = False

        # Plots call WaveThing.energy() method, which raises a warning when dt is too large because of a float overflow
        try:
            for waveDock in self._dockarea.get_waveDocks():
                waveDock.get_panel().updatePlot(self._plot_3D, Setting.lookup(self.currentConfigName))
        except RuntimeWarning as e:
            # Warning is already handled in Simulation.run()
            # We don't want to call the Error Dialog here because this method (_updateWavePanels()) is called at 30fps,
            # and therefore the Error Dialogs would pile up.
            # The Error Dialog is called in the Simulation.run() method because we can stop the simulation loop so that
            # the Error Dialog is created only once.
            pass

    def resizeEvent(self, a0: QtGui.QResizeEvent) -> None:
        """
        Overrides resizeEvent.
        If the MainWindow is resized, the spacing to the left of the Clock should dynamically adjust to be a quarter of the current window size.
        """
        super().resizeEvent(a0)
        self._label_for_spacing.setFixedWidth(self.width() // 4)