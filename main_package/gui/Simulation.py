from PyQt5 import QtCore, QtGui, QtWidgets

from main_package.gui.Setting import Setting

from main_package.simulation.DemoSetter import InitialConditionBuilder, InitialConditionDirector, FrozenMultiset

from main_package.simulation.InitialConditioners import InitialConditioner
from main_package.simulation.Integrator import LeapfrogIntegrator
from main_package.simulation.Interaction import LeapfroggableInteraction
from main_package.simulation.WaveThing import RSFWaveThing


class Simulation(QtCore.QObject):
    # Signal is emitted when the wave and interaction objects are created.
    # Sending a signal prevents race condition when checkboxes are re-built.
    # The bool carries the information if all simulation parameters should be reset (True) or not (False).
    build_finished = QtCore.pyqtSignal(bool)

    # Signal when a RuntimeWarning occurs during the simulation.
    # RuntimeWarnings occur when
    # - dt is too large (causing an overflow error)
    # - one of the lambda's is too large (causing an overflow error)
    # - one of the mSq's is negative (causing the sqrt to be taken of a negative number in the initialconditioner of a wave packet)
    sig_RuntimeWarning_occurred = QtCore.pyqtSignal(str)

    # def __init__(self, n=1001, L=60, mode: str = 'default'):
    def __init__(self, n: int = 1001, L: int = 60, dt: float = 0.01, setting: Setting = Setting.DEFAULT):
        super().__init__()

        self.n = n
        self.L = L
        self._setting = setting

        self._dt = dt
        self._t: float = 0.0
        # When set to True, only the wavepackets are set to their initial displacement at t=0
        self._restarting = False
        # When this flag is set to True in conjunction with self._restarting, all parameters (interaction lambdas, mSq and quarticTerms) are also set to their default values.
        self._reset_parameters = False
        self._stop = False

        # Is set to false when the restarting has to be paused because the GUI needs to be rebuilt.
        # This is necessary because rebuilding the GUI will synchronise the values displayed in the widget with the
        # simulation parameters (mSq, quarticTerm, lambda parameters) so that the values displayed (which the user might
        # have changed) are the ones used in the simulation.
        # Therefore, before initialising the wave packets, we have to wait until this synchronisation is finished.
        self._lock_the_restart_loop = False

        self.integrator: LeapfrogIntegrator = LeapfrogIntegrator()
        self.initialConditioners: dict[str, InitialConditioner] = {}
        # Uses internal names
        self._rsfWaveThings: dict[str, RSFWaveThing] = {}
        self._parameters: dict[FrozenMultiset, LeapfroggableInteraction] = {}

        self._build_from_setting()

    def get_dt(self) -> float:
        return self._dt

    def set_dt(self, dt: float) -> None:
        self._dt = dt

    def get_time(self) -> float:
        return self._t

    def set_setting(self, setting: Setting) -> None:
        self._setting = setting

    def stop(self):
        # Stops infinite loop in self.run()
        self._stop = True

    def unlock(self) -> None:
        self._lock_the_restart_loop = False

    def restart(self, reset_parameters: bool) -> None:
        # Pause integrator
        self._restarting = True
        # If self._reset_parameters is set to True, all parameter will also be set to their default values.
        self._reset_parameters = reset_parameters

    # Is called twice when the application is started
    def getParameterCheckBoxes(self) -> dict[FrozenMultiset, LeapfroggableInteraction]:
        return self._parameters

    # def getParameterCheckBoxes(self) -> tuple[list[str], list[bool], list[LeapfroggableInteraction]]:
    #     names: list[str] = []
    #     turn_ons: list[bool] = []
    #     for multiset in self._parameters.keys():
    #         names.append(multiset.get_name())
    #         turn_ons.append(multiset.get_on())
    #     return names, turn_ons, list(self._parameters.values())

    # Used ONLY for when the InteractionEnergyPanel is created
    def get_waveThings_and_parameters(self) -> tuple[dict[str, RSFWaveThing], dict[FrozenMultiset, LeapfroggableInteraction]]:
        return self._rsfWaveThings, self._parameters

        # res = []
        # # if self._setting is Setting.DEFAULT:
        # for multiset, interaction in self._parameters.items():
        #     label: str = multiset.get_name()
        #     val: float = interaction.get_parameter()
        #     initial_on_state: bool = multiset.get_on()
        #
        #     cb = InteractionCheckBox(initial_on_state, val, label)
        #     # NOTE: Do NOT use cb.signal_value_changed.connect(lambda val: interaction.set_parameter(val))
        #     cb.signal_value_changed.connect(interaction.set_parameter)
        #
        #     res.append(cb)
        #
        # # Synchronise Interaction objects and checkboxes
        # for cb in res:
        #     cb.value_changed()
        #
        # return res

    # Is called twice when the application is started
    def getWaveThings(self) -> tuple[list[str], list[RSFWaveThing]]:
        return list(self._rsfWaveThings.keys()), list(self._rsfWaveThings.values())


        # # print(self._rsfWaveThings)
        # res = []
        # # if self._setting is Setting.DEFAULT:
        # for name, rsfWaveThing in self._rsfWaveThings.items():
        #     res.append(WavePanel(rsfWaveThing, name))
        # return res

    # def changeName(self, old_key: str, new_key: str) -> None:
    #     self._waveNames[new_key] = self._waveNames.pop(old_key)

    def _setInitialConditions(self) -> None:
        # Set initial conditions
        for initialConditioner in self.initialConditioners.values():
            initialConditioner.setInitialCondition()

    # TODO: Speed this up with Cython?
    @QtCore.pyqtSlot()
    def run(self, infinite: bool = True, make_smooth_using_delay: bool = True, timeSteps: int = 500) -> None:
        # Initial restart. Causes the infinite loop to initialise wave packets and wait for 500ms.
        self.restart(False)

        # Simulation loop
        # NOTE: Any QtCore.QThread.sleep() statements should be inside this simulation loop.
        # sleep commands outside the simulation loop freeze the GUI (for unknown reasons).
        if infinite:
            while not self._stop:
                if not self._restarting:
                    # tic = time.perf_counter()
                    # if int(t / dt) % 900 == 0:
                    #     for name, interaction in self._parameters.items():
                    #         print(f'{name.get_name()} : {interaction.get_parameter()}')
                    #     print()

                    # If dt (or one of the lambda's) is too large, numpy will raise a RuntimeWarning because of an overflow. This warning is
                    # caught and a signal is emitted. In the class ui.py, the signal causes an ErrorDialog window to pop up.
                    try:
                        self.integrator.advance(self._dt)
                    except RuntimeWarning as e:
                        self.sig_RuntimeWarning_occurred.emit(str(e))

                    self._t += self._dt

                    if make_smooth_using_delay:
                        # 0.5ms delay
                        # QtCore.QThread.usleep(int(dt * 1e4))
                        # TODO: Automatic calculation of the optimal sleep time
                        QtCore.QThread.usleep(500)

                    # toc = time.perf_counter()
                    # with open('python_time_in_ms.log', 'a') as file:
                    #     file.write(str((toc - tic) * 1e3))
                    #     file.write('\n')
                elif self._restarting:
                    self._t = 0.0

                    # Get which wave things and interactions to display
                    self._build_from_setting()

                    # Initiate re-building of the GUI
                    self.build_finished.emit(self._reset_parameters)

                    # Wait for the GUI to call self.unlock(), which signals that the building of the wave panels and the
                    # interaction checkboxes has finished
                    self._lock_the_restart_loop = True
                    while self._lock_the_restart_loop:
                        # pass
                        QtWidgets.QApplication.processEvents()

                    try:
                        # Reset wave packets
                        self._setInitialConditions()
                    except RuntimeWarning as e:
                        # Catch RuntimWarning when sqrt is taken of a negative number (when one of the mSq's is negative)
                        self.sig_RuntimeWarning_occurred.emit(str(e))

                    # Wait for 500ms so that the initial wave packets can be seen by the user
                    QtCore.QThread.msleep(500)

                    # Resume integrator
                    self._restarting = False
        else:
            for i in range(timeSteps):
                if not self._restarting:
                    self.integrator.advance(self._dt)
                    self._t += self._dt

    def _build_from_setting(self) -> None:
        builder: InitialConditionBuilder = InitialConditionBuilder(self.n, self.L)
        director: InitialConditionDirector = InitialConditionDirector(builder)
        director.build_from_setting(self._setting)
        self.integrator, self.initialConditioners, self._rsfWaveThings, self._parameters = builder.get_result()


    # def _default(self) -> None:
    #     initialConditioners = self.initialConditioners
    #     integrator = self.integrator
    #
    #     rsfWaveThings = self._rsfWaveThings
    #
    #     # self._waveNames = {'workspace': 'A', 'heavy control': 'B', 'higgs': 'C', 'photon': 'D'}
    #     # self._parameters = {'A': (0, 0), 'B': (2.2 * 2.2, 0), 'C': (-2 * 2, 0.02 * 2 * 2), 'D': (0, 0),
    #     #                     'lambda_AC': 0.1}
    #     self._parameters = {'lambda_AC': 0.1}
    #
    #     ##################################################
    #
    #     n = self.n
    #     L = self.L
    #
    #     propSpace = SpaceProperties(n, L)
    #
    #     pCentral = 4
    #     pSpread = 1
    #     centralX = 5
    #     normScale = 6
    #
    #     ##################################################
    #
    #     rsfWaveThings['A'] = RSFWaveThing(propSpace, RSFProperties(mSq=0, quarticTerm=0))
    #
    #     initialConditioners.append(RSFInitialConditioner_PacketType(
    #         rsfWaveThings['A'],
    #         pCentral,
    #         pSpread,
    #         centralX,
    #         normScale))
    #     integrator.add(rsfWaveThings['A'])
    #
    #     ##################################################
    #
    #     rsfWaveThings['B'] = RSFWaveThing(propSpace, RSFProperties(mSq=2.2 * 2.2, quarticTerm=0))
    #
    #     initialConditioners.append(RSFInitialConditioner_PacketType(
    #         rsfWaveThings['B'],
    #         pCentral,
    #         pSpread,
    #         centralX,
    #         normScale))
    #     integrator.add(rsfWaveThings['B'])
    #
    #     ##################################################
    #
    #     narrower = 10
    #     rsfWaveThings['C'] = RSFWaveThing(propSpace, RSFProperties(mSq=-2 * 2 * (narrower ** 2), quarticTerm=0.02 * 2 * 2 * (narrower ** 2)))
    #     rsfAABBInteraction_AC = RSFAABBInteraction(rsfWaveThings['A'], rsfWaveThings['C'], 0.1)
    #
    #     initialConditioners.append(
    #         RSFInitialConditioner_HiggsTypeMatchedToRSFAABBInteraction(rsfWaveThings['C'], rsfWaveThings['A'],
    #                                                                    rsfAABBInteraction_AC))
    #
    #     integrator.add(rsfWaveThings['C'])
    #     integrator.add(rsfAABBInteraction_AC)
    #
    #     ##################################################
    #
    #     rsfWaveThings['D'] = RSFWaveThing(propSpace, RSFProperties(mSq=0, quarticTerm=0))
    #
    #     initialConditioners.append(
    #         RSFInitialConditioner_PacketType(rsfWaveThings['D'], pCentral, pSpread, centralX, normScale))
    #     integrator.add(rsfWaveThings['D'])
    #
    #     ##################################################
    #
