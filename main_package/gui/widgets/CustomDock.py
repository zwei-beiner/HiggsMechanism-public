import functools
import itertools
from typing import Iterator, Optional

import pyqtgraph as pg
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFontMetrics
from PyQt5.QtWidgets import QHBoxLayout
from pyqtgraph import dockarea

from main_package.gui import GlobalStyle
from main_package.gui.Setting import Setting
from main_package.gui.widgets.IWavePanel import IWavePanel
from main_package.gui.widgets.InteractionEnergyPanel import InteractionEnergyPanel
from main_package.gui.widgets.WavePanel import WavePanel
from main_package.simulation.DemoSetter import FrozenMultiset
from main_package.simulation.Interaction import LeapfroggableInteraction
from main_package.simulation.WaveThing import RSFWaveThing


class CustomDock(dockarea.Dock):
    """Dock which holds an IWavePanel object"""

    def __init__(self, panel: IWavePanel, name, *args, **kwargs):
        super().__init__(name, widget=panel, closable=True, *args, **kwargs)
        self._panel = panel

    def get_panel(self) -> IWavePanel:
        return self._panel

    # Overrides method in the superclass
    def close(self):
        """Remove this dock from the DockArea it lives inside."""
        self.setParent(None)
        # ----------------------------------
        # Code is copied from the superclass Dock but without the following line, which deletes the label. We want to keep the label.
        # QtGui.QLabel.close(self.label)
        # ----------------------------------
        self.label.setParent(None)
        self._container.apoptose()
        self._container = None
        self.sigClosed.emit(self)


class SmallButton(QtWidgets.QPushButton):
    def __init__(self, text: str):
        super().__init__(text)
        self.setFixedSize(100, 15)
        # Call overriden method to display elided text
        self.setText(text)

    def setText(self, a0: str) -> None:
        """Override to make text elided (=add '...' if it is too long to display)"""
        metrics: QFontMetrics = self.fontMetrics()
        elided_text: str = metrics.elidedText(a0, Qt.ElideRight, self.width())
        super().setText(elided_text)


class ClosedDockBar(QtWidgets.QWidget):
    sigRemoveButton = QtCore.pyqtSignal(CustomDock)

    def __init__(self):
        super().__init__()
        self._closed_waves: dict[str, CustomDock] = {}
        self._buttons: dict[str, SmallButton] = {}

        self._layout = QHBoxLayout()
        self._layout.addStretch()
        self.setLayout(self._layout)

    def add_dock(self, dock: CustomDock):
        text = dock.get_panel().get_name()
        # if text not in self._closed_waves: # don't allow duplicates
        self._closed_waves[text] = dock

        self._buttons[text] = SmallButton(text)
        self._buttons[text].clicked.connect(functools.partial(self._remove_button, text))
        # self._layout.addWidget(self._buttons[text])
        self._layout.insertWidget(self._layout.count() - 1, self._buttons[text])

    def _remove_button(self, text: str):
        button_to_delete = self._buttons.pop(text)
        button_to_delete.deleteLater()

        dock_to_reinsert = self._closed_waves.pop(text)
        self.sigRemoveButton.emit(dock_to_reinsert)

    def clear(self):
        # Delete all the buttons
        for button in self._buttons.values():
            button.deleteLater()
        self._buttons.clear()
        # Remove all stored docks from the dict but don't delete them
        self._closed_waves.clear()


class DockArea(pg.dockarea.DockArea):
    """Override the DockArea class to prevent the docks from collapsing when they are too small."""
    def makeContainer(self, typ):
        new = super().makeContainer(typ)
        new.setChildrenCollapsible(False)
        return new


class CustomDockArea(QtWidgets.QWidget):
    def __init__(self, currentConfig: Setting):
        super().__init__()

        self.currentConfig: Optional[Setting] = None

        self.waveDocks: list[CustomDock] = [] # stores ALL wavedocks (open and closed)
        self.wave_DockArea: DockArea = DockArea()
        # dock_layout is needed because below we will be switching the dock area with another one, and the functionality
        # of replacing a widget with another one is only given by a layout object. Hence, we need to wrap the dockarea
        # in a layout
        self.dock_layout: QtWidgets.QLayout = QtWidgets.QVBoxLayout()
        # self._dock_state: dict = None

        # Interaction energy dock is kept separate from the other wavedocks
        self._interactionEnergy_dock: CustomDock = None
        # self.wave_DockArea.addDock(self._interactionEnergy_dock)

        self._restore_layout_button = QtWidgets.QPushButton('Restore Layout')
        self._restore_layout_button.setToolTip('Reset the layout of the plots to the default order.')
        self._restore_layout_button.setStyleSheet(f'background-color: rgba{GlobalStyle.YELLOW};')
        def load():
            # 'ignore' flag ignores if docks are missing
            self.wave_DockArea.restoreState(self._make_saveState(), missing='ignore')
        self._restore_layout_button.clicked.connect(load)

        # self._closed_docks: list[CustomDock] = []
        self._closedDockBar: ClosedDockBar = ClosedDockBar()
        self._closedDockBar.sigRemoveButton.connect(self.add_to_dockArea)

        self.dock_layout.addWidget(self.wave_DockArea)
        self.dock_layout.addWidget(self._closedDockBar)
        self.setLayout(self.dock_layout)

    def get_restore_button(self) -> QtWidgets.QPushButton:
        return self._restore_layout_button

    def get_waveDocks(self) -> Iterator[CustomDock]:
        # Return generator which returns all the waveDocks, which are either in the dock or in self._closed_docks
        return itertools.chain(self.waveDocks, [self._interactionEnergy_dock])

    def add_to_dockArea(self, dock: CustomDock):
        self.wave_DockArea.addDock(dock)

    def make_wave_panels(self, reset_parameters: bool, currentConfig: Setting, wavethings: dict[str, RSFWaveThing], interactions: dict[FrozenMultiset, LeapfroggableInteraction]) -> None:
        # If self.waveDocks is empty, copy all of the new waves and return
        if len(self.waveDocks) == 0:
            for name, waveThing in wavethings.items():
                waveDock = CustomDock(WavePanel(waveThing, name), name)
                self.waveDocks.append(waveDock)
                self.wave_DockArea.addDock(waveDock)

            # Add the InteractionEnergyPanel
            self._interactionEnergy_dock = CustomDock(InteractionEnergyPanel(waves=wavethings, interactions=interactions), 'Interaction Energy')
            self.wave_DockArea.addDock(self._interactionEnergy_dock)

        else:
            # else the existing wavepanels have to be updated... (updating the existing wavePanels is better because
            # creating new ones causes the GUI to collapse and expand again which looks like it's glitching)

            # Store minimum and maximum indices of the lists
            minimum = min(len(self.waveDocks), len(wavethings))
            maximum = max(len(self.waveDocks), len(wavethings))

            # Copy new wavethings into existing wavepanels so they don't have to be re-rendered
            i = 0
            names, waves = list(wavethings.keys()), list(wavethings.values())
            for i in range(minimum):
                self.waveDocks[i].get_panel().set_waveThings(reset_parameters, wave=waves[i])
                self.waveDocks[i].setTitle(names[i])
            # Copy new wavethings into the interaction energy dock
            self._interactionEnergy_dock.get_panel().set_waveThings(reset_parameters, waves=wavethings, interactions=interactions)

            if len(self.waveDocks) < len(wavethings):
                # Need more wavepanels
                for j in range(i + 1, maximum):
                    waveDock = CustomDock(WavePanel(waves[j], names[j]), names[j])
                    self.waveDocks.append(waveDock)
                    # The interaction energy dock already exists as the last dock at the bottom, so add the new docks on top of the interaction energy dock.
                    self.wave_DockArea.addDock(waveDock, position='top', relativeTo=self._interactionEnergy_dock)
            elif len(self.waveDocks) > len(wavethings):
                # Need to delete wavepanels
                for j in range(i + 1, maximum):
                    # Remove the entry from the list
                    waveDock_to_delete = self.waveDocks.pop()
                    # Delete the widget from the GUI
                    waveDock_to_delete.close()
                    waveDock_to_delete.deleteLater()

            # If we're switching between settings, remake the dock area
            if (self.currentConfig is None) or ((self.currentConfig is not None) and (self.currentConfig != currentConfig)):
                # Clear the bar
                self._closedDockBar.clear()
                # Create a temporary dock area and fill it with the docks, then replace the existing dock area
                tempDock: DockArea = DockArea()
                for waveDock in self.waveDocks:
                    tempDock.addDock(waveDock)
                tempDock.addDock(self._interactionEnergy_dock)
                self.dock_layout.replaceWidget(self.wave_DockArea, tempDock)
                self.wave_DockArea.deleteLater()
                self.wave_DockArea = tempDock

                # Save the dock state (so that it can be set back to this "original" state at any time)
                # self._dock_state = self.wave_DockArea.saveState()
                self.currentConfig = currentConfig

        # Make all waveDocks saveable
        for waveDock in self.get_waveDocks():
            # try: waveDock.sigClosed.disconnect()
            # except: pass
            # Connect sigClosed of the Dock object to the method saveDock. The reconnect method ensures that the slot is
            # connected to the signal only once. The sigClosed signal is emitted when the close button in the upper left
            # corner is clicked.
            self.reconnect(waveDock.sigClosed, self.saveDock, self.saveDock)
            # waveDock.sigClosed.connect(self.saveDock)

    # Need to declare this slot here as a method, instead of a local function. Otherwise removing duplicates from the
    # list self._closed_docks with the reconnect method is not possible (for unknown reasons)
    def saveDock(self, dock: CustomDock):
        # self._closed_docks.append(dock)
        self._closedDockBar.add_dock(dock)
        # print(self._closed_docks)
        # for dock in self._closed_docks:
        #     print(f'id: {id(dock)}')

    def _make_saveState(self) -> dict:
        """Returns state in which docks are vertically arranged and have the same size"""
        dock_list = []
        size_list = []
        for waveDock in self.get_waveDocks():
            dock_list.append(('dock', waveDock.get_panel().get_name(), {}))
            size_list.append(1) # set vertical size to 1. the remaining space is automatically distributed among the docks
        return {'main': ('vertical', dock_list, {'sizes': size_list}), 'float': []}

    def reconnect(self, signal, newhandler=None, oldhandler=None):
        """
        Utility to disconnect/reconnect a slot from a signal
        (NB: the loop is needed for safely disconnecting a specific handler, because it may have been connected multple times, and disconnect(slot) only removes one connection at a time.).
        """
        try:
            if oldhandler is not None:
                while True:
                    signal.disconnect(oldhandler)
            else:
                signal.disconnect()
        except TypeError:
            pass
        if newhandler is not None:
            signal.connect(newhandler)