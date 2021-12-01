import warnings

from PyQt5 import QtWidgets
import pyqtgraph as pg
import sys

from PyQt5.QtWidgets import QDesktopWidget

from main_package.gui.Simulation import Simulation
from main_package.gui.ui import UI


class Controller:
    def __init__(self):
        # When an error in the calculations of the simulation occurs, numpy raises a RuntimeWarning.
        # To catch these RuntimeWarnings, we must tell Python to convert the warnings to exceptions, which can be caught
        # and handled. This is a global setting and is therefore done here at the entry point of the program.
        warnings.filterwarnings('error')

        # Set initial conditions
        simulation = Simulation()

        # Pass initial conditions to GUI
        # GUI is made in the constructor
        ui = UI(simulation)
        # (Boilerplate code)
        ui.resize(*self.get_screen_size())
        ui.show()

        # Start the repainter loop
        ui.startRepainterThread()

        # Start the simulation loop
        ui.startSimulation()

    @staticmethod
    def get_screen_size() -> tuple[int, int]:
        dw: QDesktopWidget = QDesktopWidget()
        x: int = int(dw.width())
        y: int = int(dw.height())
        return (x, y)


def main():
    # (Boilerplate code)
    app = QtWidgets.QApplication(sys.argv)
    # Sets Fusion style because the Default style does not display small checkboxes correctly
    app.setStyle('Fusion')

    # Start app
    controller = Controller()

    # Pass control to the GUI event loop (boilerplate code)
    sys.exit(app.exec_())

if __name__=='__main__':
    main()