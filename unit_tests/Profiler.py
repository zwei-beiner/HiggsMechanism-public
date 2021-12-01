from PyQt5 import QtCore

from main_package.gui.Simulation import Simulation


class Main():
    def __init__(self):
        simulation = Simulation()
        simulation.run()

        # self._thread = QtCore.QThread()
        # simulation.moveToThread(self._thread)
        # self._thread.started.connect(simulation.run)
        # self._thread.start()

        # input('Press enter')
        # simulation.stop()
        # self._thread.quit()
        # self._thread.wait()


if __name__ == '__main__':
    Main()