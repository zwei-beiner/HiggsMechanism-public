import numpy as np
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QWidget, QLabel, QHBoxLayout, QFrame

from main_package.gui.Simulation import Simulation


class Clock(QLabel):
    def __init__(self, simulation: Simulation):
        super().__init__()
        self._simulation = simulation

        self.update_time(0.0)
        self.setStyleSheet('font-size: 30pt;')

    def _sig_figs(self, x: float) -> int:
        """Returns the number of decimal places after the '.' in a floating point number x"""

        # Convert to string
        string: str = np.format_float_positional(x)
        # Get number of digits after the '.'
        number_of_digits: int = len(string.split('.')[1])
        return number_of_digits

    def update_time(self, current_dt_val: float):
        # Format time to show the same number of decimal places as in the widget dt
        decimal_places = self._sig_figs(current_dt_val)
        self.setText(f'Time: {self._simulation.get_time():.{decimal_places}f}')
