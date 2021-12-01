class PlotOptions:
    """
    Stores plot options for an RSFWaveThing.
    Is created in DemoSetter.py and read out in WavePanel.py.
    """

    def __init__(self, plot_range: tuple[float, float] = None, colour_rgba: tuple[float, float, float, float] = None, energy_plot_range: tuple[float, float] = None):
        if plot_range is None:
            self._plot_range = [-11, 11]
        else:
            self._plot_range = plot_range

        if colour_rgba is None:
            # Set to black
            self._colour_rgba = (0,0,0, 255)
        else:
            self._colour_rgba = colour_rgba

        if energy_plot_range is None:
            self._energy_plot_range = [0, 60]
        else:
            self._energy_plot_range = energy_plot_range

    def get_plot_range(self) -> tuple[float, float]:
        return self._plot_range

    def get_colour_rgba(self) -> tuple[float, float, float, float]:
        return self._colour_rgba

    def get_energy_plot_range(self) -> tuple[float, float]:
        return self._energy_plot_range