from collections import Counter

from main_package.gui.Setting import Setting
from main_package.simulation.FrozenMultiset import FrozenMultiset
from main_package.simulation.InitialConditioners import InitialConditioner, RSFInitialConditioner_PacketType, \
    RSFInitialConditioner_HiggsSillyType, RSFInitialConditioner_HiggsTypeMatchedToRSFAABBInteraction, \
    RSFInitialConditioner_HiggsTypeMatchedToFullyGeneralInteraction
from main_package.simulation.Integrator import LeapfrogIntegrator
from main_package.simulation.Interaction import RSFAABBInteraction, LeapfroggableInteraction, RSFAnBnInteraction, \
    RSFInteraction_FullGenerality
from main_package.simulation.PlotOptions import PlotOptions
from main_package.simulation.WaveThing import SpaceProperties, RSFProperties, RSFWaveThing


# # Abstract class
# class LeapfrogDemoSetter:
#     def messWith(self, integrator: LeapfrogIntegrator, initialConditioners: list[InitialConditioner]) -> None:
#         raise NotImplementedError()
#
#
# # Note: Graphics-related statements were omitted while porting this code from the Java version
# class FourWaveLeapfrogDemoSetter(LeapfrogDemoSetter):
#     # Constructors cannot be overloaded in Python so classmethods have to be used
#     def __init__(self, mode: int, n: int, L: float):
#         self.mode = mode
#         # self.narrower = narrower
#         # self.doSpecialZoomAboutNonNegVev = doSpecialZoomAboutNonNegVev
#         # self.specialZoomFactor = specialZoomFactor
#
#         self.n = n
#         self.L = L
#
#         self.initialConditioners: list[InitialConditioner] = []
#         self.integrator: LeapfrogIntegrator = LeapfrogIntegrator()
#         self._rsfWaveThings: dict[str, RSFWaveThing] = {}
#         self._parameters: dict[str, float] = {}
#
#     # @classmethod
#     # def demoSetter_without_specialZoom(cls, mode: int, narrower: float):
#     #     return cls(mode, narrower, False, 0)
#     #
#     # @classmethod
#     # def demoSetter_with_specialZoom(cls, mode: int, narrower: float, specialZoomFactor: float):
#     #     return cls(mode, narrower, True, specialZoomFactor)
#
#     # TODO: Remove code duplicates for waves A, B and C by using method
#     def messWith(self, integrator: LeapfrogIntegrator, initialConditioners: list[InitialConditioner]) -> None:
#         initialConditioners = self.initialConditioners
#         integrator = self.integrator
#
#         rsfWaveThings = self._rsfWaveThings
#
#         # self._waveNames = {'workspace': 'A', 'heavy control': 'B', 'higgs': 'C', 'photon': 'D'}
#         # self._parameters = {'A': (0, 0), 'B': (2.2 * 2.2, 0), 'C': (-2 * 2, 0.02 * 2 * 2), 'D': (0, 0),
#         #                     'lambda_AC': 0.1}
#         self._parameters = {'lambda_AC': 0.1}
#
#         ##################################################
#
#         n = self.n
#         L = self.L
#
#         propSpace = SpaceProperties(n, L)
#
#         pCentral = 4
#         pSpread = 1
#         centralX = 5
#         normScale = 6
#
#         ##################################################
#
#         rsfWaveThings['A'] = RSFWaveThing(propSpace, RSFProperties(mSq=0, quarticTerm=0))
#
#         initialConditioners.append(RSFInitialConditioner_PacketType(
#             rsfWaveThings['A'],
#             pCentral,
#             pSpread,
#             centralX,
#             normScale))
#         integrator.add(rsfWaveThings['A'])
#
#         ##################################################
#
#         rsfWaveThings['B'] = RSFWaveThing(propSpace, RSFProperties(mSq=2.2 * 2.2, quarticTerm=0))
#
#         initialConditioners.append(RSFInitialConditioner_PacketType(
#             rsfWaveThings['B'],
#             pCentral,
#             pSpread,
#             centralX,
#             normScale))
#         integrator.add(rsfWaveThings['B'])
#
#         ##################################################
#
#         narrower = 10
#         rsfWaveThings['C'] = RSFWaveThing(propSpace, RSFProperties(mSq=-2 * 2 * (narrower ** 2),
#                                                                    quarticTerm=0.02 * 2 * 2 * (narrower ** 2)))
#         # if self.doSpecialZoomAboutNonNegVev:
#         #
#         rsfAABBInteraction_AC = RSFAABBInteraction(rsfWaveThings['A'], rsfWaveThings['C'], 0.1)
#
#         if self.mode == 1 or self.mode == 2 or self.mode == 3:
#             initialConditioners.append(RSFInitialConditioner_HiggsSillyType(rsfWaveThings['C'], self.mode))
#         else:
#             initialConditioners.append(
#                 RSFInitialConditioner_HiggsTypeMatchedToRSFAABBInteraction(rsfWaveThings['C'], rsfWaveThings['A'],
#                                                                            rsfAABBInteraction_AC))
#
#         integrator.add(rsfWaveThings['C'])
#         integrator.add(rsfAABBInteraction_AC)
#
#         ##################################################
#
#         rsfWaveThings['D'] = RSFWaveThing(propSpace, RSFProperties(mSq=0, quarticTerm=0))
#
#         initialConditioners.append(
#             RSFInitialConditioner_PacketType(rsfWaveThings['D'], pCentral, pSpread, centralX, normScale))
#         integrator.add(rsfWaveThings['D'])



class InitialConditionBuilder:
    def __init__(self, n: int, L: float, p_central: float = 4, p_spread: float = 1, central_x: float = 5,
                 norm_scale: float = 6):
        self._n = n
        self._L = L
        self._propSpace = SpaceProperties(n, L)

        self._integrator: LeapfrogIntegrator = LeapfrogIntegrator()
        self._initialConditioners: dict[str, InitialConditioner] = {}
        self._rsfWaveThings: dict[str, RSFWaveThing] = {}
        self._parameters: dict[FrozenMultiset, LeapfroggableInteraction] = {}

        self._wave_params: tuple[float, ...] = (p_central, p_spread, central_x, norm_scale)

    def add_normal_wave(self, name: str, m_sq: float, q: float, plot_options: PlotOptions = None) -> None:
        self._rsfWaveThings[name] = RSFWaveThing(self._propSpace, RSFProperties(mSq=m_sq, quarticTerm=q), plot_options)
        self._initialConditioners[name] = RSFInitialConditioner_PacketType(self._rsfWaveThings[name],
                                                                           *self._wave_params)
        self._integrator.add(self._rsfWaveThings[name])

    def add_lambda_parameter(self, on: bool, val: float, *waves: str) -> None:
        key = FrozenMultiset(*waves, on=on)

        # interaction: LeapfroggableInteraction
        # # If the interaction is only between 2 waves, use the "hard-coded" interaction object for performance
        # if len(key.distinct_items()) == 2:
        #     wave_A: str = key.distinct_items()[0]
        #     n_A: int = key.frequencies()[0]
        #     wave_B: str = key.distinct_items()[1]
        #     n_B: int = key.frequencies()[1]
        #
        #     # interaction = RSFAnBnInteraction(self._rsfWaveThings[wave_A], n_A, self._rsfWaveThings[wave_B], n_B, val)
        #     interaction = RSFAABBInteraction(self._rsfWaveThings[wave_A], self._rsfWaveThings[wave_B], val)
        # else: # full generality
        waves: dict[str, tuple[int, RSFWaveThing]] = {string: (integer, self._rsfWaveThings[string]) for string, integer in key.to_dict().items()}
        interaction = RSFInteraction_FullGenerality(val, **waves)

        self._parameters[key] = interaction
        self._integrator.add(self._parameters[key])

    def set_higgs_wave(self, higgs_name: str, setting: Setting, other_name: str = None) -> None:
        if setting in (Setting.DEFAULT, Setting.INTERMEDIATE, Setting.ODD):

            # interaction = self._parameters[FrozenMultiset(higgs_name, higgs_name, other_name, other_name)]
            # if isinstance(interaction, RSFAABBInteraction):
            #     interaction = self._parameters[FrozenMultiset(higgs_name, higgs_name, other_name, other_name)]
            #     self._initialConditioners[higgs_name] = RSFInitialConditioner_HiggsTypeMatchedToRSFAABBInteraction(
            #         self._rsfWaveThings[higgs_name], self._rsfWaveThings[other_name],
            #         interaction
            #     )
            # else:
            self._initialConditioners[higgs_name] = RSFInitialConditioner_HiggsTypeMatchedToFullyGeneralInteraction(
                higgs_name, self._rsfWaveThings, self._parameters
            )

        elif setting in (Setting.SILLY1, Setting.SILLY2, Setting.SILLY3):
            num: int = {Setting.SILLY1: 1,
                        Setting.SILLY2: 2,
                        Setting.SILLY3: 3}[setting]
            self._initialConditioners[higgs_name] = RSFInitialConditioner_HiggsSillyType(self._rsfWaveThings[higgs_name], num)

    # def set_waves(self, *args) -> None:
    #     if len(args) % 3 != 0:
    #         raise Exception(f'Number of args not divisible by 3: {len(args)}')
    #
    #     number_of_waves: int = len(args) // 3
    #     for i in range(number_of_waves):
    #         name: str = args[3 * i]
    #         m_sq: float = args[3 * i + 1]
    #         q: float = args[3 * i + 2]
    #
    #         self._normal_wave(name, m_sq, q)

    # def set_params(self, *args) -> None:
    #     if len(args) % 2 != 0:
    #         raise Exception(f'Number of args not divisible by 2: {len(args)}')
    #
    #     number_of_params: int = len(args) // 2
    #     for i in range(number_of_params):
    #         name: str = args[2 * i]
    #         val: float = args[2 * i + 1]
    #         wave_A, wave_B = name.split('_')
    #         print(name.split('_'))
    #
    #         self._param(name, val, wave_A, wave_B)

    # def set_higgs(self, *args) -> None:
    #     if len(args) % 2 != 0:
    #         raise Exception(f'Number of args not divisible by 2: {len(args)}')
    #
    #     number_of_higgs_waves: int = len(args) // 2
    #     for i in range(number_of_higgs_waves):
    #         name: str = args[2 * i]
    #         init_mode: str = args[2 * i + 1]
    #
    #         self._higgs_wave(name, init_mode)

    def get_result(self) -> tuple[LeapfrogIntegrator,
                                  dict[str, InitialConditioner],
                                  dict[str, RSFWaveThing],
                                  dict[FrozenMultiset, LeapfroggableInteraction]]:
        return (self._integrator, self._initialConditioners, self._rsfWaveThings, self._parameters)


class InitialConditionDirector:
    def __init__(self, builder: InitialConditionBuilder):
        self._builder = builder

    def build_from_setting(self, setting):
        if setting == Setting.DEFAULT:
            self._default()
        elif setting == Setting.INTERMEDIATE:
            self._intermediate()
        elif setting == Setting.ODD:
            self._non_dispersive()
        elif setting == Setting.SILLY1:
            self._silly1()
        elif setting == Setting.SILLY2:
            self._silly2()
        elif setting == Setting.SILLY3:
            self._silly3()
        else:
            raise Exception(f'Setting not found: {setting}')

    # def _set_default_mass_terms(self, narrower: float):
    #     A = 'Test'
    #     B = 'Massive'
    #     C = 'Higgs'
    #     D = 'Photon'
    #
    #     self._builder.add_normal_wave(A, 0, 0)
    #     self._builder.add_normal_wave(B, 2.2 * 2.2, 0)
    #     self._builder.add_normal_wave(C, -2 * 2 * narrower ** 2, 0.02 * 2 * 2 * narrower ** 2)
    #     self._builder.add_normal_wave(D, 0, 0)
    #
    #     self._builder.add_lambda_parameter(True, 0.05, A, A, C, C)
    #     self._builder.add_lambda_parameter(False, 0.05, B, B, C, C)
    #     self._builder.add_lambda_parameter(False, 0.35, A, A, C)
    #     self._builder.add_lambda_parameter(False, 0.35, B, B, C)
    #     self._builder.add_lambda_parameter(False, 0.01, A, B, C)
    #
    #     self._builder.set_higgs_wave(C, Setting.DEFAULT)

    def _default(self):
        narrower = 10
        A = 'Test'
        B = 'Massive'
        C = 'Higgs'
        D = 'Photon'

        self._builder.add_normal_wave(A, 0, 0)
        self._builder.add_normal_wave(B, 2.2 * 2.2, 0)
        self._builder.add_normal_wave(C, -2 * 2 * narrower ** 2, 0.02 * 2 * 2 * narrower ** 2, PlotOptions(plot_range=(7.02, 7.08), energy_plot_range=(0, 0.05)))
        self._builder.add_normal_wave(D, 0, 0)

        self._builder.add_lambda_parameter(True, 0.05, A, A, C, C)
        self._builder.add_lambda_parameter(False, 0.05, B, B, C, C)
        self._builder.add_lambda_parameter(False, -0.35, A, A, C)
        self._builder.add_lambda_parameter(False, -0.35, B, B, C)
        self._builder.add_lambda_parameter(False, -0.01, A, B, C)

        self._builder.set_higgs_wave(C, Setting.DEFAULT)

    def _intermediate(self):
        narrower = 3
        A = 'Test'
        B = 'Massive'
        C = 'Higgs'
        D = 'Photon'

        self._builder.add_normal_wave(A, 0, 0)
        self._builder.add_normal_wave(B, 2.2 * 2.2, 0)
        self._builder.add_normal_wave(C, -2 * 2 * narrower ** 2, 0.02 * 2 * 2 * narrower ** 2, PlotOptions(plot_range=(6.3, 7.2), energy_plot_range=(0, 0.6)))
        self._builder.add_normal_wave(D, 0, 0)

        self._builder.add_lambda_parameter(True, 0.05, A, A, C, C)
        self._builder.add_lambda_parameter(False, 0.05, B, B, C, C)
        self._builder.add_lambda_parameter(False, -0.35, A, A, C)
        self._builder.add_lambda_parameter(False, -0.35, B, B, C)
        self._builder.add_lambda_parameter(False, -0.01, A, B, C)

        self._builder.set_higgs_wave(C, Setting.INTERMEDIATE)

    def _non_dispersive(self):
        narrower = 1
        A = 'Test'
        B = 'Massive'
        C = 'Higgs'
        D = 'Photon'

        self._builder.add_normal_wave(A, 0, 0)
        self._builder.add_normal_wave(B, 2.2 * 2.2, 0)
        self._builder.add_normal_wave(C, -2 * 2 * narrower ** 2, 0.02 * 2 * 2 * narrower ** 2, PlotOptions(plot_range=(-10, 10), energy_plot_range=(0, 50)))
        self._builder.add_normal_wave(D, 0, 0)

        self._builder.add_lambda_parameter(True, 0.05, A, A, C, C)
        self._builder.add_lambda_parameter(False, 0.05, B, B, C, C)
        self._builder.add_lambda_parameter(False, -0.35, A, A, C)
        self._builder.add_lambda_parameter(False, -0.35, B, B, C)
        self._builder.add_lambda_parameter(False, -0.01, A, B, C)

        self._builder.set_higgs_wave(C, Setting.INTERMEDIATE)

    def _silly1(self):
        narrower = 1
        A = 'Test'
        B = 'Massive'
        C = 'Higgs'
        D = 'Photon'

        self._builder.add_normal_wave(A, 0, 0)
        self._builder.add_normal_wave(B, 2.2 * 2.2, 0)
        self._builder.add_normal_wave(C, -2 * 2 * narrower ** 2, 0.02 * 2 * 2 * narrower ** 2)
        self._builder.add_normal_wave(D, 0, 0)

        self._builder.add_lambda_parameter(True, 0.05, A, A, C, C)
        self._builder.add_lambda_parameter(False, 0.05, B, B, C, C)
        self._builder.add_lambda_parameter(False, -0.35, A, A, C)
        self._builder.add_lambda_parameter(False, -0.35, B, B, C)
        self._builder.add_lambda_parameter(False, -0.01, A, B, C)

        self._builder.set_higgs_wave(C, Setting.SILLY1)

    def _silly2(self):
        narrower = 1
        A = 'Test'
        B = 'Massive'
        C = 'Higgs'
        D = 'Photon'

        self._builder.add_normal_wave(A, 0, 0)
        self._builder.add_normal_wave(B, 2.2 * 2.2, 0)
        self._builder.add_normal_wave(C, -2 * 2 * narrower ** 2, 0.02 * 2 * 2 * narrower ** 2, PlotOptions(energy_plot_range=(0, 12)))
        self._builder.add_normal_wave(D, 0, 0)

        self._builder.add_lambda_parameter(True, 0.05, A, A, C, C)
        self._builder.add_lambda_parameter(False, 0.05, B, B, C, C)
        self._builder.add_lambda_parameter(False, -0.35, A, A, C)
        self._builder.add_lambda_parameter(False, -0.35, B, B, C)
        self._builder.add_lambda_parameter(False, -0.01, A, B, C)

        self._builder.set_higgs_wave(C, Setting.SILLY2)

    def _silly3(self):
        narrower = 1
        A = 'Test'
        B = 'Massive'
        C = 'Higgs'
        D = 'Photon'

        self._builder.add_normal_wave(A, 0, 0)
        self._builder.add_normal_wave(B, 2.2 * 2.2, 0)
        self._builder.add_normal_wave(C, -2 * 2 * narrower ** 2, 0.02 * 2 * 2 * narrower ** 2, PlotOptions(energy_plot_range=(0, 12)))
        self._builder.add_normal_wave(D, 0, 0)

        self._builder.add_lambda_parameter(True, 0.05, A, A, C, C)
        self._builder.add_lambda_parameter(False, 0.05, B, B, C, C)
        self._builder.add_lambda_parameter(False, -0.35, A, A, C)
        self._builder.add_lambda_parameter(False, -0.35, B, B, C)
        self._builder.add_lambda_parameter(False, -0.01, A, B, C)

        self._builder.set_higgs_wave(C, Setting.SILLY3)