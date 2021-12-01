import math

import numpy as np
from numpy.polynomial.polynomial import Polynomial

from main_package.simulation.FrozenMultiset import FrozenMultiset
from main_package.simulation.Interaction import RSFAABBInteraction, LeapfroggableInteraction
from main_package.simulation.WaveThing import RSFWaveThing


# Abstract class
class InitialConditioner:
    def setInitialCondition(self) -> None:
        raise NotImplementedError()


class RSFInitialConditioner_PacketType(InitialConditioner):
    def __init__(self, rsfWaveThing: RSFWaveThing, pCentral: float, pSpread: float, centralX: float, normScale: float):
        self.rsfWaveThing = rsfWaveThing
        self.pCentral = pCentral
        self.pSpread = pSpread
        self.centralX = centralX
        self.normScale = normScale

    def setInitialCondition(self) -> None:
        """
        Method which initialises the wave packet by calling the function __putPacket.
        :return: None
        """
        self.rsfWaveThing.data.setToZero()
        self.__putPacket(self.rsfWaveThing, self.pCentral, self.pSpread, self.centralX, self.normScale)

    @staticmethod
    def __putPacket(wave: RSFWaveThing, centralP: float, pSpread: float, centralX: float, norm: float) -> None:
        """
        Function which initialises the position and initial velocity of a wave packet.
        Initialisation works by Fourier transforming a Gaussian in momentum space.
        See documentation for a derivation.
        :param wave: the wave to be initialised
        :param centralP: centre of the Gaussian in p-space
        :param pSpread: standard deviation of the Gaussian in p-space
        :param centralX: centre of the wave packet in x-space
        :param norm: scale factor by which the final wave packet is multiplied. Without it, the wave packet has a height of 1
        :return: None
        """

        # Calculate x_m, which has to be Fourier transformed to get the initial position wave.data.x

        # PhiTilde is the Gaussian envelope in p-space
        def PhiTilde(p: np.ndarray) -> np.ndarray:
            return 1 / (math.sqrt(2 * math.pi) * pSpread) * np.exp(-0.5 * ((p - centralP) / pSpread) ** 2)

        # modifiedPhiTilde is the modified function, which is fed into the FFT algorithm
        # Input: momentum p
        def modifiedPhiTilde(p: np.ndarray) -> np.ndarray:
            p_lower = centralP - (wave.propSpace.n // 2) * wave.propSpace.dp
            return PhiTilde(p + p_lower) * wave.propSpace.dp * np.exp(
                1.j * 2 * math.pi / wave.propSpace.n * (p * wave.propSpace.oneOverDp) * (
                        centralX // wave.propSpace.delta))

        # N-dimensional complex vector which is fed into the FFT algorithm
        x_m = modifiedPhiTilde(wave.propSpace.dp * np.arange(wave.propSpace.n))

        # Calculate v_m, which has to be Fourier transformed to get the initial velocity wave.data.v

        # Omega(m) is the dispersion relation
        # m: int array, which labels the lattice point by an integer (p = m * dp)
        def omega(m: np.ndarray) -> np.ndarray:
            return np.sqrt((np.sin(
                math.pi * m / wave.propSpace.n) * 2 / wave.propSpace.delta) ** 2 + wave.propRSF.mSq)

        # modifiedPhiTildeForV is the function which is fed into the FFT algorithm to get v
        def modifiedPhiTildeForV(p: np.ndarray) -> np.ndarray:
            p_lower = centralP - (wave.propSpace.n // 2) * wave.propSpace.dp
            return 1.j * omega((p + p_lower) / wave.propSpace.dp) * modifiedPhiTilde(p)

        # N-dimensional complex vector which is fed into the FFT algorithm
        v_m = modifiedPhiTildeForV(wave.propSpace.dp * np.arange(wave.propSpace.n))

        # pre_factor has to be multiplied with the result of the Fourier transforms
        # Shifts the wave packet to centralX in x-space. See documentation for details.
        pre_factor = norm * np.exp(
            -1.j * 2 * math.pi / wave.propSpace.n * (centralP * wave.propSpace.oneOverDp - wave.propSpace.n // 2) * (
                    np.arange(wave.propSpace.n) - centralX // wave.propSpace.delta))

        # Finally, the real part of the Fourier transforms is taken because we are working with real scalar fields
        wave.data.x = np.real(pre_factor * np.fft.fft(x_m))
        wave.data.v = np.real(pre_factor * np.fft.fft(v_m))


class RSFInitialConditioner_HiggsTypeMatchedToRSFAABBInteraction(InitialConditioner):
    def __init__(self, rsfWaveThing_higgs: RSFWaveThing, rsfWaveThing_other: RSFWaveThing,
                 rsfAABBInteraction: RSFAABBInteraction):
        self.rsfWaveThing_higgs = rsfWaveThing_higgs
        self.rsfWaveThing_other = rsfWaveThing_other
        self.rsfAABBInteraction = rsfAABBInteraction

    def setInitialCondition(self) -> None:
        # Assign values to the parameters here.
        # Reason: When the GUI is built, the values of mSq and quarticTerm will be modified to match the values in the GUI.
        # Only after that, the simulation is started and setInitialConditioner() is called.
        # Hence, self._alpha and self._beta have to be set in the setInitialConditioner() method.
        alpha: float = 0.25 * self.rsfWaveThing_higgs.get_propRSF().quarticTerm
        beta: float = -0.5 * self.rsfWaveThing_higgs.get_propRSF().mSq
        lambdaAABB: float = 0.5 * self.rsfAABBInteraction.get_parameter()

        phi: np.ndarray = self.rsfWaveThing_other.data.x
        phiDot: np.ndarray = self.rsfWaveThing_other.data.v

        inside = (beta - lambdaAABB * phi ** 2) / (2 * alpha)

        safeInside = np.where(inside >= 0, inside, 0)
        hTilde = np.sqrt(safeInside)

        self.rsfWaveThing_higgs.data.x = hTilde
        self.rsfWaveThing_higgs.data.v = - (lambdaAABB * phi * phiDot) / (2 * alpha * hTilde)


class RSFInitialConditioner_HiggsTypeMatchedToFullyGeneralInteraction(InitialConditioner):
    def __init__(self, higgs_name: str, rsfWaveThings: dict[str, RSFWaveThing], parameters: dict[FrozenMultiset, LeapfroggableInteraction]):
        """
        Input format: ( h, list[ (lambda_i, n_i, list[ (k_0, phi_0), (k_1, phi_1), ... ] ),... ] )
        """

        # Store REFERENCES to the wave things and the parameters.
        # If the parameters (mSq, quarticTerm, lambda parameters) are modified by the user in the GUI part of the code,
        # the objects stored here will show these changes. Hence, we can always read out the currently displayed values
        # here.
        # Instead, if we would only store the values (which are set in DemoSetter.py) we wouldn't be able to get the values
        # that are currently displayed in the GUI.
        # Therefore, we must store REFERENCES to rsfWaveThins and parameters.
        self._higgs_name: str = higgs_name
        self._rsfWaveThings: dict[str, RSFWaveThing] = rsfWaveThings
        self._parameters: dict[FrozenMultiset, LeapfroggableInteraction] = parameters

        # Here we declare variables which will be initialised in setInitialConditions().
        self._h: RSFWaveThing
        self._interaction_terms: list[tuple[float, int, list[tuple[int, RSFWaveThing],...]],...]
        self._N: int
        self._alpha: float
        self._beta: float

    def _find_all_interactions_with_higgs_field(self, higgs_name: str) -> tuple[RSFWaveThing, list[tuple[float, int, list[tuple[int, RSFWaveThing],...]],...]]:
        """
        Brings self._rsfWaveThings and self._parameters in the correct data format to be used for later processing.
        Here we search through all interactions and collect only those which interact with the higgs field.

        Returns all interaction terms which interact with the higgs field in the following format:
        \sum_{i=0}^{m-1} \lambda_i [\prod_{k=0}^{n-1} \phi_k^{k_i}] h^{n_i}
        then this method returns
        ( h, list[ (lambda_i, n_i, list[ (k_0, phi_0), (k_1, phi_1), ... ] ),... ] )
        """

        outer_list: list[tuple[float, int, list[tuple[int, RSFWaveThing],...]], ...] = []

        h: RSFWaveThing = self._rsfWaveThings[higgs_name]
        # Search through all interaction terms for those which have h in them AND are turned on (= have their checkbox ticked in the GUI)
        for frozenmultiset, interaction in self._parameters.items():
            if (higgs_name in frozenmultiset.distinct_items()) and frozenmultiset.get_on():
                lambda_i: float = interaction.get_parameter()
                n_i: int = 0

                other_names_and_exponents: list[tuple[int, str], ...] = []
                for name, exponent in frozenmultiset.to_dict().items():
                    if name != higgs_name:
                        other_names_and_exponents.append((exponent, name))
                    else:
                        n_i = exponent

                assert len(other_names_and_exponents) != 0
                assert n_i >= 1

                other_names_and_wavethings: list[tuple[int, RSFWaveThing],...] = [(exponent, self._rsfWaveThings[name]) for exponent, name in other_names_and_exponents]

                outer_list.append((lambda_i, n_i, other_names_and_wavethings))

        return (h, outer_list)

    def setInitialCondition(self) -> None:
        # Read out the data we need
        information_about_higgs_interaction_terms: tuple[RSFWaveThing, list[tuple[float, int, list[tuple[int, RSFWaveThing], ...]], ...]] = self._find_all_interactions_with_higgs_field(self._higgs_name)
        self._h: RSFWaveThing = information_about_higgs_interaction_terms[0]
        self._interaction_terms = information_about_higgs_interaction_terms[1]

        self._N: int = self._h.get_N()

        self._alpha: float = 0.25 * self._h.get_propRSF().quarticTerm
        self._beta: float = -0.5 * self._h.get_propRSF().mSq

        # Set the initial position and velocity of the Higgs field
        hTilde: np.ndarray = self._calculate_position()
        self._h.setX(hTilde)
        self._h.setV(self._calculate_velocity(hTilde))

    def _calculate_position(self) -> np.ndarray:
        # Pre-compute the expressions [\lambda_i * \prod_{k=0}^{n-1} \phi_k^{k_i}]
        interaction_terms_with_products_computed: list[tuple[int, np.ndarray],...] = []
        for lambda_i, n_i, list_of_waves in self._interaction_terms:
            prod: np.ndarray = np.ones(self._N)
            for k_i, phi_k in list_of_waves:
                prod = prod * (phi_k.getX() ** k_i)
            prod = lambda_i * prod
            interaction_terms_with_products_computed.append((n_i, prod))

        # Build the polynomial coefficients
        # Append the self-interaction terms
        interaction_terms_with_products_computed.append((4, self._alpha * np.ones(self._N)))
        interaction_terms_with_products_computed.append((2, -self._beta * np.ones(self._N)))
        # "Collapse" the coefficient list by grouping terms with the same power in h
        # Do this by walking through the list of terms and adding to a dictionary
        d: dict[int, np.ndarray] = {i: np.zeros(self._N) for i in range(max(interaction_terms_with_products_computed, key=lambda term: term[0])[0] + 1)}
        for n_i, prod in interaction_terms_with_products_computed:
            d[n_i] += prod

        # Calculate the minimum energy h(x) by minimising the potential for each value of x. Here, the index i plays the role of x
        hTilde: list[float] = [0 for _ in range(self._N)]
        for i in range(self._N):
            # Sort the dictionary according to its keys and take the ith element in the numpy array (which is the value of the wave phi at position x = i * deltaX)
            # We need to sort the dict because the numpy implementation of Polynomial takes the list of coefficients in ASCENDING powers
            coeff: list[float] = [d[k][i] for k in sorted(d)]
            p = Polynomial(coeff)
            # Choose a domain [-30, 30] without special reason
            hTilde[i] = self.find_global_min_of_polynomial(p, [-30, 30])

        return np.array(hTilde)

    def _calculate_velocity(self, hTilde: np.ndarray) -> np.ndarray:
        # Calculate numerator -------------------------------------------------------------------------

        def take_derivative(terms: list[tuple[int, RSFWaveThing],...]) -> np.ndarray:
            # Apply the product rule to differentiate a term of type d/dt (\prod \phi_k^{k_i})
            derivative_of_prod: np.ndarray = np.zeros(self._N)
            for i in range(len(terms)):
                prod: np.ndarray = np.ones(self._N)
                for j in range(len(terms)):
                    current_term: tuple[int, RSFWaveThing] = terms[i]
                    k_i: int = current_term[0]
                    phi: np.ndarray = current_term[1].getX()
                    phiDot: np.ndarray = current_term[1].getV()

                    if j == i:
                        prod = prod * k_i * (phi ** (k_i - 1)) * phiDot
                    else:
                        prod = prod * (phi ** k_i)

                derivative_of_prod = derivative_of_prod + prod
            return derivative_of_prod

        numerator: np.ndarray = np.zeros(self._N)
        for lambda_i, n_i, l in self._interaction_terms:
            derivative_of_prod: np.ndarray = take_derivative(l)
            numerator = numerator + n_i * lambda_i * derivative_of_prod * (hTilde ** (n_i - 1))

        # Calculate denominator -----------------------------------------------------------------------

        def compute_prod(terms: list[tuple[int, RSFWaveThing],...]) -> np.ndarray:
            prod: np.ndarray = np.ones(self._h.get_N())
            for k_i, phi_k in terms:
                prod = prod * (phi_k.getX() ** k_i)
            return prod

        sum: np.ndarray = np.zeros(self._N)
        for lambda_i , n_i, l in self._interaction_terms:
            if n_i > 1:
                sum = sum + n_i * (n_i - 1) * lambda_i * compute_prod(l) * (hTilde ** (n_i - 2))

        denominator: np.ndarray = 12 * self._alpha * (hTilde ** 2) - 2 * self._beta + sum

        # Calculate entire expression -----------------------------------------------------------------------
        res: np.ndarray = -numerator / denominator
        return res

    @staticmethod
    def find_global_min_of_polynomial(p: Polynomial, domain: list[float, float]) -> np.float64:
        """
        Finds argmin(p) where p is a polynomial, i.e. the value of x at which p(x) attains its global minimum.
        If multiple such argmin(p) are found, the largest one is returned.
        """

        # Derivative
        d: Polynomial = p.deriv()
        # Critical points
        r: np.ndarray = d.roots()
        # Take real roots. Tolerance for imaginary part is 1e-5
        real_valued: np.ndarray = r.real[abs(r.imag) < 1e-5]

        # Also consider the domain boundaries in the following, so append them to the list of critical points
        real_valued = np.concatenate((real_valued, np.array(domain)))

        # Find the minimum
        minimum: np.float = np.amin(p(real_valued))
        # Return the x value for which the polynomial p(x) is minimised. If there are multiple such values, take the largest one
        x_min = np.amax(real_valued[np.where(abs(p(real_valued) - minimum) < 1e-5)])

        # Silence the type warning by ensuring that x_min is a float (instead of array/int/complex)
        assert isinstance(x_min, np.float64)
        return x_min




class RSFInitialConditioner_HiggsSillyType(InitialConditioner):
    def __init__(self, rsfWaveThing_higgs: RSFWaveThing, mode: int):
        self.wave = rsfWaveThing_higgs
        self.mode = mode

    def __setZeroAndABit(self) -> None:
        """Manton instantons"""
        self.wave.data.x = np.full((self.wave.propSpace.n,), 0.00001)
        self.wave.data.v = np.zeros(self.wave.propSpace.n)

    def __setFlop(self) -> None:
        n = self.wave.propSpace.n
        i = np.arange(n)
        d = 2 * (np.cos(2 * np.pi * i / n) ** 2) - 1
        self.wave.data.x = self.wave.propRSF.nonNegVevOrZero() * d
        self.wave.data.v = np.zeros(n)

    def __setStaticFlop(self) -> None:
        w = self.wave.propRSF.nonNegVevOrZero()
        higgsMassOrOne = self.wave.propRSF.getQuadraticMassAboutVevOrZero()
        xPos = self.wave.propSpace.delta * np.arange(self.wave.propSpace.n)
        self.wave.data.x = w * np.tanh((xPos - 1. * self.wave.propSpace.L / 4.) * higgsMassOrOne / 2.) - w * np.tanh(
            (xPos - 3. * self.wave.propSpace.L / 4.) * higgsMassOrOne / 2.) - w
        self.wave.data.v = np.zeros(self.wave.propSpace.n)

    def setInitialCondition(self) -> None:
        if self.mode == 1:
            self.__setZeroAndABit()
        elif self.mode == 2:
            self.__setFlop()
        else:
            self.__setStaticFlop()
