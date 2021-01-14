# This code is part of Qiskit.
#
# (C) Copyright IBM 2018, 2020.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

""" Fermionic Operator """

from collections import defaultdict

import itertools
import logging
import sys

import numpy as np
from qiskit.quantum_info import Pauli
from qiskit.tools import parallel_map
from qiskit.tools.events import TextProgressBar

from qiskit.aqua import aqua_globals
from qiskit.aqua.operators import WeightedPauliOperator
from .qiskit_chemistry_error import QiskitChemistryError
from .bksf import bksf_mapping
from .particle_hole import particle_hole_transformation

logger = logging.getLogger(__name__)


class FermionicOperator:
    """
    A set of functions to map fermionic Hamiltonians to qubit Hamiltonians.

    References:

    - *E. Wigner and P. Jordan., Über das Paulische Äquivalenzverbot,
      Z. Phys., 47:631 (1928).*
    - *S. Bravyi and A. Kitaev. Fermionic quantum computation,
      Ann. of Phys., 298(1):210–226 (2002).*
    - *A. Tranter, S. Sofia, J. Seeley, M. Kaicher, J. McClean, R. Babbush,
      P. Coveney, F. Mintert, F. Wilhelm, and P. Love. The Bravyi–Kitaev
      transformation: Properties and applications. Int. Journal of Quantum
      Chemistry, 115(19):1431–1441 (2015).*
    - *S. Bravyi, J. M. Gambetta, A. Mezzacapo, and K. Temme,
      arXiv e-print arXiv:1701.08213 (2017).*
    - *K. Setia, J. D. Whitfield, arXiv:1712.00446 (2017)*

    """

    def __init__(self, h1, h2=None, ph_trans_shift=None):
        """
        This class requires the integrals stored in the '*chemist*' notation

             h2(i,j,k,l) --> adag_i adag_k a_l a_j

        and the integral values are used for the coefficients of the second-quantized
        Hamiltonian that is built. The integrals input here should be in block spin
        format and also have indexes reordered as follows 'ijkl->ljik'

        There is another popular notation, the '*physicist*' notation

             h2(i,j,k,l) --> adag_i adag_j a_k a_l

        If you are using the '*physicist*' notation, you need to convert it to
        the '*chemist*' notation. E.g. h2=numpy.einsum('ikmj->ijkm', h2)

        The :class:`~qiskit.chemistry.QMolecule` class has
        :attr:`~qiskit.chemistry.QMolecule.one_body_integrals` and
        :attr:`~qiskit.chemistry.QMolecule.two_body_integrals` properties that can be
        directly supplied to the `h1` and `h2` parameters here respectively.

        Args:
            h1 (numpy.ndarray): second-quantized fermionic one-body operator, a 2-D (NxN) tensor
            h2 (numpy.ndarray): second-quantized fermionic two-body operator,
                                a 4-D (NxNxNxN) tensor
            ph_trans_shift (float): energy shift caused by particle hole transformation
        """
        self._h1 = h1
        if h2 is None:
            #h2 = np.zeros((h1.shape[0], h1.shape[0], h1.shape[0], h1.shape[0]), dtype=h1.dtype)
            h2 = None
        self._h2 = h2
        self._ph_trans_shift = ph_trans_shift
        self._modes = self._h1.shape[0]
        self._map_type = None

    @property
    def modes(self):
        """Getter of modes."""
        return self._modes

    @property
    def h1(self):  # pylint: disable=invalid-name
        """Getter of one body integral tensor."""
        return self._h1

    @h1.setter
    def h1(self, new_h1):  # pylint: disable=invalid-name
        """Setter of one body integral tensor."""
        self._h1 = new_h1

    @property
    def h2(self):  # pylint: disable=invalid-name
        """Getter of two body integral tensor."""
        return self._h2

    @h2.setter
    def h2(self, new_h2):  # pylint: disable=invalid-name
        """Setter of two body integral tensor."""
        self._h2 = new_h2

    def __eq__(self, other):
        """Overload == ."""
        ret = np.all(self._h1 == other._h1)
        if not ret:
            return ret
        ret = np.all(self._h2 == other._h2)
        return ret

    def __ne__(self, other):
        """Overload != ."""
        return not self.__eq__(other)

    def transform(self, unitary_matrix):
        """Transform the one and two body term based on unitary_matrix."""
        self._h1_transform(unitary_matrix)
        self._h2_transform(unitary_matrix)

    def _h1_transform(self, unitary_matrix):
        """
        Transform h1 based on unitary matrix, and overwrite original property.

        Args:
            unitary_matrix (numpy.ndarray): A 2-D unitary matrix for h1 transformation.
        """
        self._h1 = unitary_matrix.T.conj().dot(self._h1.dot(unitary_matrix))

    def _h2_transform(self, unitary_matrix):
        """
        Transform h2 to get fermionic hamiltonian, and overwrite original property.

        Args:
            unitary_matrix (numpy.ndarray): A 2-D unitary matrix for h1 transformation.
        """
        num_modes = unitary_matrix.shape[0]
        temp_ret = np.zeros((num_modes, num_modes, num_modes, num_modes),
                            dtype=unitary_matrix.dtype)
        unitary_matrix_dagger = np.conjugate(unitary_matrix)

        # option 3: temp1 is a 3-D tensor, temp2 and temp3 are 2-D tensors
        # pylint: disable=unsubscriptable-object
        for a_i in range(num_modes):
            temp1 = np.einsum('i,i...->...', unitary_matrix_dagger[:, a_i], self._h2)
            for b_i in range(num_modes):
                temp2 = np.einsum('j,j...->...', unitary_matrix[:, b_i], temp1)
                temp3 = np.einsum('kc,k...->...c', unitary_matrix_dagger, temp2)
                temp_ret[a_i, b_i, :, :] = np.einsum('ld,l...c->...cd', unitary_matrix, temp3)

        self._h2 = temp_ret

    def _jordan_wigner_mode(self, n):
        r"""
        Jordan_Wigner mode.

        Each Fermionic Operator is mapped to 2 Pauli Operators, added together with the
        appropriate phase, i.e.:

        a_i\^\\dagger = Z\^i (X + iY) I\^(n-i-1) = (Z\^i X I\^(n-i-1)) + i (Z\^i Y I\^(n-i-1))
        a_i = Z\^i (X - iY) I\^(n-i-1)

        This is implemented by creating an array of tuples, each including two operators.
        The phase between two elements in a tuple is implicitly assumed, and added calculated at the
        appropriate time (see for example _one_body_mapping).

        Args:
            n (int): number of modes

        Returns:
            list[Tuple]: Pauli
        """
        a_list = []
        for i in range(n):
            a_z = np.asarray([1] * i + [0] + [0] * (n - i - 1), dtype=np.bool)
            a_x = np.asarray([0] * i + [1] + [0] * (n - i - 1), dtype=np.bool)
            b_z = np.asarray([1] * i + [1] + [0] * (n - i - 1), dtype=np.bool)
            b_x = np.asarray([0] * i + [1] + [0] * (n - i - 1), dtype=np.bool)
            a_list.append((Pauli(a_z, a_x), Pauli(b_z, b_x)))
        return a_list

    def _parity_mode(self, n):
        """
        Parity mode.

        Args:
            n (int): number of modes

        Returns:
            list[Tuple]: Pauli
        """
        a_list = []
        for i in range(n):
            a_z = [0] * (i - 1) + [1] if i > 0 else []
            a_x = [0] * (i - 1) + [0] if i > 0 else []
            b_z = [0] * (i - 1) + [0] if i > 0 else []
            b_x = [0] * (i - 1) + [0] if i > 0 else []
            a_z = np.asarray(a_z + [0] + [0] * (n - i - 1), dtype=np.bool)
            a_x = np.asarray(a_x + [1] + [1] * (n - i - 1), dtype=np.bool)
            b_z = np.asarray(b_z + [1] + [0] * (n - i - 1), dtype=np.bool)
            b_x = np.asarray(b_x + [1] + [1] * (n - i - 1), dtype=np.bool)
            a_list.append((Pauli(a_z, a_x), Pauli(b_z, b_x)))
        return a_list

    def _bravyi_kitaev_mode(self, n):
        """
        Bravyi-Kitaev mode.

        Args:
            n (int): number of modes

         Returns:
             numpy.ndarray: Array of mode indexes
        """

        def parity_set(j, n):
            """
            Computes the parity set of the j-th orbital in n modes.

            Args:
                j (int) : the orbital index
                n (int) : the total number of modes

            Returns:
                numpy.ndarray: Array of mode indexes
            """
            indexes = np.array([])
            if n % 2 != 0:
                return indexes

            if j < n / 2:
                indexes = np.append(indexes, parity_set(j, n / 2))
            else:
                indexes = np.append(indexes, np.append(
                    parity_set(j - n / 2, n / 2) + n / 2, n / 2 - 1))
            return indexes

        def update_set(j, n):
            """
            Computes the update set of the j-th orbital in n modes.

            Args:
                j (int) : the orbital index
                n (int) : the total number of modes

            Returns:
                numpy.ndarray: Array of mode indexes
            """
            indexes = np.array([])
            if n % 2 != 0:
                return indexes
            if j < n / 2:
                indexes = np.append(indexes, np.append(
                    n - 1, update_set(j, n / 2)))
            else:
                indexes = np.append(indexes, update_set(j - n / 2, n / 2) + n / 2)
            return indexes

        def flip_set(j, n):
            """
            Computes the flip set of the j-th orbital in n modes.

            Args:
                j (int) : the orbital index
                n (int) : the total number of modes

            Returns:
                numpy.ndarray: Array of mode indexes
            """
            indexes = np.array([])
            if n % 2 != 0:
                return indexes
            if j < n / 2:
                indexes = np.append(indexes, flip_set(j, n / 2))
            elif j >= n / 2 and j < n - 1:  # pylint: disable=chained-comparison
                indexes = np.append(indexes, flip_set(j - n / 2, n / 2) + n / 2)
            else:
                indexes = np.append(np.append(indexes, flip_set(
                    j - n / 2, n / 2) + n / 2), n / 2 - 1)
            return indexes

        a_list = []
        # FIND BINARY SUPERSET SIZE
        bin_sup = 1
        # pylint: disable=comparison-with-callable
        while n > np.power(2, bin_sup):
            bin_sup += 1
        # DEFINE INDEX SETS FOR EVERY FERMIONIC MODE
        update_sets = []
        update_pauli = []

        parity_sets = []
        parity_pauli = []

        flip_sets = []

        remainder_sets = []
        remainder_pauli = []
        for j in range(n):

            update_sets.append(update_set(j, np.power(2, bin_sup)))
            update_sets[j] = update_sets[j][update_sets[j] < n]

            parity_sets.append(parity_set(j, np.power(2, bin_sup)))
            parity_sets[j] = parity_sets[j][parity_sets[j] < n]

            flip_sets.append(flip_set(j, np.power(2, bin_sup)))
            flip_sets[j] = flip_sets[j][flip_sets[j] < n]

            remainder_sets.append(np.setdiff1d(parity_sets[j], flip_sets[j]))

            update_pauli.append(Pauli(np.zeros(n, dtype=np.bool), np.zeros(n, dtype=np.bool)))
            parity_pauli.append(Pauli(np.zeros(n, dtype=np.bool), np.zeros(n, dtype=np.bool)))
            remainder_pauli.append(Pauli(np.zeros(n, dtype=np.bool), np.zeros(n, dtype=np.bool)))
            for k in range(n):
                if np.in1d(k, update_sets[j]):
                    update_pauli[j].update_x(True, k)
                if np.in1d(k, parity_sets[j]):
                    parity_pauli[j].update_z(True, k)
                if np.in1d(k, remainder_sets[j]):
                    remainder_pauli[j].update_z(True, k)

            x_j = Pauli(np.zeros(n, dtype=np.bool), np.zeros(n, dtype=np.bool))
            x_j.update_x(True, j)
            y_j = Pauli(np.zeros(n, dtype=np.bool), np.zeros(n, dtype=np.bool))
            y_j.update_z(True, j)
            y_j.update_x(True, j)
            a_list.append((update_pauli[j] * x_j * parity_pauli[j],
                           update_pauli[j] * y_j * remainder_pauli[j]))
        return a_list

    def mapping(self, map_type, threshold=0.00000001):
        """
        Map fermionic operator to qubit operator.

        Using multiprocess to speedup the mapping, the improvement can be
        observed when h2 is a non-sparse matrix.

        Args:
            map_type (str): case-insensitive mapping type.
                            "jordan_wigner", "parity", "bravyi_kitaev", "bksf"
            threshold (float): threshold for Pauli simplification

        Returns:
            WeightedPauliOperator: create an Operator object in Paulis form.

        Raises:
            QiskitChemistryError: if the `map_type` can not be recognized.
        """

        # ###################################################################
        # ###########   DEFINING MAPPED FERMIONIC OPERATORS    ##############
        # ###################################################################

        self._map_type = map_type
        n = self._modes  # number of fermionic modes / qubits
        map_type = map_type.lower()
        if map_type == 'jordan_wigner':
            a_list = self._jordan_wigner_mode(n)
        elif map_type == 'parity':
            a_list = self._parity_mode(n)
        elif map_type == 'bravyi_kitaev':
            a_list = self._bravyi_kitaev_mode(n)
        elif map_type == 'bksf':
            return bksf_mapping(self)
        else:
            raise QiskitChemistryError('Please specify the supported modes: '
                                       'jordan_wigner, parity, bravyi_kitaev, bksf')

        # ###################################################################
        # ###########    BUILDING THE MAPPED HAMILTONIAN     ################
        # ###################################################################

        pauli_list = WeightedPauliOperator(paulis=[])
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Mapping one-body terms to Qubit Hamiltonian:")
            TextProgressBar(output_handler=sys.stderr)
        results = parallel_map(FermionicOperator._one_body_mapping,
                               [(self._h1[i, j], a_list[i], a_list[j])
                                for i, j in itertools.product(range(n), repeat=2)
                                if self._h1[i, j] != 0],
                               task_args=(threshold,), num_processes=aqua_globals.num_processes)
        for result in results:
            pauli_list += result
        pauli_list.chop(threshold=threshold)

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Mapping two-body terms to Qubit Hamiltonian:")
            TextProgressBar(output_handler=sys.stderr)
        results = parallel_map(FermionicOperator._two_body_mapping,
                               [(self._h2[i, j, k, m], a_list[i], a_list[j], a_list[k], a_list[m])
                                for i, j, k, m in itertools.product(range(n), repeat=4)
                                if self._h2[i, j, k, m] != 0],
                               task_args=(threshold,), num_processes=aqua_globals.num_processes)
        for result in results:
            pauli_list += result
        pauli_list.chop(threshold=threshold)

        if self._ph_trans_shift is not None:
            pauli_term = [self._ph_trans_shift, Pauli.from_label('I' * self._modes)]
            pauli_list += WeightedPauliOperator(paulis=[pauli_term])

        return pauli_list

    @staticmethod
    def _one_body_mapping(h1_ij_aij, threshold):
        """
        Subroutine for one body mapping.

        Args:
            h1_ij_aij (tuple): value of h1 at index (i,j), pauli at index i, pauli at index j
            threshold (float): threshold to remove a pauli

        Returns:
            WeightedPauliOperator: Operator for those paulis
        """
        h1_ij, a_i, a_j = h1_ij_aij
        pauli_list = []
        for alpha in range(2):
            for beta in range(2):
                pauli_prod = Pauli.sgn_prod(a_i[alpha], a_j[beta])
                coeff = h1_ij / 4 * pauli_prod[1] * np.power(-1j, alpha) * np.power(1j, beta)
                pauli_term = [coeff, pauli_prod[0]]
                if np.absolute(pauli_term[0]) > threshold:
                    pauli_list.append(pauli_term)
        return WeightedPauliOperator(paulis=pauli_list)

    @staticmethod
    def _two_body_mapping(h2_ijkm_a_ijkm, threshold):
        """
        Subroutine for two body mapping. We use the chemists notation
        for the two-body term, h2(i,j,k,m) adag_i adag_k a_m a_j.

        Args:
            h2_ijkm_a_ijkm (tuple): value of h2 at index (i,j,k,m),
                                    pauli at index i, pauli at index j,
                                    pauli at index k, pauli at index m
            threshold (float): threshold to remove a pauli

        Returns:
            WeightedPauliOperator: Operator for those paulis
        """
        h2_ijkm, a_i, a_j, a_k, a_m = h2_ijkm_a_ijkm
        pauli_list = []
        for alpha in range(2):
            for beta in range(2):
                for gamma in range(2):
                    for delta in range(2):
                        pauli_prod_1 = Pauli.sgn_prod(a_i[alpha], a_k[beta])
                        pauli_prod_2 = Pauli.sgn_prod(pauli_prod_1[0], a_m[gamma])
                        pauli_prod_3 = Pauli.sgn_prod(pauli_prod_2[0], a_j[delta])

                        phase1 = pauli_prod_1[1] * pauli_prod_2[1] * pauli_prod_3[1]
                        phase2 = np.power(-1j, alpha + beta) * np.power(1j, gamma + delta)
                        pauli_term = [h2_ijkm / 16 * phase1 * phase2, pauli_prod_3[0]]
                        if np.absolute(pauli_term[0]) > threshold:
                            pauli_list.append(pauli_term)
        return WeightedPauliOperator(paulis=pauli_list)

    def _convert_to_interleaved_spins(self):
        """
        Converting the spin order from block to interleaved.

        From up-up-up-up-down-down-down-down
        to up-down-up-down-up-down-up-down
        """
        # pylint: disable=unsubscriptable-object
        matrix = np.zeros((self._h1.shape), self._h1.dtype)
        n = matrix.shape[0]
        j = np.arange(n // 2)
        matrix[j, 2 * j] = 1.0
        matrix[j + n // 2, 2 * j + 1] = 1.0
        self.transform(matrix)

    def _convert_to_block_spins(self):
        """
        Converting the spin order from interleaved to block.

        From up-down-up-down-up-down-up-down
        to up-up-up-up-down-down-down-down
        """
        # pylint: disable=unsubscriptable-object
        matrix = np.zeros((self._h1.shape), self._h1.dtype)
        n = matrix.shape[0]
        j = np.arange(n // 2)
        matrix[2 * j, j] = 1.0
        matrix[2 * j + 1, n // 2 + j] = 1.0
        self.transform(matrix)

    # Modified for Open-Shell : 17.07.2019 by iso and bpa
    def particle_hole_transformation(self, num_particles):
        """
        The 'standard' second quantized Hamiltonian can be transformed in the
        particle-hole (p/h) picture, which makes the expansion of the trail wavefunction
        from the HF reference state more natural. In fact, for both trail wavefunctions
        implemented in q-lib ('heuristic' hardware efficient and UCCSD) the p/h Hamiltonian
        improves the speed of convergence of the VQE algorithm for the calculation of
        the electronic ground state properties.
        For more information on the p/h formalism see:
        P. Barkoutsos, arXiv:1805.04340(https://arxiv.org/abs/1805.04340).

        Args:
            num_particles (list, int): number of particles, if it is a list,
                                       the first number is alpha
                                       and the second number is beta.
        Returns:
            tuple: new_fer_op, energy_shift
        """

        self._convert_to_interleaved_spins()
        h_1, h_2, energy_shift = particle_hole_transformation(self._modes, num_particles,
                                                              self._h1, self._h2)
        new_fer_op = FermionicOperator(h1=h_1, h2=h_2, ph_trans_shift=energy_shift)
        new_fer_op._convert_to_block_spins()
        return new_fer_op, energy_shift

    def fermion_mode_elimination(self, fermion_mode_array):
        """Eliminate modes.

        Generate a new fermionic operator with the modes in fermion_mode_array deleted

        Args:
            fermion_mode_array (list): orbital index for elimination

        Returns:
            FermionicOperator: Fermionic Hamiltonian
        """
        # ROY INJECT
        return self.fermion_mode_elimination_fast(fermion_mode_array)
        # ROY INJECT

        fermion_mode_array = np.sort(fermion_mode_array)
        n_modes_old = self._modes
        n_modes_new = n_modes_old - fermion_mode_array.size
        mode_set_diff = np.setdiff1d(np.arange(n_modes_old), fermion_mode_array)
        h1_id_i, h1_id_j = np.meshgrid(mode_set_diff, mode_set_diff, indexing='ij')
        h1_new = self._h1[h1_id_i, h1_id_j].copy()
        if np.count_nonzero(self._h2) > 0:
            h2_id_i, h2_id_j, h2_id_k, h2_id_l = np.meshgrid(
                mode_set_diff, mode_set_diff, mode_set_diff, mode_set_diff, indexing='ij')
            h2_new = self._h2[h2_id_i, h2_id_j, h2_id_k, h2_id_l].copy()
        else:
            h2_new = np.zeros((n_modes_new, n_modes_new, n_modes_new, n_modes_new))
        return FermionicOperator(h1_new, h2_new)

    def fermion_mode_elimination_fast(self, fermion_mode_array):
        fermion_mode_array = np.sort(fermion_mode_array)
        n_modes_old = self._modes
        n_modes_new = n_modes_old - fermion_mode_array.size
        mode_set_diff = np.setdiff1d(np.arange(n_modes_old), fermion_mode_array)
        h1_id_i, h1_id_j = np.meshgrid(mode_set_diff, mode_set_diff, indexing='ij')
        h1_new = self._h1[h1_id_i, h1_id_j].copy()
        if np.count_nonzero(self._h2) > 0:
            #h2_id_i, h2_id_j, h2_id_k, h2_id_l = np.meshgrid(
            #    mode_set_diff, mode_set_diff, mode_set_diff, mode_set_diff, indexing='ij')
            #h2_new = self._h2[h2_id_i, h2_id_j, h2_id_k, h2_id_l].copy()
            #h2_new = defaultdict(float)
            #h2_new = np.zeros([mode_set_diff, mode_set_diff, mode_set_diff, mode_set_diff])
            h2_new = np.zeros((n_modes_new, n_modes_new, n_modes_new, n_modes_new))
            for ni, i in enumerate(mode_set_diff):
                for nj, j in enumerate(mode_set_diff):
                    for nk, k in enumerate(mode_set_diff):
                        for nl, l in enumerate(mode_set_diff):
                            h2_new[(ni, nj, nk, nl)] = self._h2[(i, j, k, l)]
        else:
            h2_new = np.zeros((n_modes_new, n_modes_new, n_modes_new, n_modes_new))
        return FermionicOperator(h1_new, h2_new)


    def fermion_mode_freezing(self, fermion_mode_array):
        """
        Freezing modes and extracting its energy.

        Generate a fermionic operator with the modes in fermion_mode_array deleted and
        provide the shifted energy after freezing.

        Args:
            fermion_mode_array (list): orbital index for freezing

        Returns:
            tuple(FermionicOperator, float):
                Fermionic Hamiltonian and energy of frozen modes
        """
        # ROY INJECT
        return self.fermion_mode_freezing_fast(fermion_mode_array)
        # ROY INJECT

        fermion_mode_array = np.sort(fermion_mode_array)
        n_modes_old = self._modes
        n_modes_new = n_modes_old - fermion_mode_array.size
        mode_set_diff = np.setdiff1d(np.arange(n_modes_old), fermion_mode_array)

        h_1 = self._h1.copy()
        h2_new = np.zeros((n_modes_new, n_modes_new, n_modes_new, n_modes_new))

        energy_shift = 0.0
        if np.count_nonzero(self._h2) > 0:
            # First simplify h2 and renormalize original h1
            for __i, __j, __l, __k in itertools.product(range(n_modes_old), repeat=4):
                # Untouched terms
                h2_ijlk = self._h2[__i, __j, __l, __k]
                if h2_ijlk == 0.0:
                    continue
                if __i in mode_set_diff and __j in mode_set_diff \
                        and __l in mode_set_diff and __k in mode_set_diff:
                    h2_new[__i - np.where(fermion_mode_array < __i)[0].size,
                           __j - np.where(fermion_mode_array < __j)[0].size,
                           __l - np.where(fermion_mode_array < __l)[0].size,
                           __k - np.where(fermion_mode_array < __k)[0].size] = h2_ijlk
                else:
                    if __i in fermion_mode_array:
                        if __l not in fermion_mode_array:
                            if __i == __k and __j not in fermion_mode_array:
                                h_1[__l, __j] -= h2_ijlk
                            elif __i == __j and __k not in fermion_mode_array:
                                h_1[__l, __k] += h2_ijlk
                        elif __i != __l:
                            if __j in fermion_mode_array and __i == __k and __l == __j:
                                energy_shift -= h2_ijlk
                            elif __l in fermion_mode_array and __i == __j and __l == __k:
                                energy_shift += h2_ijlk
                    elif __i not in fermion_mode_array and __l in fermion_mode_array:
                        if __l == __k and __j not in fermion_mode_array:
                            h_1[__i, __j] += h2_ijlk
                        elif __l == __j and __k not in fermion_mode_array:
                            h_1[__i, __k] -= h2_ijlk

        # now simplify h1
        energy_shift += np.sum(np.diagonal(h_1)[fermion_mode_array])
        h1_id_i, h1_id_j = np.meshgrid(mode_set_diff, mode_set_diff, indexing='ij')
        h1_new = h_1[h1_id_i, h1_id_j]

        return FermionicOperator(h1_new, h2_new), energy_shift

    def fermion_mode_freezing_fast(self, fermion_mode_array):
        """
        Freezing modes and extracting its energy.

        Generate a fermionic operator with the modes in fermion_mode_array deleted and
        provide the shifted energy after freezing.

        Args:
            fermion_mode_array (list): orbital index for freezing

        Returns:
            tuple(FermionicOperator, float):
                Fermionic Hamiltonian and energy of frozen modes
        """
        fermion_mode_array = np.sort(fermion_mode_array)
        n_modes_old = self._modes
        n_modes_new = n_modes_old - fermion_mode_array.size
        mode_set_diff = np.setdiff1d(np.arange(n_modes_old), fermion_mode_array)

        h_1 = self._h1.copy()
        #h2_new = np.zeros((n_modes_new, n_modes_new, n_modes_new, n_modes_new))
        h2_new = defaultdict(float)

        energy_shift = 0.0
        #if np.count_nonzero(self._h2) > 0:
        if len(self._h2.keys()) > 0:
            # First simplify h2 and renormalize original h1
            for __i, __j, __l, __k in itertools.product(range(n_modes_old), repeat=4):
                # Untouched terms
                h2_ijlk = self._h2[__i, __j, __l, __k]
                if h2_ijlk == 0.0:
                    continue
                if __i in mode_set_diff and __j in mode_set_diff \
                        and __l in mode_set_diff and __k in mode_set_diff:
                    h2_new[(__i - np.where(fermion_mode_array < __i)[0].size,
                           __j - np.where(fermion_mode_array < __j)[0].size,
                           __l - np.where(fermion_mode_array < __l)[0].size,
                           __k - np.where(fermion_mode_array < __k)[0].size)] = h2_ijlk
                else:
                    if __i in fermion_mode_array:
                        if __l not in fermion_mode_array:
                            if __i == __k and __j not in fermion_mode_array:
                                h_1[__l, __j] -= h2_ijlk
                            elif __i == __j and __k not in fermion_mode_array:
                                h_1[__l, __k] += h2_ijlk
                        elif __i != __l:
                            if __j in fermion_mode_array and __i == __k and __l == __j:
                                energy_shift -= h2_ijlk
                            elif __l in fermion_mode_array and __i == __j and __l == __k:
                                energy_shift += h2_ijlk
                    elif __i not in fermion_mode_array and __l in fermion_mode_array:
                        if __l == __k and __j not in fermion_mode_array:
                            h_1[__i, __j] += h2_ijlk
                        elif __l == __j and __k not in fermion_mode_array:
                            h_1[__i, __k] -= h2_ijlk

        # now simplify h1
        energy_shift += np.sum(np.diagonal(h_1)[fermion_mode_array])
        h1_id_i, h1_id_j = np.meshgrid(mode_set_diff, mode_set_diff, indexing='ij')
        h1_new = h_1[h1_id_i, h1_id_j]

        return FermionicOperator(h1_new, h2_new), energy_shift

    @staticmethod
    def construct_operator(molecule, freeze_list, remove_list):
        """ Returns a reduced fermionic operator

        Args:
            freeze_list: list of indices to freeze (in spin-space)
            remove_list: list of indices to remove (in spin-space)
        """
        # ROY INJECT

        # Sanity check; freeze and remove lists should be unique
        assert len(freeze_list + remove_list) == len(set(freeze_list + remove_list))
        freeze_list = np.array(list(freeze_list))
        remove_list = np.array(list(remove_list))

        # Print molecular info
        print("Alpha orbitals:")
        print("=" * 20)
        print(molecule.orbital_energies)
        print()

        print("Beta orbitals:")
        print("=" * 20)
        print(molecule.orbital_energies_b)
        print()

        # ================
        # ---- REMOVE ----
        # ================

        N = molecule.num_orbitals  # length of integral arrays (spinless)

        # Convert remove_list to spin-less masks
        remove_list_a = remove_list[remove_list < N]
        remove_list_b = remove_list[remove_list >= N] - N

        # Invert to get non-removed lists
        nonremoved_list_a = np.setdiff1d(np.arange(0, N), remove_list_a)
        nonremoved_list_b = np.setdiff1d(np.arange(0, N), remove_list_b)

        # Prep h1 and h2 arrays

        print("Preparing h1 arrays...")

        # H1
        ints_a = molecule.mo_onee_ints

        if molecule.mo_onee_ints_b is None:
            ints_b = ints_a
        else:
            ints_b = molecule.mo_onee_ints_b

        print("Preparing h2 arrays...")
        
        # H2
        ints_aa = np.einsum('ijkl->ljik', molecule.mo_eri_ints)

        if molecule.mo_eri_ints_bb is None or molecule.mo_eri_ints_ba is None:
            ints_bb = ints_ba = ints_ab = ints_aa
        else:
            ints_bb = np.einsum('ijkl->ljik', molecule.mo_eri_ints_bb)
            ints_ba = np.einsum('ijkl->ljik', molecule.mo_eri_ints_ba)
            ints_ab = np.einsum('ijkl->ljik', molecule.mo_eri_ints_ba.transpose())

        new_N = len(nonremoved_list_a) + len(nonremoved_list_b)

        ######## REMOVE STEP FOR H1 ########

        print("Starting remove step for h1...")

        mask_h1_a_i, mask_h1_a_j = np.meshgrid(nonremoved_list_a, nonremoved_list_a, indexing='ij', sparse=True)
        mask_h1_b_i, mask_h1_b_j = np.meshgrid(nonremoved_list_b, nonremoved_list_b, indexing='ij', sparse=True)

        h1_a = ints_a[mask_h1_a_i, mask_h1_a_j]
        h1_b = ints_b[mask_h1_b_i, mask_h1_b_j]

        # Merge h1_a and h1_b

        #               ---------------------
        #               |         |xxxxxxxxx|
        #               |  h1_a   |xxxxxxxxx|
        #               |         |xxxxxxxxx|
        #               |         |xxxxxxxxx|
        #    h_1 =      |-------------------|
        #               |xxxxxxxxx|         |
        #               |xxxxxxxxx|  h1_b   |
        #               |xxxxxxxxx|         |
        #               |xxxxxxxxx|         |
        #               |-------------------|

        # Note that h1_a = h1_b for RHF since the Hamiltonian is spin-agnostic

        h1 = np.zeros((new_N, new_N))
        h1[:len(nonremoved_list_a), :len(nonremoved_list_a)] = h1_a
        h1[len(nonremoved_list_a):, len(nonremoved_list_a):] = h1_b

        ######## REMOVE STEP FOR H2 ########

        print("Starting remove step for h2...")

        mask_h2_aa_i, mask_h2_aa_j, mask_h2_aa_k, mask_h2_aa_l = np.meshgrid(nonremoved_list_a, nonremoved_list_a, nonremoved_list_a, nonremoved_list_a, indexing='ij', sparse=True)
        mask_h2_ab_i, mask_h2_ab_j, mask_h2_ab_k, mask_h2_ab_l = np.meshgrid(nonremoved_list_a, nonremoved_list_b, nonremoved_list_b, nonremoved_list_a, indexing='ij', sparse=True)
        mask_h2_ba_i, mask_h2_ba_j, mask_h2_ba_k, mask_h2_ba_l = np.meshgrid(nonremoved_list_b, nonremoved_list_a, nonremoved_list_a, nonremoved_list_b, indexing='ij', sparse=True)
        mask_h2_bb_i, mask_h2_bb_j, mask_h2_bb_k, mask_h2_bb_l = np.meshgrid(nonremoved_list_b, nonremoved_list_b, nonremoved_list_b, nonremoved_list_b, indexing='ij', sparse=True)

        h2_aa = ints_aa[mask_h2_aa_i, mask_h2_aa_j, mask_h2_aa_k, mask_h2_aa_l]
        h2_ba = ints_ba[mask_h2_ba_i, mask_h2_ba_j, mask_h2_ba_k, mask_h2_ba_l]
        h2_ab = ints_ab[mask_h2_ab_i, mask_h2_ab_j, mask_h2_ab_k, mask_h2_ab_l]
        h2_bb = ints_bb[mask_h2_bb_i, mask_h2_bb_j, mask_h2_bb_k, mask_h2_bb_l]

        # Merge h2_aa, h2_bb, h2_ab, h2_ba

        #               ---------------------
        #               |         |         |
        #               |  h2_aa  |  h2_ab  |
        #               |         |         |
        #               |         |         |
        #    h_2 =      |-------------------|
        #               |         |         |
        #               |  h2_ba  |  h2_bb  |
        #               |         |         |
        #               |         |         |
        #               |-------------------|

        # Note that in h_pqrs, p == s, q == r, else or the value is 0

        h2 = np.zeros((new_N, new_N, new_N, new_N))
        num_a = len(nonremoved_list_a)
        h2[:num_a, :num_a, :num_a, :num_a] = h2_aa
        h2[:num_a, num_a:, num_a:, :num_a] = h2_ab
        h2[num_a:, :num_a, :num_a, num_a:] = h2_ba
        h2[num_a:, num_a:, num_a:, num_a:] = h2_bb

        h2 *= -0.5  # Normalization to match qiskit's twoe_to_spin function

        #return FermionicOperator(h1, h2), 0

        ######## CLEANUP BEFORE FREEZING ########

        freeze_list = [x - np.where(remove_list < x)[0].size for x in freeze_list]

        # At this point, h1/h2 should be spin-agnostic
        # (as in, we can think of them as generic MO arrays without considering spin)
        # If it is not then you did something wrong

        # ================
        # ---- FREEZE ----
        # ================

        print("Starting freeze steps...")
        N = h1.shape[0]  # Use new h1 size

        # Invert to get non-frozen lists
        nonfreeze_list = np.setdiff1d(np.arange(0, N), freeze_list)

        # In order to remove frozen orbitals, there are effectively three subspaces we care about:
        #
        #
        #
        #               ---------------------------
        #               |               |         |
        #               |  frozen       |  mid-   |
        #               |               |  frozen |
        #               |               |         |
        #               |               |         |
        #               |               |         |
        #    h_2 =      |-------------------------|
        #               |               |         |
        #               |  mid-         |  non-   |
        #               |  frozen       |  frozen |
        #               |               |         |
        #               |-------------------------|

        #   (Note that this is a 2d projection of a 4d matrix, but it generalizes)
        #
        #   The "non-frozen" subspace is effectively moved to the new array untouched
        #   The "frozen" subspace is removed and contributes to the energy_shift constant
        #   The "mid-frozen" subspace is removed and contributes to various rows in h_1

        energy_shift = 0

        ######## HANDLE NON-FROZEN SUBSPACE ########
        from time import time
        print("Start non-frozen...")
        start_time = time()
        mask_f_h2_i, mask_f_h2_j, mask_f_h2_k, mask_f_h2_l = np.meshgrid(nonfreeze_list, nonfreeze_list, nonfreeze_list, nonfreeze_list, indexing='ij', sparse=True)
        h2_new = h2[mask_f_h2_i, mask_f_h2_j, mask_f_h2_k, mask_f_h2_l]
        print("Done. Time: " + str(time() - start_time))

        ######## HANDLE FROZEN SUBSPACE ########
        print("Start frozen...")
        start_time = time()
        for __i, __l in itertools.product(freeze_list, repeat=2):
            # i == k, j == l, i != l:   energy -= h[ijlk]
            # i == j, k == l, i != l:   energy += h[ijlk]
            if __i == __l:
                continue

            energy_shift -= h2[__i, __l, __l, __i]
            energy_shift += h2[__i, __i, __l, __l]
        print("Done. Time: " + str(time() - start_time))

        ######## HANDLE MID-FROZEN SUBSPACE ########

        print("Start mid-frozen...")
        start_time = time()
        for x, y in itertools.product(nonfreeze_list, repeat=2):
            for z in freeze_list:
                h1[x, y] -= h2[z, x, y, z]
                h1[x, y] += h2[z, z, x, y]
                h1[x, y] += h2[x, y, z, z]
                h1[x, y] -= h2[x, z, z, y]
        print("Done. Time: " + str(time() - start_time))

        ######## HANDLE H1 ########
        energy_shift += np.sum(np.diagonal(h1)[freeze_list])

        h1_id_i, h1_id_j = np.meshgrid(nonfreeze_list, nonfreeze_list, indexing='ij', sparse=True)
        h1_new = h1[h1_id_i, h1_id_j]

        h1 = h1_new
        h2 = h2_new

        return FermionicOperator(h1, h2), energy_shift



        # ---- Calculate integrals ----
        h1 = molecule.one_body_integrals
        #h2 = molecule.two_body_integrals

        N = molecule.mo_eri_ints.size
        mohijkl = molecule.mo_eri_ints
        mohijkl_bb = molecule.mo_eri_ints_bb
        mohijkl_ba = molecule.mo_eri_ints_ba

        # Perform reduction
        all_list = np.arange(N)
        nonremove_list = np.setdiff1d(all_list, remove_list)

        mohijkl = mohijkl[nonremove_list]
        mohijkl_bb = mohijkl_bb[nonremove_list]
        mohijkl_ba = mohijkl_ba[nonremove_list]
        1/0

        # Update freeze_list to account for missing orbitals
        for i in freeze_list:
            ni = i - np.where(remove_list < i)[0].size
            freeze_list.append(ni)

        ints_aa = np.einsum('ijkl->ljik', mohijkl)

        if mohijkl_bb is None or mohijkl_ba is None:
            ints_bb = ints_ba = ints_ab = ints_aa
        else:
            ints_bb = np.einsum('ijkl->ljik', mohijkl_bb)
            ints_ba = np.einsum('ijkl->ljik', mohijkl_ba)
            ints_ab = np.einsum('ijkl->ljik', mohijkl_ba.transpose())

        # The number of spin orbitals is twice the number of orbitals
        norbs = mohijkl.shape[0]
        nspin_orbs = 2*norbs

        # Two electron terms
        moh2_qubit = defaultdict(float)
        for p in range(nspin_orbs):  # pylint: disable=invalid-name
            for q in range(nspin_orbs):
                for r in range(nspin_orbs):
                    for s in range(nspin_orbs):  # pylint: disable=invalid-name
                        spinp = int(p/norbs)
                        spinq = int(q/norbs)
                        spinr = int(r/norbs)
                        spins = int(s/norbs)
                        if spinp != spins:
                            continue
                        if spinq != spinr:
                            continue
                        if spinp == 0:
                            ints = ints_aa if spinq == 0 else ints_ba
                        else:
                            ints = ints_ab if spinq == 0 else ints_bb
                        orbp = int(p % norbs)
                        orbq = int(q % norbs)
                        orbr = int(r % norbs)
                        orbs = int(s % norbs)
                        if abs(ints[orbp, orbq, orbr, orbs]) > threshold:
                            moh2_qubit[(p, q, r, s)] = -0.5*ints[orbp, orbq, orbr, orbs]

        return moh2_qubit




        # ---- Full thing (probably remove) ----
        fermop = FermionicOperator(h1=h1)

        # ---- My sad attempt at a rewrite ----
        freeze_list = np.sort(freeze_list)
        remove_list = np.sort(remove_list)


        active_list = np.setdiff1d(np.arange(fermop._modes), freeze_list, remove_list)
        N = len(active_list)

        h1_new = np.zeroes((N, N))
        h2_new = np.zeros((N, N, N, N))

        for __i, __l in itertools.product(active_list, repeat=2):
            for __j, __k in itertools.product(active_list, repeat=2):
                h2_new[(__i - np.where(freeze_list < __i)[0].size,
                       __j - np.where(freeze_list < __j)[0].size,
                       __l - np.where(freeze_list < __l)[0].size,
                       __k - np.where(freeze_list < __k)[0].size)] = get_h2(__i, __j, __l, __k)
        print("Done. Time: " + str(time() - start_time))

        print("Start frozen...")
        start_time = time()
        # Everything is frozen
        for __i, __l in itertools.product(freeze_list, repeat=2):
            # i == k, j == l, i != l:   energy -= h[ijlk]
            # i == j, k == l, i != l:   energy += h[ijlk]
            if __i == __l:
                continue

            energy_shift -= get_h2(__i, __l, __l, __i)
            energy_shift += get_h2(__i, __i, __l, __l)
        print("Done. Time: " + str(time() - start_time))

        print("Start else-frozen...")
        start_time = time()
        for x, y in itertools.product(active_list, repeat=2):
            for z in freeze_list:
                h_1[x, y] -= get_h2(z, x, y, z)
                h_1[x, y] += get_h2(z, z, x, y)
                h_1[x, y] += get_h2(x, y, z, z)
                h_1[x, y] -= get_h2(x, z, z, y)
        print("Done. Time: " + str(time() - start_time))

        # now simplify h1
        energy_shift += np.sum(np.diagonal(h_1)[freeze_list])
        h1_id_i, h1_id_j = np.meshgrid(active_list, active_list, indexing='ij')
        h1_new = h_1[h1_id_i, h1_id_j]

        fermop = FermionicOperator(h1_new, h2_new)



        # ---- Freezing ----
        #freeze_list = np.sort(freeze_list)

        #n_modes_old = fermop._modes
        #n_modes_new = n_modes_old - freeze_list.size
        #active_list = np.setdiff1d(np.arange(n_modes_old), freeze_list, remove_list)

        #h_1 = fermop._h1.copy()
        #h2_new = np.zeros((n_modes_new, n_modes_new, n_modes_new, n_modes_new))

        ## Precomputation step
        #from time import time

        #print("Start precompute...")
        #start_time = time()
        #molecule.two_body_integral_precompute()
        #print("Done. Time: " + str(time() - start_time))

        #def get_h2(i, j, l, k):
        #    """ Get a single point of the two_body_integrals tensor """
        #    return molecule.two_body_integral_single(i, j, l, k)

        #energy_shift = 0.0
        ##if np.count_nonzero(fermop._h2) > 0:
        ##if len(fermop._h2.keys()) > 0:
        #if True:
        #    # First simplify h2 and renormalize original h1

        #    # Everything is non-frozen
        #    print("Start non-frozen...")
        #    start_time = time()
        #    for __i, __l in itertools.product(active_list, repeat=2):
        #        for __j, __k in itertools.product(active_list, repeat=2):
        #            h2_new[(__i - np.where(freeze_list < __i)[0].size,
        #                   __j - np.where(freeze_list < __j)[0].size,
        #                   __l - np.where(freeze_list < __l)[0].size,
        #                   __k - np.where(freeze_list < __k)[0].size)] = get_h2(__i, __j, __l, __k)
        #    print("Done. Time: " + str(time() - start_time))

        #    print("Start frozen...")
        #    start_time = time()
        #    # Everything is frozen
        #    for __i, __l in itertools.product(freeze_list, repeat=2):
        #        # i == k, j == l, i != l:   energy -= h[ijlk]
        #        # i == j, k == l, i != l:   energy += h[ijlk]
        #        if __i == __l:
        #            continue

        #        energy_shift -= get_h2(__i, __l, __l, __i)
        #        energy_shift += get_h2(__i, __i, __l, __l)
        #    print("Done. Time: " + str(time() - start_time))

        #    print("Start else-frozen...")
        #    start_time = time()
        #    for x, y in itertools.product(active_list, repeat=2):
        #        for z in freeze_list:
        #            h_1[x, y] -= get_h2(z, x, y, z)
        #            h_1[x, y] += get_h2(z, z, x, y)
        #            h_1[x, y] += get_h2(x, y, z, z)
        #            h_1[x, y] -= get_h2(x, z, z, y)
        #    print("Done. Time: " + str(time() - start_time))

        ## now simplify h1
        #energy_shift += np.sum(np.diagonal(h_1)[freeze_list])
        #h1_id_i, h1_id_j = np.meshgrid(active_list, active_list, indexing='ij')
        #h1_new = h_1[h1_id_i, h1_id_j]

        #fermop = FermionicOperator(h1_new, h2_new)

        ## ---- ELIMINATION ----
        #print("Start elimination...")
        #start_time = time()

        #remove_list = np.sort(remove_list)
        #n_modes_old = fermop._modes
        #n_modes_new = n_modes_old - remove_list.size
        #nonremove_list = np.setdiff1d(np.arange(n_modes_old), remove_list)
        #h1_id_i, h1_id_j = np.meshgrid(nonremove_list, nonremove_list, indexing='ij')
        #h1_new = fermop._h1[h1_id_i, h1_id_j].copy()
        #if np.count_nonzero(fermop._h2) > 0:
        #    h2_id_i, h2_id_j, h2_id_k, h2_id_l = np.meshgrid(
        #        nonremove_list, nonremove_list, nonremove_list, nonremove_list, indexing='ij')
        #    h2_new = fermop._h2[h2_id_i, h2_id_j, h2_id_k, h2_id_l].copy()
        #else:
        #    h2_new = np.zeros((n_modes_new, n_modes_new, n_modes_new, n_modes_new))

        #print("Done. Time: " + str(time() - start_time))
        #fermop = FermionicOperator(h1_new, h2_new)

        return fermop, energy_shift


    def total_particle_number(self):
        """
        A data_preprocess_helper fermionic operator which can be used to evaluate the number of
        particle of the given eigenstate.

        Returns:
            FermionicOperator: Fermionic Hamiltonian
        """
        modes = self._modes
        h_1 = np.eye(modes, dtype=np.complex)
        h_2 = np.zeros((modes, modes, modes, modes))
        return FermionicOperator(h_1, h_2)

    def total_magnetization(self):
        """
        A data_preprocess_helper fermionic operator which can be used to
        evaluate the magnetization of the given eigenstate.

        Returns:
            FermionicOperator: Fermionic Hamiltonian
        """
        modes = self._modes
        h_1 = np.eye(modes, dtype=np.complex) * 0.5
        h_1[modes // 2:, modes // 2:] *= -1.0
        h_2 = np.zeros((modes, modes, modes, modes))
        return FermionicOperator(h_1, h_2)

    def _s_x_squared(self):
        """

        Returns:
            FermionicOperator: Fermionic Hamiltonian
        """
        num_modes = self._modes
        num_modes_2 = num_modes // 2
        h_1 = np.zeros((num_modes, num_modes))
        h_2 = np.zeros((num_modes, num_modes, num_modes, num_modes))

        for p, q in itertools.product(range(num_modes_2), repeat=2):  # pylint: disable=invalid-name
            if p != q:
                h_2[p, p + num_modes_2, q, q + num_modes_2] += 1.0
                h_2[p + num_modes_2, p, q, q + num_modes_2] += 1.0
                h_2[p, p + num_modes_2, q + num_modes_2, q] += 1.0
                h_2[p + num_modes_2, p, q + num_modes_2, q] += 1.0
            else:
                h_2[p, p + num_modes_2, p, p + num_modes_2] -= 1.0
                h_2[p + num_modes_2, p, p + num_modes_2, p] -= 1.0
                h_2[p, p, p + num_modes_2, p + num_modes_2] -= 1.0
                h_2[p + num_modes_2, p + num_modes_2, p, p] -= 1.0

                h_1[p, p] += 1.0
                h_1[p + num_modes_2, p + num_modes_2] += 1.0

        h_1 *= 0.25
        h_2 *= 0.25
        return h_1, h_2

    def _s_y_squared(self):
        """

        Returns:
            FermionicOperator: Fermionic Hamiltonian
        """
        num_modes = self._modes
        num_modes_2 = num_modes // 2
        h_1 = np.zeros((num_modes, num_modes))
        h_2 = np.zeros((num_modes, num_modes, num_modes, num_modes))

        for p, q in itertools.product(range(num_modes_2), repeat=2):  # pylint: disable=invalid-name
            if p != q:
                h_2[p, p + num_modes_2, q, q + num_modes_2] -= 1.0
                h_2[p + num_modes_2, p, q, q + num_modes_2] += 1.0
                h_2[p, p + num_modes_2, q + num_modes_2, q] += 1.0
                h_2[p + num_modes_2, p, q + num_modes_2, q] -= 1.0
            else:
                h_2[p, p + num_modes_2, p, p + num_modes_2] += 1.0
                h_2[p + num_modes_2, p, p + num_modes_2, p] += 1.0
                h_2[p, p, p + num_modes_2, p + num_modes_2] -= 1.0
                h_2[p + num_modes_2, p + num_modes_2, p, p] -= 1.0

                h_1[p, p] += 1.0
                h_1[p + num_modes_2, p + num_modes_2] += 1.0

        h_1 *= 0.25
        h_2 *= 0.25
        return h_1, h_2

    def _s_z_squared(self):
        """

        Returns:
            FermionicOperator: Fermionic Hamiltonian
        """
        num_modes = self._modes
        num_modes_2 = num_modes // 2
        h_1 = np.zeros((num_modes, num_modes))
        h_2 = np.zeros((num_modes, num_modes, num_modes, num_modes))

        for p, q in itertools.product(range(num_modes_2), repeat=2):  # pylint: disable=invalid-name
            if p != q:
                h_2[p, p, q, q] += 1.0
                h_2[p + num_modes_2, p + num_modes_2, q, q] -= 1.0
                h_2[p, p, q + num_modes_2, q + num_modes_2] -= 1.0
                h_2[p + num_modes_2, p + num_modes_2,
                    q + num_modes_2, q + num_modes_2] += 1.0
            else:
                h_2[p, p + num_modes_2, p + num_modes_2, p] += 1.0
                h_2[p + num_modes_2, p, p, p + num_modes_2] += 1.0

                h_1[p, p] += 1.0
                h_1[p + num_modes_2, p + num_modes_2] += 1.0

        h_1 *= 0.25
        h_2 *= 0.25
        return h_1, h_2

    def total_angular_momentum(self):
        """
        Total angular momentum.

        A data_preprocess_helper fermionic operator which can be used to evaluate the total
        angular momentum of the given eigenstate.

        Returns:
            FermionicOperator: Fermionic Hamiltonian
        """
        x_h1, x_h2 = self._s_x_squared()
        y_h1, y_h2 = self._s_y_squared()
        z_h1, z_h2 = self._s_z_squared()
        h_1 = x_h1 + y_h1 + z_h1
        h_2 = x_h2 + y_h2 + z_h2

        return FermionicOperator(h1=h_1, h2=h_2)
