import cirq
import math
import numpy as np
from itertools import product
from rewriter.simple_circuit_rewrites import get_circuit_unitaries
# , get_random_state, apply_unitaries_to_circuit

def get_unitary_of_operation(op):
    circ = cirq.Circuit()
    circ.append(op)
    matrix = get_circuit_unitaries(circ)
    return matrix

def apply_unitaries_to_circuit(circuit, initial_state, dtype=np.complex128, noise_model=None):
    state = np.copy(initial_state)
    qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
    min_val = 100000
    for a in qd_map.keys():
        qd_map[a] = a.x
        min_val = min(min_val, a.x)
    for a in qd_map.keys():
        qd_map[a] = a.x - min_val
    buffer = np.zeros(state.shape, dtype=dtype)
    for m_index, moment in enumerate(circuit):
        for op in moment:
            qudit_indices = [qd_map[qd] for qd in op.qubits]
            mat = cirq.unitary(op)
            # get_unitary_of_operation(op)
            matrix = mat.astype(dtype).reshape((3,)*(2*len(op.qubits)))
            # print("compare the two operatiins ")
            # print(np.allclose(cirq.unitary(op).reshape((3,)*(2*len(op.qubits))), matrix))
            # print("matrix shape ", matrix.shape, matrix.size)
            # print("target ", qudit_indices)
            # print(qd_map)
            a = matrix.dot(matrix.conj().T)
            cirq.linalg.targeted_left_multiply(matrix, state, target_axes=qudit_indices, out=buffer)
            state, buffer = buffer, state
            # print(state.shape)
            # print("current operation")
            # print(op)
            if noise_model != None:
                # print("--- matrix size and shape", matrix.shape, matrix.size)
                if matrix.size == 9:
                    gate_error_matrix = noise_model.pick_single_qutrit_gate_channel()
                elif matrix.size == 81:
                    gate_error_matrix = noise_model.pick_two_qutrit_gate_channel()
                gate_error_matrix = gate_error_matrix.astype(dtype).reshape((3,) * (2*len(op.qubits)))
                cirq.linalg.targeted_left_multiply(gate_error_matrix, state, target_axes=qudit_indices, out=buffer)
                state, buffer = buffer, state

        if noise_model != None:
            for index in range(len(qd_map)):
                if any([len(op.qubits) == 2 for op in moment]):  # apply long idle channel
                    gate_error_matrix = noise_model.pick_long_idle_channel(state, buffer, index)
                else:
                    gate_error_matrix = noise_model.pick_short_idle_channel(state, buffer, index)
                gate_error_matrix = gate_error_matrix.astype(dtype).reshape((3,) * 2)
                cirq.linalg.targeted_left_multiply(gate_error_matrix, state, [index], out=buffer)
                buffer /= np.linalg.norm(buffer)  # idle errors may be incoherent (non-unitary) so need to renormalize
                state, buffer = buffer, state
    return state

X3 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
Z3 = np.array([[1, 0, 0], [0, math.e ** (math.pi * 1j * -2.0/3.0), 0], [0, 0, math.e ** (math.pi * 1j * -4.0/3.0)]])
Y3 = X3 @ Z3
V3 = X3 @ Z3 @ Z3


X2 = np.array([[0, 1], [1, 0]])
Z2 = np.array([[1, 0], [0, -1]])
Y2 = X2 @ Z2


def weighted_draw(weights):
    """Returns random index, with draw probabilities weighted by the given weights."""
    assert sum(weights) > .999 and sum(weights) < 1.001, 'sum is %s for weights: %s' % (sum(weights), weights)
    total = sum(weights)
    weights = [weight / total for weight in weights]  # normalize anyways to make sum exactly 1
    return np.random.choice(len(weights), p=weights)


class NoiseChannel(object):
    single_qutrit_kraus_operators = []
    single_qutrit_kraus_operator_weights = []
    two_qutrit_kraus_operators = []
    two_qutrit_kraus_operator_weights = []

    idle_channel_operators = []
    long_idle_channel_weights = []
    short_idle_channel_weights = []

    def pick_single_qutrit_gate_channel(self):
        index = weighted_draw(self.single_qutrit_kraus_operator_weights)
        return self.single_qutrit_kraus_operators[index]

    def pick_two_qutrit_gate_channel(self):
        index = weighted_draw(self.two_qutrit_kraus_operator_weights)
        return self.two_qutrit_kraus_operators[index]

    def pick_long_idle_channel(self, state, buffer, index):
        index = weighted_draw(self.long_idle_channel_weights)
        return self.idle_channel_operators[index]

    def pick_short_idle_channel(self, state, buffer, index):
        index = weighted_draw(self.short_idle_channel_weights)
        return self.idle_channel_operators[index]

class GenericQutritErrors(NoiseChannel):
    X3 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    Z3 = np.array([[1, 0, 0], [0, math.e ** (math.pi * 1j * -2.0/3.0), 0], [0, 0, math.e ** (math.pi * 1j * -4.0/3.0)]])
    Y3 = X3 @ Z3
    V3 = X3 @ Z3 @ Z3

    single_qutrit_kraus_operators = [np.eye(3), Z3, Z3 @ Z3, X3, X3 @ Z3, X3 @ Z3 @ Z3, X3 @ X3, X3 @ X3 @ Z3, X3 @ X3 @ Z3 @ Z3]
    two_qutrit_kraus_operators = [np.kron(np.eye(3), np.eye(3)), np.kron(np.eye(3), Z3), np.kron(np.eye(3), Z3 @ Z3), np.kron(np.eye(3), X3), np.kron(np.eye(3), X3 @ Z3), np.kron(np.eye(3), X3 @ Z3 @ Z3), np.kron(np.eye(3), X3 @ X3), np.kron(np.eye(3), X3 @ X3 @ Z3), np.kron(np.eye(3), X3 @ X3 @ Z3 @ Z3),
            np.kron(Z3, np.eye(3)), np.kron(Z3, Z3), np.kron(Z3, Z3 @ Z3), np.kron(Z3, X3), np.kron(Z3, X3 @ Z3), np.kron(Z3, X3 @ Z3 @ Z3), np.kron(Z3, X3 @ X3), np.kron(Z3, X3 @ X3 @ Z3), np.kron(Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(Z3 @ Z3, np.eye(3)), np.kron(Z3 @ Z3, Z3), np.kron(Z3 @ Z3, Z3 @ Z3), np.kron(Z3 @ Z3, X3), np.kron(Z3 @ Z3, X3 @ Z3), np.kron(Z3 @ Z3, X3 @ Z3 @ Z3), np.kron(Z3 @ Z3, X3 @ X3), np.kron(Z3 @ Z3, X3 @ X3 @ Z3), np.kron(Z3 @ Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3, np.eye(3)), np.kron(X3, Z3), np.kron(X3, Z3 @ Z3), np.kron(X3, X3), np.kron(X3, X3 @ Z3), np.kron(X3, X3 @ Z3 @ Z3), np.kron(X3, X3 @ X3), np.kron(X3, X3 @ X3 @ Z3), np.kron(X3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ Z3, np.eye(3)), np.kron(X3 @ Z3, Z3), np.kron(X3 @ Z3, Z3 @ Z3), np.kron(X3 @ Z3, X3), np.kron(X3 @ Z3, X3 @ Z3), np.kron(X3 @ Z3, X3 @ Z3 @ Z3), np.kron(X3 @ Z3, X3 @ X3), np.kron(X3 @ Z3, X3 @ X3 @ Z3), np.kron(X3 @ Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ Z3 @ Z3, np.eye(3)), np.kron(X3 @ Z3 @ Z3, Z3), np.kron(X3 @ Z3 @ Z3, Z3 @ Z3), np.kron(X3 @ Z3 @ Z3, X3), np.kron(X3 @ Z3 @ Z3, X3 @ Z3), np.kron(X3 @ Z3 @ Z3, X3 @ Z3 @ Z3), np.kron(X3 @ Z3 @ Z3, X3 @ X3), np.kron(X3 @ Z3 @ Z3, X3 @ X3 @ Z3), np.kron(X3 @ Z3 @ Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ X3, np.eye(3)), np.kron(X3 @ X3, Z3), np.kron(X3 @ X3, Z3 @ Z3), np.kron(X3 @ X3, X3), np.kron(X3 @ X3, X3 @ Z3), np.kron(X3 @ X3, X3 @ Z3 @ Z3), np.kron(X3 @ X3, X3 @ X3), np.kron(X3 @ X3, X3 @ X3 @ Z3), np.kron(X3 @ X3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ X3 @ Z3, np.eye(3)), np.kron(X3 @ X3 @ Z3, Z3), np.kron(X3 @ X3 @ Z3, Z3 @ Z3), np.kron(X3 @ X3 @ Z3, X3), np.kron(X3 @ X3 @ Z3, X3 @ Z3), np.kron(X3 @ X3 @ Z3, X3 @ Z3 @ Z3), np.kron(X3 @ X3 @ Z3, X3 @ X3), np.kron(X3 @ X3 @ Z3, X3 @ X3 @ Z3), np.kron(X3 @ X3 @ Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ X3 @ Z3 @ Z3, np.eye(3)), np.kron(X3 @ X3 @ Z3 @ Z3, Z3), np.kron(X3 @ X3 @ Z3 @ Z3, Z3 @ Z3), np.kron(X3 @ X3 @ Z3 @ Z3, X3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ Z3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ Z3 @ Z3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ X3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ X3 @ Z3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ X3 @ Z3 @ Z3)]

    def pick_long_idle_channel(self, state, buffer, index):
        weights = []
        #print("QUTRIT LONG ", self.long_idle_channel_operators)
        for gate_error_matrix in self.long_idle_channel_operators:
            gate_error_matrix = gate_error_matrix.astype(np.complex128).reshape((3,) * 2)
            cirq.linalg.targeted_left_multiply(gate_error_matrix, state[:], [index], out=buffer)
            weights.append(np.linalg.norm(buffer) ** 2)

        index = weighted_draw(weights)
        return self.short_idle_channel_operators[index]

    def pick_short_idle_channel(self, state, buffer, index):
        weights = []
        #print("QUTRIT short ", self.short_idle_channel_operators)
        for gate_error_matrix in self.short_idle_channel_operators:
            gate_error_matrixn = gate_error_matrix.astype(np.complex128).reshape((3,) * 2)
            cirq.linalg.targeted_left_multiply(gate_error_matrix, state[:], [index], out=buffer)
            weights.append(np.linalg.norm(buffer) ** 2)

        index = weighted_draw(weights)
        return self.long_idle_channel_operators[index]



class GenericQubitErrors(NoiseChannel):

    single_qutrit_kraus_operators = [np.eye(2), Z2, X2, X2 @ Z2]
    two_qutrit_kraus_operators = [np.kron(np.eye(2), np.eye(2)), np.kron(np.eye(2), Z2), np.kron(np.eye(2), X2 @ Z2), np.kron(np.eye(2), X2),
            np.kron(Z2, np.eye(2)), np.kron(Z2, Z2), np.kron(Z2, X2 @ Z2), np.kron(Z2, X2),
            np.kron(X2 @ Z2, np.eye(2)), np.kron(X2 @ Z2, X2), np.kron(X2 @ Z2, X2 @ Z2), np.kron(X2 @ Z2, Z2),
            np.kron(X2, np.eye(2)), np.kron(X2, Z2), np.kron(X2, X2 @ Z2), np.kron(X2, X2)]

    def pick_long_idle_channel(self, state, buffer, index):
        weights = []
        #print("QUBIT LONG ", self.long_idle_channel_operators)
        for gate_error_matrix in self.long_idle_channel_operators:
            gate_error_matrix = gate_error_matrix.astype(np.complex128).reshape((2,) * 2)
            cirq.linalg.targeted_left_multiply(gate_error_matrix, state[:], [index], out=buffer)
            weights.append(np.linalg.norm(buffer) ** 2)

        index = weighted_draw(weights)
        return self.short_idle_channel_operators[index]

    def pick_short_idle_channel(self, state, buffer, index):
        weights = []
        #print("QUBIT SHORT ", self.short_idle_channel_operators)
        for gate_error_matrix in self.short_idle_channel_operators:
            gate_error_matrixn = gate_error_matrix.astype(np.complex128).reshape((2,) * 2)
            cirq.linalg.targeted_left_multiply(gate_error_matrix, state[:], [index], out=buffer)
            weights.append(np.linalg.norm(buffer) ** 2)

        index = weighted_draw(weights)
        return self.long_idle_channel_operators[index]


class AssymetricQutritErrors(NoiseChannel):
    X3 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    Z3 = np.array([[1, 0, 0], [0, math.e ** (math.pi * 1j * -2.0/3.0), 0], [0, 0, math.e ** (math.pi * 1j * -4.0/3.0)]])
    Y3 = X3 @ Z3
    V3 = X3 @ Z3 @ Z3

    single_qutrit_kraus_operators = [np.eye(3), Z3, Z3 @ Z3, X3, X3 @ Z3, X3 @ Z3 @ Z3, X3 @ X3, X3 @ X3 @ Z3, X3 @ X3 @ Z3 @ Z3]
    two_qutrit_kraus_operators = [np.kron(np.eye(3), np.eye(3)), np.kron(np.eye(3), Z3), np.kron(np.eye(3), Z3 @ Z3), np.kron(np.eye(3), X3), np.kron(np.eye(3), X3 @ Z3), np.kron(np.eye(3), X3 @ Z3 @ Z3), np.kron(np.eye(3), X3 @ X3), np.kron(np.eye(3), X3 @ X3 @ Z3), np.kron(np.eye(3), X3 @ X3 @ Z3 @ Z3),
            np.kron(Z3, np.eye(3)), np.kron(Z3, Z3), np.kron(Z3, Z3 @ Z3), np.kron(Z3, X3), np.kron(Z3, X3 @ Z3), np.kron(Z3, X3 @ Z3 @ Z3), np.kron(Z3, X3 @ X3), np.kron(Z3, X3 @ X3 @ Z3), np.kron(Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(Z3 @ Z3, np.eye(3)), np.kron(Z3 @ Z3, Z3), np.kron(Z3 @ Z3, Z3 @ Z3), np.kron(Z3 @ Z3, X3), np.kron(Z3 @ Z3, X3 @ Z3), np.kron(Z3 @ Z3, X3 @ Z3 @ Z3), np.kron(Z3 @ Z3, X3 @ X3), np.kron(Z3 @ Z3, X3 @ X3 @ Z3), np.kron(Z3 @ Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3, np.eye(3)), np.kron(X3, Z3), np.kron(X3, Z3 @ Z3), np.kron(X3, X3), np.kron(X3, X3 @ Z3), np.kron(X3, X3 @ Z3 @ Z3), np.kron(X3, X3 @ X3), np.kron(X3, X3 @ X3 @ Z3), np.kron(X3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ Z3, np.eye(3)), np.kron(X3 @ Z3, Z3), np.kron(X3 @ Z3, Z3 @ Z3), np.kron(X3 @ Z3, X3), np.kron(X3 @ Z3, X3 @ Z3), np.kron(X3 @ Z3, X3 @ Z3 @ Z3), np.kron(X3 @ Z3, X3 @ X3), np.kron(X3 @ Z3, X3 @ X3 @ Z3), np.kron(X3 @ Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ Z3 @ Z3, np.eye(3)), np.kron(X3 @ Z3 @ Z3, Z3), np.kron(X3 @ Z3 @ Z3, Z3 @ Z3), np.kron(X3 @ Z3 @ Z3, X3), np.kron(X3 @ Z3 @ Z3, X3 @ Z3), np.kron(X3 @ Z3 @ Z3, X3 @ Z3 @ Z3), np.kron(X3 @ Z3 @ Z3, X3 @ X3), np.kron(X3 @ Z3 @ Z3, X3 @ X3 @ Z3), np.kron(X3 @ Z3 @ Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ X3, np.eye(3)), np.kron(X3 @ X3, Z3), np.kron(X3 @ X3, Z3 @ Z3), np.kron(X3 @ X3, X3), np.kron(X3 @ X3, X3 @ Z3), np.kron(X3 @ X3, X3 @ Z3 @ Z3), np.kron(X3 @ X3, X3 @ X3), np.kron(X3 @ X3, X3 @ X3 @ Z3), np.kron(X3 @ X3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ X3 @ Z3, np.eye(3)), np.kron(X3 @ X3 @ Z3, Z3), np.kron(X3 @ X3 @ Z3, Z3 @ Z3), np.kron(X3 @ X3 @ Z3, X3), np.kron(X3 @ X3 @ Z3, X3 @ Z3), np.kron(X3 @ X3 @ Z3, X3 @ Z3 @ Z3), np.kron(X3 @ X3 @ Z3, X3 @ X3), np.kron(X3 @ X3 @ Z3, X3 @ X3 @ Z3), np.kron(X3 @ X3 @ Z3, X3 @ X3 @ Z3 @ Z3),
            np.kron(X3 @ X3 @ Z3 @ Z3, np.eye(3)), np.kron(X3 @ X3 @ Z3 @ Z3, Z3), np.kron(X3 @ X3 @ Z3 @ Z3, Z3 @ Z3), np.kron(X3 @ X3 @ Z3 @ Z3, X3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ Z3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ Z3 @ Z3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ X3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ X3 @ Z3), np.kron(X3 @ X3 @ Z3 @ Z3, X3 @ X3 @ Z3 @ Z3)]

    def pick_long_idle_channel(self, state, buffer, index):
        weights = []
        #print("QUTRIT LONG ", self.long_idle_channel_operators)
        for gate_error_matrix in self.long_idle_channel_operators:
            gate_error_matrix = gate_error_matrix.astype(np.complex128).reshape((3,) * 2)
            cirq.linalg.targeted_left_multiply(gate_error_matrix, state[:], [index], out=buffer)
            weights.append(np.linalg.norm(buffer) ** 2)

        index = weighted_draw(weights)
        return self.short_idle_channel_operators[index]

    def pick_short_idle_channel(self, state, buffer, index):
        weights = []
        #print("QUTRIT short ", self.short_idle_channel_operators)
        for gate_error_matrix in self.short_idle_channel_operators:
            gate_error_matrixn = gate_error_matrix.astype(np.complex128).reshape((3,) * 2)
            cirq.linalg.targeted_left_multiply(gate_error_matrix, state[:], [index], out=buffer)
            weights.append(np.linalg.norm(buffer) ** 2)

        index = weighted_draw(weights)

class CurrentSuperconductingQCErrors(GenericQutritErrors):
    p_1 = .001 / 3  # single qubit gate error  (https://arxiv.org/pdf/1702.01852.pdf and www.research.ibm.com/ibm-q/technology/devices)
                    # division by 3 because 3 possible error channels for qubits
    p_2 = .01 / 15  # same references, division by 15 because 15 error channels for qubits
    single_qutrit_kraus_operator_weights = [1 - 8*p_1] + 8 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 80*p_2] + 80 * [p_2]
    # gamma_m = 1 - e^(-m * dt / T1), where dt is gate time duration (https://web.physics.ucsb.edu/~martinisgroup/papers/Ghosh2013b.pdf)
    # Here, I picked T1 = 100 microseconds (representative from https://www.research.ibm.com/ibm-q/technology/devices/)
    # and dt = 100 ns for single-qudit gates and 300 ns for two-qudit gates (https://arxiv.org/pdf/1702.01852.pdf)
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 100000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 100000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 100000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 100000.0)
    short_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_short)**.5,0],[0,0,(1-gamma_2_short)**.5]]), np.array([[0,gamma_1_short**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_short**0.5],[0,0,0],[0,0,0]])]
    long_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_long)**.5,0],[0,0,(1-gamma_2_long)**.5]]), np.array([[0,gamma_1_long**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_long**0.5],[0,0,0],[0,0,0]])]


class ModifiedCurrSuperconductingQCErrors(GenericQutritErrors):
    p_1 = .001 / 8  # single qubit gate error  (https://arxiv.org/pdf/1702.01852.pdf and www.research.ibm.com/ibm-q/technology/devices)
                    # division by 3 because 3 possible error channels for qubits
    p_2 = .01 / 80  # same references, division by 15 because 15 error channels for qubits
    single_qutrit_kraus_operator_weights = [1 - 8*p_1] + 8 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 80*p_2] + 80 * [p_2]
    # gamma_m = 1 - e^(-m * dt / T1), where dt is gate time duration (https://web.physics.ucsb.edu/~martinisgroup/papers/Ghosh2013b.pdf)
    # Here, I picked T1 = 100 microseconds (representative from https://www.research.ibm.com/ibm-q/technology/devices/)
    # and dt = 100 ns for single-qudit gates and 300 ns for two-qudit gates (https://arxiv.org/pdf/1702.01852.pdf)
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 100000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 100000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 100000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 100000.0)
    short_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_short)**.5,0],[0,0,(1-gamma_2_short)**.5]]), np.array([[0,gamma_1_short**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_short**0.5],[0,0,0],[0,0,0]])]
    long_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_long)**.5,0],[0,0,(1-gamma_2_long)**.5]]), np.array([[0,gamma_1_long**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_long**0.5],[0,0,0],[0,0,0]])]

class ModifiedFutureSuperconductingQCErrors(GenericQutritErrors):
    p_1 = .0001 / 8  # single qubit gate error  (https://arxiv.org/pdf/1702.01852.pdf and www.research.ibm.com/ibm-q/technology/devices)
                    # division by 3 because 3 possible error channels for qubits
    p_2 = .001 / 80  # same references, division by 15 because 15 error channels for qubits

    single_qutrit_kraus_operator_weights = [1 - 8*p_1] + 8 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 80*p_2] + 80 * [p_2]
    # gamma_m = 1 - e^(-m * dt / T1), where dt is gate time duration (https://web.physics.ucsb.edu/~martinisgroup/papers/Ghosh2013b.pdf)
    # Here, I picked T1 = 100 microseconds (representative from https://www.research.ibm.com/ibm-q/technology/devices/)
    # and dt = 100 ns for single-qudit gates and 300 ns for two-qudit gates (https://arxiv.org/pdf/1702.01852.pdf)
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 100000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 100000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 100000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 100000.0)
    short_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_short)**.5,0],[0,0,(1-gamma_2_short)**.5]]), np.array([[0,gamma_1_short**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_short**0.5],[0,0,0],[0,0,0]])]
    long_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_long)**.5,0],[0,0,(1-gamma_2_long)**.5]]), np.array([[0,gamma_1_long**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_long**0.5],[0,0,0],[0,0,0]])]


class CurrentSuperconductingQubitErrors(GenericQubitErrors):
    p_1 = .001 / 3  # single qubit gate error  (https://arxiv.org/pdf/1702.01852.pdf and www.research.ibm.com/ibm-q/technology/devices)
                    # division by 3 because 3 possible error channels for qubits
    p_2 = .01 / 15  # same references, division by 15 because 15 error channels for qubits
    single_qutrit_kraus_operator_weights = [1 - 3*p_1] + 3 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 15*p_2] + 15 * [p_2]
    # gamma_m = 1 - e^(-m * dt / T1), where dt is gate time duration (https://web.physics.ucsb.edu/~martinisgroup/papers/Ghosh2013b.pdf)
    # Here, I picked T1 = 100 microseconds (representative from https://www.research.ibm.com/ibm-q/technology/devices/)
    # and dt = 100 ns for single-qudit gates and 300 ns for two-qudit gates (https://arxiv.org/pdf/1702.01852.pdf)
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 100000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 100000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 100000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 100000.0)
    short_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_short)**.5]]), np.array([[0,gamma_1_short**0.5],[0,0]])]
    long_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_long)**.5]]), np.array([[0,gamma_1_long**0.5],[0,0]])]
    #long_idle_channel_operators = [np.array([[1,0],[0,0]]), np.array([[0,1],[0,0]])]


class FutureSuperconductingQCErrors(GenericQutritErrors):
    # Here, I 10x'ed T1 and 1/10x'ed the p's
    p_1 = .0001 / 3
    p_2 = .001 / 15
    single_qutrit_kraus_operator_weights = [1 - 8*p_1] + 8 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 80*p_2] + 80 * [p_2]
    # Here, I picked T1 = 1000 microseconds 
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 1000000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 1000000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 1000000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 1000000.0)
    short_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_short)**.5,0],[0,0,(1-gamma_2_short)**.5]]), np.array([[0,gamma_1_short**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_short**0.5],[0,0,0],[0,0,0]])]
    long_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_long)**.5,0],[0,0,(1-gamma_2_long)**.5]]), np.array([[0,gamma_1_long**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_long**0.5],[0,0,0],[0,0,0]])]


class FutureSuperconductingQubitErrors(GenericQubitErrors):
    # Here, I 10x'ed T1 and 1/10x'ed the p's
    p_1 = .0001 / 3
    p_2 = .001 / 15
    single_qutrit_kraus_operator_weights = [1 - 3*p_1] + 3 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 15*p_2] + 15 * [p_2]
    # Here, I picked T1 = 1000 microseconds 
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 1000000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 1000000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 1000000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 1000000.0)
    short_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_short)**.5]]), np.array([[0,gamma_1_short**0.5,],[0,0]])]
    long_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_long)**.5]]), np.array([[0,gamma_1_long**0.5],[0,0]])]


class FutureSuperconductingQCErrorsBetterT1(GenericQutritErrors):
    p_1 = .0001 / 3
    p_2 = .001 / 15
    single_qutrit_kraus_operator_weights = [1 - 8*p_1] + 8 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 80*p_2] + 80 * [p_2]
    # Here, I picked T1 = 10000 microseconds 
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 10000000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 10000000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 10000000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 10000000.0)
    short_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_short)**.5,0],[0,0,(1-gamma_2_short)**.5]]), np.array([[0,gamma_1_short**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_short**0.5],[0,0,0],[0,0,0]])]
    long_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_long)**.5,0],[0,0,(1-gamma_2_long)**.5]]), np.array([[0,gamma_1_long**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_long**0.5],[0,0,0],[0,0,0]])]


class FutureSuperconductingQCErrorsBetterGates(GenericQutritErrors):
    p_1 = .00001 / 3
    p_2 = .0001 / 15
    single_qutrit_kraus_operator_weights = [1 - 8*p_1] + 8 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 80*p_2] + 80 * [p_2]
    # Here, I picked T1 = 1000 microseconds 
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 1000000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 1000000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 1000000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 1000000.0)
    short_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_short)**.5,0],[0,0,(1-gamma_2_short)**.5]]), np.array([[0,gamma_1_short**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_short**0.5],[0,0,0],[0,0,0]])]
    long_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_long)**.5,0],[0,0,(1-gamma_2_long)**.5]]), np.array([[0,gamma_1_long**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_long**0.5],[0,0,0],[0,0,0]])]


class FutureSuperconductingQCErrorsBetterT1AndGates(GenericQutritErrors):
    p_1 = .00001 / 3
    p_2 = .0001 / 15
    single_qutrit_kraus_operator_weights = [1 - 8*p_1] + 8 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 80*p_2] + 80 * [p_2]
    # Here, I picked T1 = 10000 microseconds 
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 10000000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 10000000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 10000000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 10000000.0)
    short_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_short)**.5,0],[0,0,(1-gamma_2_short)**.5]]), np.array([[0,gamma_1_short**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_short**0.5],[0,0,0],[0,0,0]])]
    long_idle_channel_operators = [np.array([[1,0,0],[0,(1-gamma_1_long)**.5,0],[0,0,(1-gamma_2_long)**.5]]), np.array([[0,gamma_1_long**0.5,0],[0,0,0],[0,0,0]]), np.array([[0,0,gamma_2_long**0.5],[0,0,0],[0,0,0]])]


class FutureSuperconductingQubitErrorsBetterT1(GenericQubitErrors):
    p_1 = .0001 / 3
    p_2 = .001 / 15
    single_qutrit_kraus_operator_weights = [1 - 3*p_1] + 3 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 15*p_2] + 15 * [p_2]
    # Here, I picked T1 = 10000 microseconds 
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 10000000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 10000000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 10000000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 10000000.0)
    short_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_short)**.5]]), np.array([[0,gamma_1_short**0.5],[0,0]])]
    long_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_long)**.5]]), np.array([[0,gamma_1_long**0.5],[0,0]])]


class FutureSuperconductingQubitErrorsBetterGates(GenericQubitErrors):
    p_1 = .00001 / 3
    p_2 = .0001 / 15
    single_qutrit_kraus_operator_weights = [1 - 3*p_1] + 3 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 15*p_2] + 15 * [p_2]
    # Here, I picked T1 = 1000 microseconds 
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 1000000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 1000000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 1000000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 1000000.0)
    short_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_short)**.5]]), np.array([[0,gamma_1_short**0.5],[0,0]])]
    long_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_long)**.5]]), np.array([[0,gamma_1_long**0.5],[0,0]])]


class FutureSuperconductingQubitErrorsBetterT1AndGates(GenericQubitErrors):
    p_1 = .00001 / 3
    p_2 = .0001 / 15
    single_qutrit_kraus_operator_weights = [1 - 3*p_1] + 3 * [p_1]
    two_qutrit_kraus_operator_weights = [1 - 15*p_2] + 15 * [p_2]
    # Here, I picked T1 = 10000 microseconds 
    gamma_1_short = 1 - math.exp(-1 *  100.0 / 10000000.0)
    gamma_1_long = 1 - math.exp(-1 * 300.0 / 10000000.0)
    gamma_2_short = 1 - math.exp(-2 *  100.0 / 10000000.0)
    gamma_2_long = 1 - math.exp(-2 * 300.0 / 10000000.0)
    short_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_short)**.5]]), np.array([[0,gamma_1_short**0.5],[0,0]])]
    long_idle_channel_operators = [np.array([[1,0],[0,(1-gamma_1_long)**.5]]), np.array([[0,gamma_1_long**0.5],[0,0]])]


class DressedQutritErrors(NoiseChannel):
    idle_channel_operators = [Z3, Z3 @ Z3, np.eye(3)]
    long_idle_channel_weights = [0.0000228617, 0.0000228617, 0.9999542766]
    short_idle_channel_weights = [5.71579E-10, 5.71579E-10, 0.99999999885]

    single_qutrit_kraus_operators = [X3, X3 @ X3, Z3, Z3 @ Z3, Y3, Y3 @ Y3, V3, V3 @ V3, np.eye(3)]
    single_qutrit_kraus_operator_weights = [8.07695E-7, 8.35653E-7, 4.49109E-7, 4.49109E-7, 8.07695E-7, 8.35653E-7, 8.07695E-7, 8.35653E-7, 0.999846624680123]

    XX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    X2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    ZX_weights = [2.512069977868293e-10,2.5990231483578203e-10,1.3968041859275328e-10,1.3968041859275328e-10,2.512069977868293e-10,2.5990231483578203e-10,2.512069977868293e-10,2.5990231483578203e-10,0.000011816834382389788]
    Z2X_weights = [2.512069977868293e-10,2.5990231483578203e-10,1.3968041859275328e-10,1.3968041859275328e-10,2.512069977868293e-10,2.5990231483578203e-10,2.512069977868293e-10,2.5990231483578203e-10,0.000011816834382389788]
    YX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    Y2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    VX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    V2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    TX_weights = [0.000021251879958915915,0.00002198749574891023,0.000011816834382389788,0.000011816834382389788,0.000021251879958915915,0.00002198749574891023,0.000021251879958915915,0.00002198749574891023,0.9996932728842347]

    XX_operators = [np.kron(X3, X3), np.kron(X3, X3 @ X3), np.kron(X3, Z3), np.kron(X3, Z3 @ Z3), np.kron(X3, Y3), np.kron(X3, Y3 @ Y3), np.kron(X3, V3), np.kron(X3, V3 @ V3), np.kron(X3, np.eye(3))]
    X2X_operators = [np.kron(X3 @ X3, X3), np.kron(X3 @ X3, X3 @ X3), np.kron(X3 @ X3, Z3), np.kron(X3 @ X3, Z3 @ Z3), np.kron(X3 @ X3, Y3), np.kron(X3 @ X3, Y3 @ Y3), np.kron(X3 @ X3, V3), np.kron(X3 @ X3, V3 @ V3), np.kron(X3 @ X3, np.eye(3))]
    ZX_operators = [np.kron(Z3, X3), np.kron(Z3, X3 @ X3), np.kron(Z3, Z3), np.kron(Z3, Z3 @ Z3), np.kron(Z3, Y3), np.kron(Z3, Y3 @ Y3), np.kron(Z3, V3), np.kron(Z3, V3 @ V3), np.kron(Z3, np.eye(3))]
    Z2X_operators = [np.kron(Z3 @ Z3, X3), np.kron(Z3 @ Z3, X3 @ X3), np.kron(Z3 @ Z3, Z3), np.kron(Z3 @ Z3, Z3 @ Z3), np.kron(Z3 @ Z3, Y3), np.kron(Z3 @ Z3, Y3 @ Y3), np.kron(Z3 @ Z3, V3), np.kron(Z3 @ Z3, V3 @ V3), np.kron(Z3 @ Z3, np.eye(3))]
    YX_operators = [np.kron(Y3, X3), np.kron(Y3, X3 @ X3), np.kron(Y3, Z3), np.kron(Y3, Z3 @ Z3), np.kron(Y3, Y3), np.kron(Y3, Y3 @ Y3), np.kron(Y3, V3), np.kron(Y3, V3 @ V3), np.kron(Y3, np.eye(3))]
    Y2X_operators = [np.kron(Y3 @ Y3, X3), np.kron(Y3 @ Y3, X3 @ X3), np.kron(Y3 @ Y3, Z3), np.kron(Y3 @ Y3, Z3 @ Z3), np.kron(Y3 @ Y3, Y3), np.kron(Y3 @ Y3, Y3 @ Y3), np.kron(Y3 @ Y3, V3), np.kron(Y3 @ Y3, V3 @ V3), np.kron(Y3 @ Y3, np.eye(3))]
    VX_operators = [np.kron(V3, X3), np.kron(V3, X3 @ X3), np.kron(V3, Z3), np.kron(V3, Z3 @ Z3), np.kron(V3, Y3), np.kron(V3, Y3 @ Y3), np.kron(V3, V3), np.kron(V3, V3 @ V3), np.kron(V3, np.eye(3))]
    V2X_operators = [np.kron(V3 @ V3, X3), np.kron(V3 @ V3, X3 @ X3), np.kron(V3 @ V3, Z3), np.kron(V3 @ V3, Z3 @ Z3), np.kron(V3 @ V3, Y3), np.kron(V3 @ V3, Y3 @ Y3), np.kron(V3 @ V3, V3), np.kron(V3 @ V3, V3 @ V3), np.kron(V3 @ V3, np.eye(3))]
    TX_operators = [np.kron(np.eye(3), X3), np.kron(np.eye(3), X3 @ X3), np.kron(np.eye(3), Z3), np.kron(np.eye(3), Z3 @ Z3), np.kron(np.eye(3), Y3), np.kron(np.eye(3), Y3 @ Y3), np.kron(np.eye(3), V3), np.kron(np.eye(3), V3 @ V3), np.kron(np.eye(3), np.eye(3))]

    two_qutrit_kraus_operators = XX_operators + X2X_operators + ZX_operators + Z2X_operators + YX_operators + Y2X_operators + VX_operators + V2X_operators + TX_operators
    two_qutrit_kraus_operator_weights =  XX_weights + X2X_weights + ZX_weights + Z2X_weights + YX_weights + Y2X_weights + VX_weights + V2X_weights + TX_weights


class BareQutritErrors(NoiseChannel):
    idle_channel_operators = [Z3, Z3 @ Z3, np.eye(3)]
    long_idle_channel_weights = [0, 0.00007751, 0.99992249] 
    short_idle_channel_weights = [0, 0.0000000019379, 0.99999999806]

    single_qutrit_kraus_operators = [X3, X3 @ X3, Z3, Z3 @ Z3, Y3, Y3 @ Y3, V3, V3 @ V3, np.eye(3)]
    single_qutrit_kraus_operator_weights = [8.076953189667142e-7,8.356530070058442e-7,4.4910858870426395e-7,4.4910858870426395e-7,8.076953189667142e-7,8.356530070058442e-7,8.076953189667142e-7,8.356530070058442e-7,0.999846624680123]

    XX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    X2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    ZX_weights = [2.512069977868293e-10,2.5990231483578203e-10,1.3968041859275328e-10,1.3968041859275328e-10,2.512069977868293e-10,2.5990231483578203e-10,2.512069977868293e-10,2.5990231483578203e-10,0.000011816834382389788]
    Z2X_weights = [2.512069977868293e-10,2.5990231483578203e-10,1.3968041859275328e-10,1.3968041859275328e-10,2.512069977868293e-10,2.5990231483578203e-10,2.512069977868293e-10,2.5990231483578203e-10,0.000011816834382389788]
    YX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    Y2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    VX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    V2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    TX_weights = [0.000021251879958915915,0.00002198749574891023,0.000011816834382389788,0.000011816834382389788,0.000021251879958915915,0.00002198749574891023,0.000021251879958915915,0.00002198749574891023,0.9996932728842347]

    XX_operators = [np.kron(X3, X3), np.kron(X3, X3 @ X3), np.kron(X3, Z3), np.kron(X3, Z3 @ Z3), np.kron(X3, Y3), np.kron(X3, Y3 @ Y3), np.kron(X3, V3), np.kron(X3, V3 @ V3), np.kron(X3, np.eye(3))]
    X2X_operators = [np.kron(X3 @ X3, X3), np.kron(X3 @ X3, X3 @ X3), np.kron(X3 @ X3, Z3), np.kron(X3 @ X3, Z3 @ Z3), np.kron(X3 @ X3, Y3), np.kron(X3 @ X3, Y3 @ Y3), np.kron(X3 @ X3, V3), np.kron(X3 @ X3, V3 @ V3), np.kron(X3 @ X3, np.eye(3))]
    ZX_operators = [np.kron(Z3, X3), np.kron(Z3, X3 @ X3), np.kron(Z3, Z3), np.kron(Z3, Z3 @ Z3), np.kron(Z3, Y3), np.kron(Z3, Y3 @ Y3), np.kron(Z3, V3), np.kron(Z3, V3 @ V3), np.kron(Z3, np.eye(3))]
    Z2X_operators = [np.kron(Z3 @ Z3, X3), np.kron(Z3 @ Z3, X3 @ X3), np.kron(Z3 @ Z3, Z3), np.kron(Z3 @ Z3, Z3 @ Z3), np.kron(Z3 @ Z3, Y3), np.kron(Z3 @ Z3, Y3 @ Y3), np.kron(Z3 @ Z3, V3), np.kron(Z3 @ Z3, V3 @ V3), np.kron(Z3 @ Z3, np.eye(3))]
    YX_operators = [np.kron(Y3, X3), np.kron(Y3, X3 @ X3), np.kron(Y3, Z3), np.kron(Y3, Z3 @ Z3), np.kron(Y3, Y3), np.kron(Y3, Y3 @ Y3), np.kron(Y3, V3), np.kron(Y3, V3 @ V3), np.kron(Y3, np.eye(3))]
    Y2X_operators = [np.kron(Y3 @ Y3, X3), np.kron(Y3 @ Y3, X3 @ X3), np.kron(Y3 @ Y3, Z3), np.kron(Y3 @ Y3, Z3 @ Z3), np.kron(Y3 @ Y3, Y3), np.kron(Y3 @ Y3, Y3 @ Y3), np.kron(Y3 @ Y3, V3), np.kron(Y3 @ Y3, V3 @ V3), np.kron(Y3 @ Y3, np.eye(3))]
    VX_operators = [np.kron(V3, X3), np.kron(V3, X3 @ X3), np.kron(V3, Z3), np.kron(V3, Z3 @ Z3), np.kron(V3, Y3), np.kron(V3, Y3 @ Y3), np.kron(V3, V3), np.kron(V3, V3 @ V3), np.kron(V3, np.eye(3))]
    V2X_operators = [np.kron(V3 @ V3, X3), np.kron(V3 @ V3, X3 @ X3), np.kron(V3 @ V3, Z3), np.kron(V3 @ V3, Z3 @ Z3), np.kron(V3 @ V3, Y3), np.kron(V3 @ V3, Y3 @ Y3), np.kron(V3 @ V3, V3), np.kron(V3 @ V3, V3 @ V3), np.kron(V3 @ V3, np.eye(3))]
    TX_operators = [np.kron(np.eye(3), X3), np.kron(np.eye(3), X3 @ X3), np.kron(np.eye(3), Z3), np.kron(np.eye(3), Z3 @ Z3), np.kron(np.eye(3), Y3), np.kron(np.eye(3), Y3 @ Y3), np.kron(np.eye(3), V3), np.kron(np.eye(3), V3 @ V3), np.kron(np.eye(3), np.eye(3))]

    two_qutrit_kraus_operators = XX_operators + X2X_operators + ZX_operators + Z2X_operators + YX_operators + Y2X_operators + VX_operators + V2X_operators + TX_operators
    two_qutrit_kraus_operator_weights =  XX_weights + X2X_weights + ZX_weights + Z2X_weights + YX_weights + Y2X_weights + VX_weights + V2X_weights + TX_weights


class DressedQubitErrors(NoiseChannel):
    idle_channel_operators = [Z2, np.eye(2)]
    long_idle_channel_weights = [0.0000228617, 0.9999771383]
    short_idle_channel_weights = [5.71579E-10, 0.99999999942]

    single_qutrit_kraus_operators = [X2, Z2, Y2, np.eye(2)]
    single_qutrit_kraus_operator_weights = [8.07695E-7, 4.49109E-7, 8.07695E-7,  0.999979355]


    #XX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    XX_weights = [4.517809752636722e-10,2.512069977868293e-10,4.517809752636722e-10,0.000021251879958915915]

    #X2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    ZX_weights = [2.512069977868293e-10,1.3968041859275328e-10,2.512069977868293e-10,0.000011816834382389788]
    
    #Z2X_weights = [2.512069977868293e-10,2.5990231483578203e-10,1.3968041859275328e-10,1.3968041859275328e-10,2.512069977868293e-10,2.5990231483578203e-10,2.512069977868293e-10,2.5990231483578203e-10,0.000011816834382389788]
    YX_weights = [4.517809752636722e-10,2.512069977868293e-10,4.517809752636722e-10,0.000021251879958915915]
    TX_weights = [0.000021251879958915915,0.000011816834382389788,0.000021251879958915915, 0.99995871042] #   0.9996932728842347]

    XX_operators = [np.kron(X2, X2), np.kron(X2, Z2), np.kron(X2, Y2), np.kron(X2, np.eye(2))]
    ZX_operators = [np.kron(Z2, X2), np.kron(Z2, Z2), np.kron(Z2, Y2), np.kron(Z2, np.eye(2))]
    #Z2X_operators = [np.kron(Z3 @ Z3, X3), np.kron(Z3 @ Z3, X3 @ X3), np.kron(Z3 @ Z3, Z3), np.kron(Z3 @ Z3, Z3 @ Z3), np.kron(Z3 @ Z3, Y3), np.kron(Z3 @ Z3, Y3 @ Y3), np.kron(Z3 @ Z3, V3), np.kron(Z3 @ Z3, V3 @ V3), np.kron(Z3 @ Z3, np.eye(3))]
    YX_operators = [np.kron(Y2, X2), np.kron(Y2, Z2), np.kron(Y2, Y2), np.kron(Y2, np.eye(2))]
    #Y2X_operators = [np.kron(Y3 @ Y3, X3), np.kron(Y3 @ Y3, X3 @ X3), np.kron(Y3 @ Y3, Z3), np.kron(Y3 @ Y3, Z3 @ Z3), np.kron(Y3 @ Y3, Y3), np.kron(Y3 @ Y3, Y3 @ Y3), np.kron(Y3 @ Y3, V3), np.kron(Y3 @ Y3, V3 @ V3), np.kron(Y3 @ Y3, np.eye(3))]
    #VX_operators = [np.kron(V3, X3), np.kron(V3, X3 @ X3), np.kron(V3, Z3), np.kron(V3, Z3 @ Z3), np.kron(V3, Y3), np.kron(V3, Y3 @ Y3), np.kron(V3, V3), np.kron(V3, V3 @ V3), np.kron(V3, np.eye(3))]
    #V2X_operators = [np.kron(V3 @ V3, X3), np.kron(V3 @ V3, X3 @ X3), np.kron(V3 @ V3, Z3), np.kron(V3 @ V3, Z3 @ Z3), np.kron(V3 @ V3, Y3), np.kron(V3 @ V3, Y3 @ Y3), np.kron(V3 @ V3, V3), np.kron(V3 @ V3, V3 @ V3), np.kron(V3 @ V3, np.eye(3))]
    TX_operators = [np.kron(np.eye(2), X2), np.kron(np.eye(2),Z2), np.kron(np.eye(2), Y2), np.kron(np.eye(2), np.eye(2))]

    two_qutrit_kraus_operators = XX_operators + ZX_operators + YX_operators + TX_operators
    two_qutrit_kraus_operator_weights =  XX_weights + ZX_weights + YX_weights + TX_weights


class BareQubitErrors(NoiseChannel):
    idle_channel_operators = [Z3, Z3 @ Z3, np.eye(3)]
    long_idle_channel_weights = [0, 0.00007751, 0.99992249] 
    short_idle_channel_weights = [0, 0.0000000019379, 0.99999999806]

    single_qutrit_kraus_operators = [X3, X3 @ X3, Z3, Z3 @ Z3, Y3, Y3 @ Y3, V3, V3 @ V3, np.eye(3)]
    single_qutrit_kraus_operator_weights = [8.076953189667142e-7,8.356530070058442e-7,4.4910858870426395e-7,4.4910858870426395e-7,8.076953189667142e-7,8.356530070058442e-7,8.076953189667142e-7,8.356530070058442e-7,0.999846624680123]

    XX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    X2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    ZX_weights = [2.512069977868293e-10,2.5990231483578203e-10,1.3968041859275328e-10,1.3968041859275328e-10,2.512069977868293e-10,2.5990231483578203e-10,2.512069977868293e-10,2.5990231483578203e-10,0.000011816834382389788]
    Z2X_weights = [2.512069977868293e-10,2.5990231483578203e-10,1.3968041859275328e-10,1.3968041859275328e-10,2.512069977868293e-10,2.5990231483578203e-10,2.512069977868293e-10,2.5990231483578203e-10,0.000011816834382389788]
    YX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    Y2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    VX_weights = [4.517809752636722e-10,4.674189903317726e-10,2.512069977868293e-10,2.512069977868293e-10,4.517809752636722e-10,4.674189903317726e-10,4.517809752636722e-10,4.674189903317726e-10,0.000021251879958915915]
    V2X_weights = [4.674189903317726e-10,4.835983020207133e-10,2.5990231483578203e-10,2.5990231483578203e-10,4.674189903317726e-10,4.835983020207133e-10,4.674189903317726e-10,4.835983020207133e-10,0.00002198749574891023]
    TX_weights = [0.000021251879958915915,0.00002198749574891023,0.000011816834382389788,0.000011816834382389788,0.000021251879958915915,0.00002198749574891023,0.000021251879958915915,0.00002198749574891023,0.9996932728842347]

    XX_operators = [np.kron(X3, X3), np.kron(X3, X3 @ X3), np.kron(X3, Z3), np.kron(X3, Z3 @ Z3), np.kron(X3, Y3), np.kron(X3, Y3 @ Y3), np.kron(X3, V3), np.kron(X3, V3 @ V3), np.kron(X3, np.eye(3))]
    X2X_operators = [np.kron(X3 @ X3, X3), np.kron(X3 @ X3, X3 @ X3), np.kron(X3 @ X3, Z3), np.kron(X3 @ X3, Z3 @ Z3), np.kron(X3 @ X3, Y3), np.kron(X3 @ X3, Y3 @ Y3), np.kron(X3 @ X3, V3), np.kron(X3 @ X3, V3 @ V3), np.kron(X3 @ X3, np.eye(3))]
    ZX_operators = [np.kron(Z3, X3), np.kron(Z3, X3 @ X3), np.kron(Z3, Z3), np.kron(Z3, Z3 @ Z3), np.kron(Z3, Y3), np.kron(Z3, Y3 @ Y3), np.kron(Z3, V3), np.kron(Z3, V3 @ V3), np.kron(Z3, np.eye(3))]
    Z2X_operators = [np.kron(Z3 @ Z3, X3), np.kron(Z3 @ Z3, X3 @ X3), np.kron(Z3 @ Z3, Z3), np.kron(Z3 @ Z3, Z3 @ Z3), np.kron(Z3 @ Z3, Y3), np.kron(Z3 @ Z3, Y3 @ Y3), np.kron(Z3 @ Z3, V3), np.kron(Z3 @ Z3, V3 @ V3), np.kron(Z3 @ Z3, np.eye(3))]
    YX_operators = [np.kron(Y3, X3), np.kron(Y3, X3 @ X3), np.kron(Y3, Z3), np.kron(Y3, Z3 @ Z3), np.kron(Y3, Y3), np.kron(Y3, Y3 @ Y3), np.kron(Y3, V3), np.kron(Y3, V3 @ V3), np.kron(Y3, np.eye(3))]
    Y2X_operators = [np.kron(Y3 @ Y3, X3), np.kron(Y3 @ Y3, X3 @ X3), np.kron(Y3 @ Y3, Z3), np.kron(Y3 @ Y3, Z3 @ Z3), np.kron(Y3 @ Y3, Y3), np.kron(Y3 @ Y3, Y3 @ Y3), np.kron(Y3 @ Y3, V3), np.kron(Y3 @ Y3, V3 @ V3), np.kron(Y3 @ Y3, np.eye(3))]
    VX_operators = [np.kron(V3, X3), np.kron(V3, X3 @ X3), np.kron(V3, Z3), np.kron(V3, Z3 @ Z3), np.kron(V3, Y3), np.kron(V3, Y3 @ Y3), np.kron(V3, V3), np.kron(V3, V3 @ V3), np.kron(V3, np.eye(3))]
    V2X_operators = [np.kron(V3 @ V3, X3), np.kron(V3 @ V3, X3 @ X3), np.kron(V3 @ V3, Z3), np.kron(V3 @ V3, Z3 @ Z3), np.kron(V3 @ V3, Y3), np.kron(V3 @ V3, Y3 @ Y3), np.kron(V3 @ V3, V3), np.kron(V3 @ V3, V3 @ V3), np.kron(V3 @ V3, np.eye(3))]
    TX_operators = [np.kron(np.eye(3), X3), np.kron(np.eye(3), X3 @ X3), np.kron(np.eye(3), Z3), np.kron(np.eye(3), Z3 @ Z3), np.kron(np.eye(3), Y3), np.kron(np.eye(3), Y3 @ Y3), np.kron(np.eye(3), V3), np.kron(np.eye(3), V3 @ V3), np.kron(np.eye(3), np.eye(3))]

    two_qutrit_kraus_operators = XX_operators + X2X_operators + ZX_operators + Z2X_operators + YX_operators + Y2X_operators + VX_operators + V2X_operators + TX_operators
    two_qutrit_kraus_operator_weights =  XX_weights + X2X_weights + ZX_weights + Z2X_weights + YX_weights + Y2X_weights + VX_weights + V2X_weights + TX_weights


class QutritPlusGate(cirq.Gate):
    """A gate that adds one in the computational basis of a qutrit.

    This gate acts on three-level systems. In the computational basis of
    this system it enacts the transformation U|x〉 = |x + 1 mod 3〉, or
    in other words U|0〉 = |1〉, U|1〉 = |2〉, and U|2> = |0〉.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix.
        return np.array([[0, 0, 1],
                         [1, 0, 0],
                         [0, 1, 0]])

    def _circuit_diagram_info_(self, args):
        return '[+1]'


class QutritControlPlusGate(cirq.Gate):
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3, 3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix.
        return np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1]])

    def _circuit_diagram_info_(self, args):
        return '[+1]'

def get_random_state(n, dim=3):
    rand_state2 = np.zeros(2 ** n, np.complex64)
    rand_state2.real = np.random.randn(2 ** n)
    rand_state2.imag = np.random.randn(2 ** n)
    rand_state2 /= np.linalg.norm(rand_state2)

    rand_state3 = np.zeros(dim ** n, np.complex64)
    for i, val in enumerate(rand_state2):
        rand_state3[int(bin(i)[2:], dim)] = val
    new_state = np.zeros((dim,)*n, dtype=np.complex128)
    if dim == 3:
        for c, s in zip(product("012", repeat=n), rand_state3):
            current_ndarray = new_state
            for i in c[:-1]:
                index = int(i)
                current_ndarray = current_ndarray[index]
            current_ndarray[int(c[-1])] = s
    else:
        for c, s in zip(product("01", repeat=n), rand_state3):
            current_ndarray = new_state
            for i in c[:-1]:
                index = int(i)
                current_ndarray = current_ndarray[index]
            current_ndarray[int(c[-1])] = s
    return new_state

def qutrits_decomposed(num=4):
    qutrit = cirq.LineQid.range(num, dimension=3)
    circ = cirq.Circuit()
    if num == 4:
        circ.append(apply_C1C1PlusOne(qutrit, 0, 2, 1))
        circ.append(QutritC2PlusGate().on(qutrit[1], qutrit[3]))
        circ.append(apply_C1C1MinusOne(qutrit, 0, 2, 1))
    if num == 8:
        circ.append(apply_C1C1PlusOne(qutrit, 0, 2, 1)) 
        circ.append(apply_C1C1PlusOne(qutrit, 4, 6, 5))

        circ.append(apply_C2C2PlusOne(qutrit, 1, 5, 3))
        circ.append(QutritC2PlusGate().on(qutrit[3], qutrit[7]))
        circ.append(apply_C2C2MinusOne(qutrit, 1, 5, 3))
        
        circ.append(apply_C1C1MinusOne(qutrit, 0, 2, 1)) 
        circ.append(apply_C1C1MinusOne(qutrit, 4, 6, 5))
    if num == 16:
        circ.append(apply_C1C1PlusOne(qutrit, 0, 2, 1))
        circ.append(apply_C1C1PlusOne(qutrit, 4, 6, 5))
        circ.append(apply_C1C1PlusOne(qutrit, 8, 10, 9))
        circ.append(apply_C1C1PlusOne(qutrit, 12, 14, 13))

        circ.append(apply_C2C2PlusOne(qutrit, 1, 5, 3))
        circ.append(apply_C2C2PlusOne(qutrit, 9, 13, 11))
        circ.append(apply_C2C2PlusOne(qutrit, 3, 11, 7))
        circ.append(QutritC2PlusGate().on(qutrit[7], qutrit[15]))
        circ.append(apply_C2C2MinusOne(qutrit, 1, 5, 3))
        circ.append(apply_C2C2MinusOne(qutrit, 9, 13, 11))
        circ.append(apply_C2C2MinusOne(qutrit, 3, 11, 7))

        circ.append(apply_C1C1MinusOne(qutrit, 0, 2, 1))
        circ.append(apply_C1C1MinusOne(qutrit, 4, 6, 5))
        circ.append(apply_C1C1MinusOne(qutrit, 8, 10, 9))
        circ.append(apply_C1C1MinusOne(qutrit, 12, 14, 13))


    if False:
        circ.append(generate_incrementer(qutrit, 0, 2, 1, 1))
        circ.append(generate_incrementer(qutrit, 4, 6, 5, 1))
        circ.append(generate_incrementer(qutrit, 8, 10, 9, 1))
        circ.append(generate_incrementer(qutrit, 12, 14, 13, 1))


        circ.append(generate_incrementer(qutrit, 1, 5, 3, 2))
        circ.append(generate_incrementer(qutrit, 9, 13, 11, 2))
        circ.append(generate_incrementer(qutrit, 3, 11, 7, 2))
        circ.append(QutritC2PlusGate().on(q[7], q[15]))
        circ.append(generate_decrementer(qutrit, 3, 11, 7, 2))
        circ.append(generate_decrementer(qutrit, 9, 13, 11, 2))
        circ.append(generate_decrementer(qutrit, 1, 5, 3, 2))
        
        circ.append(generate_decrementer(qutrit, 12, 14, 13, 1))
        circ.append(generate_decrementer(qutrit, 8, 10, 9, 1))
        circ.append(generate_decrementer(qutrit, 4, 6, 5, 1))
        circ.append(generate_decrementer(qutrit, 0, 2, 1, 1))
    return circ

# cirq.unitary(op)

if __name__ == "__main__":
    example_noise_model = FutureSuperconductingQCErrors()  # CurrentSuperconductingQCErrors()
    # FutureSuperconductingQCErrorsBetterT1 FutureSuperconductingQCErrorsBetterGates  FutureSuperconductingQCErrorsBetterT1AndGates
    nums = [4] # , 8]
    while len(nums) != 0:
        num  = nums.pop(0)
        q = cirq.LineQid.range(num, dimension=3)
        if num == 4:
             circuit = cirq.Circuit(

                        QutritCC1PlusGate().on(q[0], q[2], q[1]), 
                        QutritC2PlusGate().on(q[1], q[3]),
                        QutritCC1MinusGate().on(q[0], q[2], q[1]), 
            )
        elif num == 8:
              circuit = cirq.Circuit(

                        QutritCC1PlusGate().on(q[0], q[2], q[1]),
                        QutritCC1PlusGate().on(q[4], q[6], q[5]),
                        
                        QutritCC2PlusGate().on(q[1], q[5], q[3]),
                        QutritC2PlusGate().on(q[3], q[7]),
                        QutritCC2MinusGate().on(q[1], q[5], q[3]), 
                        
                        QutritCC1MinusGate().on(q[0], q[2], q[1]),
                        QutritCC1MinusGate().on(q[4], q[6], q[5]),
            )
        elif num == 16:
            circuit = cirq.Circuit(

                        QutritCC1PlusGate().on(q[0], q[2], q[1]),
                        QutritCC1PlusGate().on(q[4], q[6], q[5]),
                        QutritCC1PlusGate().on(q[8], q[10], q[9]),
                        QutritCC1PlusGate().on(q[12], q[14], q[13]),
                        QutritCC2PlusGate().on(q[1], q[5], q[3]),
                        QutritCC2PlusGate().on(q[9], q[13], q[11]),
                        QutritCC2PlusGate().on(q[3], q[11], q[7]),
                        QutritC2PlusGate().on(q[7], q[15]),

                        QutritCC2MinusGate().on(q[3], q[11], q[7]), 
                        
                        QutritCC2MinusGate().on(q[1], q[5], q[3]),
                        QutritCC2MinusGate().on(q[9], q[13], q[11]),

                        QutritCC1MinusGate().on(q[0], q[2], q[1]),
                        QutritCC1MinusGate().on(q[4], q[6], q[5]),
                        QutritCC1MinusGate().on(q[8], q[10], q[9]),
                        QutritCC1MinusGate().on(q[12], q[14], q[13]),
            )
        circ2 =  cirq.Circuit()
        circ2.append(qutrits_decomposed(num))
        results = []

        
        for i in range(100):
            if i % 100 == 0:
                print(i)
            random_state = get_random_state(num)
            state = np.zeros(random_state.shape, dtype ='complex')
            state[1, 1, 1, 0] = 1
            # random_state = state
            ideal_state = apply_unitaries_to_circuit(circ2, random_state)
            noisy_state = apply_unitaries_to_circuit(circ2, random_state, noise_model=example_noise_model)
            result = np.linalg.norm(np.vdot(ideal_state, noisy_state)) ** 2
            results.append(result)
        print(sum(results) / len(results))

