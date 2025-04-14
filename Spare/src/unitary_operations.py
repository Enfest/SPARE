import numpy as np
import cirq
from itertools import product
import re

def get_circuit_unitaries(circuit, dimension=3, debug=False, dtype=np.complex128):
    dim = dimension
    qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
    min_val = 100000000000
    max_val = -1
    if len(circuit.all_qubits()) > 8 and dim == 3:
        # This is a uncheckable case simply return
        return np.eye(3, dtype=np.complex128)
    if len(circuit.all_qubits()) > 12 and dim == 2:
        return np.eye(2, dtype=np.complex128) 

    for a in qd_map.keys():
        if isinstance(a, cirq.LineQid):
            qd_map[a] = a.x
        elif isinstance(a, cirq.NamedQubit):
            qd_map[a] = int(re.findall(r'\d+', a.name)[-1])
        min_val = min(min_val, qd_map[a])
        max_val = max(max_val, qd_map[a])

    keys = list(qd_map.keys())
    keys.sort()
    for i, k in enumerate(keys):
        qd_map[k] = i
    
    num_qudits = len(qd_map) # max_val - min_val + 1 # len(qd_map.keys())
    if num_qudits < 0:
        # Raise Exception
        return np.eye(dim)
    num_qudits = len(qd_map.keys())
    n = num_qudits
    n_ = 0
    new_state = np.zeros((dim,)*n, dtype=np.complex128)
    unitary_matrix = np.zeros((dim ** n, dim ** n), dtype=np.complex128)
    j = 0
    for i in range(dim ** n) :
        # print(i, unitary_matrix.shape)
        temp_state = np.zeros(dim ** n, dtype=np.complex128)
        temp_state[i] = 1.0
        state = temp_state.reshape(new_state.shape)
        #if i > 15:
        #    print(temp_state)
        buffer = np.zeros(state.shape, dtype=dtype)
        for m_index, moment in enumerate(circuit):
            if i == 1:
                j += 1
            for op in moment:
                qudit_indices = [qd_map[qd] for qd in op.qubits]
                mat = cirq.unitary(op)
                matrix =  mat.astype(dtype).reshape((dim,)*(2*len(op.qubits)))
                cirq.linalg.targeted_left_multiply(matrix, state, target_axes=qudit_indices, out=buffer)
                state, buffer = buffer, state
        state = state.flatten()
        unitary_matrix[:, i] = state[:]
    # for i in range(unitary_matrix.shape[0] // 27):
    #     if not np.allclose(unitary_matrix.dot(unitary_matrix.T.conj())[i*27:27*(i+1), i*27:27*(i+1)], np.eye(27)):
    #         print("range is ", i*27, (i+1)*27)
    #         print(unitary_matrix.dot(unitary_matrix.T.conj())[i*27:27*(i+1), i*27:27*(i+1)])
    assert np.linalg.norm(unitary_matrix.dot(unitary_matrix.T.conj()) - np.eye(unitary_matrix.shape[0])) <= 1e-3
    # assert np.allclose(np.eye(len(unitary_matrix)), unitary_matrix.dot(unitary_matrix.T.conj()), rtol=1e-3)
    return unitary_matrix

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

# def get_random_state(n):
#     rand_state2 = np.zeros(2 ** n, np.complex64)
#     rand_state2.real = np.random.randn(2 ** n)
#     rand_state2.imag = np.random.randn(2 ** n)
#     rand_state2 /= np.linalg.norm(rand_state2)
#     rand_state3 = np.zeros(3 ** n, np.complex64)
#     for i, val in enumerate(rand_state2):
#         rand_state3[int(bin(i)[2:], 3)] = val
#     new_state = np.zeros((3,)*n, dtype=np.complex128)
#     for c, s in zip(product("012", repeat=n), rand_state3):
#         current_ndarray = new_state
#         for i in c[:-1]:
#             index = int(i)
#             current_ndarray = current_ndarray[index]
#         current_ndarray[int(c[-1])] = s
#     return new_state

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
