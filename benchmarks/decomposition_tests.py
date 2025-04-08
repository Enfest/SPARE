import numpy as np
from scipy.stats import unitary_group
from quantumcircuitbenchmarks.quantumcircuitbenchmarks.cirq import *
import cirq
from src.toffoli import *
from src.controlled_gate_impl import *
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from rewriter.utils.check_subspace_equality import return_qubit_circuit
from src.create_qutrit_circ import *
import random
from ripl.decomp.util.circuit_utils import save_circuit
from replacement_lib.replace_qubit_circ import replace_qubit_circuit
from benchmarks.get_square_benchmarks import get_simple_circuits, get_errorcorr_circuits, get_twooffive_circuits, get_sixsim_circuits, get_rd53_circuits, get_adder_circuits
from src.dotConvertor import dotConvertor
# from test_qubit_decomposition import run_qubit_sim

def check_if_block_diagonal(matrix, D=3):
    assert matrix.shape[0] % 3 == 0
    dim_shape = matrix.shape[0] // D
    if dim_shape == 1:
        return np.allclose(matrix, np.diag(np.diag(matrix)))
    result = []
    for i in range(D):
        for j in range(D):
            if i != j and not np.allclose(matrix[i * dim_shape : (i+1) * dim_shape,
                                                 j * dim_shape : (j+1) * dim_shape], np.zeros((dim_shape, dim_shape))):
                return False, [None, None, None]
            elif i == j:
                result.append(matrix[i * dim_shape: (i+1) * dim_shape, i * dim_shape: (i+1) * dim_shape].copy())
    # print("blk diag")
    return True, result

def compose_block_diagonals(matrix, bds, D=3):
    assert check_if_block_diagonal(matrix)
    assert D == len(bds)
    result = np.zeros(matrix.shape, dtype=np.complex128)
    dim_shape = matrix.shape[0] // D
    for i in range(D):
        assert bds[i].shape[0] == dim_shape
        result[i * dim_shape: (i+1)*dim_shape, i * dim_shape: (i+1)*dim_shape] = bds[i].copy()
    return result.copy()

def attain_symmetry(matrix):
    if matrix.shape[0] == 3:
        return matrix.copy()
    result, [bd1, bd2, bd3] = check_if_block_diagonal(matrix)
    if result:
        if np.allclose(bd1, bd2) and np.allclose(bd2, bd3):
            return matrix
        if np.allclose(bd1, bd2):
            bd3 = bd1.copy()
            result_bd1 = attain_symmetry(bd1)
            result_bd2 = attain_symmetry(bd2)
            result_bd3 = attain_symmetry(bd3)
            matrix_result = compose_block_diagonals(matrix, [result_bd1, result_bd2, result_bd3])
        else:
            bd3 = bd1.copy()
            result_bd1 = attain_symmetry(bd1)
            result_bd2 = attain_symmetry(bd2)
            result_bd3 = attain_symmetry(bd3)
            matrix_result = compose_block_diagonals(matrix, [result_bd1, result_bd2, result_bd3])    
        return matrix_result.copy()
    else:
        return matrix.copy()
    

def get_circuit_unitaries(circuit, dim=3, dtype=np.complex128):
    old_dim = dim
    dim = 3
    qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
    min_val = 100000000000
    for a in qd_map.keys():
        qd_map[a] = a.x
        min_val = min(min_val, a.x)
    for a in qd_map.keys():
        qd_map[a] = a.x - min_val

    num_qudits = len(qd_map.keys())
    n = num_qudits
    new_state = np.zeros((dim,)*n, dtype=np.complex128)
    unitary_matrix = np.zeros((dim ** n, dim ** n), dtype=np.complex128)
    for i in range(dim ** n) :
        temp_state = np.zeros(dim ** n, dtype=np.complex128)
        temp_state[i] = 1.0
        state = temp_state.reshape(new_state.shape)
        buffer = np.zeros(state.shape, dtype=dtype)
        for m_index, moment in enumerate(circuit):
            for op in moment:
                qudit_indices = [qd_map[qd] for qd in op.qubits]
                mat = cirq.unitary(op)
                if old_dim == 2 and mat.shape[0] == 2:
                    matrix2 = np.eye(3, dtype=np.complex128)
                    matrix2[0:2, 0:2] = mat
                    mat = matrix2
                elif old_dim == 2 and mat.shape[0] == 4:
                    matrix2 = np.eye(9, dtype=np.complex128)
                    matrix2[0:2, 0:2] = mat[0:2, 0:2]
                    matrix2[3:5, 3:5] = mat[2:4, 2:4]
                    
                    mat = matrix2
                elif old_dim == 2 and mat.shape[0] == 8:
                    matrix2 = np.eye(3 ** len(op.qubits), dtype=np.complex128)
                    # for i in range()
                    matrix2[0:2, 0:2] = mat[0:2, 0:2]
                    matrix2[3:5, 3:5] = mat[2:4, 2:4]
                    matrix2[9:11, 9:11] = mat[4:6, 4:6]
                    matrix2[12:14, 12:14] = mat[6:8, 6:8]
                    mat = matrix2
                elif old_dim == 3:
                    pass
                else:
                    assert False
                # matrix =  mat.astype(dtype).reshape((dim,)*((3)*len(op.qubits)))
                matrix =  mat.astype(dtype).reshape((dim,)*(2*len(op.qubits)))
                cirq.linalg.targeted_left_multiply(matrix, state, target_axes=qudit_indices, out=buffer)
                state, buffer = buffer, state
        state = state.flatten()
        unitary_matrix[:, i] = state[:]
    
    assert np.allclose(np.eye(len(unitary_matrix)), unitary_matrix.dot(unitary_matrix.T.conj()))
    return unitary_matrix

def get_circuit_unitaries_2(circuit, dim=2, dtype=np.complex128):
    old_dim = dim
    dim = 2
    qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
    min_val = 100000000000
    for a in qd_map.keys():
        qd_map[a] = a.x
        min_val = min(min_val, a.x)
    for a in qd_map.keys():
        qd_map[a] = a.x - min_val
    num_qudits = len(qd_map.keys())
    n = num_qudits
    new_state = np.zeros((dim,)*n, dtype=np.complex128)
    unitary_matrix = np.zeros((dim ** n, dim ** n), dtype=np.complex128)
    for i in range(dim ** n) :
        temp_state = np.zeros(dim ** n, dtype=np.complex128)
        temp_state[i] = 1.0
        state = temp_state.reshape(new_state.shape)
        buffer = np.zeros(state.shape, dtype=dtype)
        a = 0
        for m_index, moment in enumerate(circuit):
            for op in moment:
                qudit_indices = [qd_map[qd] for qd in op.qubits]
                mat = cirq.unitary(op)

                matrix =  mat.astype(dtype).reshape((dim,)*((2)*len(op.qubits)))
                cirq.linalg.targeted_left_multiply(matrix, state, target_axes=qudit_indices, out=buffer)
                state, buffer = buffer, state
        state = state.flatten()
        unitary_matrix[:, i] = state[:]
    assert np.allclose(np.eye(len(unitary_matrix)), unitary_matrix.dot(unitary_matrix.T.conj()))
    return unitary_matrix

def return_toffoli_3(size=3, dimension=3, return_base=False, return_qubit=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    circ.append(QutritX12Gate().on(q[1]).controlled_by(q[0]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritCX101Gate().on(q[1], q[2]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[1]).controlled_by(q[0]))
    print(circ)
    matrix = get_circuit_unitaries(circ)
    circ2 = cirq.Circuit()
    circ2.append(cirq.MatrixGate(matrix, qid_shape=[dimension]*size).on(* q))
    if return_base:
        return circ, q
    if return_qubit:
        qubit = cirq.LineQid.range(size-1, dimension=2)
        mat = np.zeros((2, 2), dtype=np.complex128)
        mat[0, 1] = 1
        mat[1, 0] = 1
        circ.append(cirq.MatrixGate(mat).on(qubit[-1]).controlled_by(* qubit[:-1]))

        return
    return circ2, q

def return_toffoli_4(size=4, dimension=3, return_base=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    
    circ.append(QutritC1PlusGate().on(q[0], q[1]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX01Gate().on(q[3]).controlled_by(q[0], q[1], q[2]))
    
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritC1MinusGate().on(q[0], q[1]))
    
    circ.append(cirq.MatrixGate(np.eye(3), qid_shape=[dimension]).on(q[-1]))
    print(circ)

    matrix = get_circuit_unitaries(circ)
    circ2 =  cirq.Circuit()

    circ2.append(cirq.MatrixGate(matrix, qid_shape=[dimension]*(size)).on(* q))
    if return_base:
        return circ, q
    return circ2, q

def return_toffoli_5(size=5, dimension=3, return_base=False):
    # change this to 5
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    circ.append(QutritC1PlusGate().on(q[0], q[1]))
    circ.append(QutritC1PlusGate().on(q[2], q[3]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[3]))
    circ.append(QutritX01Gate().on(q[4]).controlled_by(q[1], q[3]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[3]))
    circ.append(QutritC1MinusGate().on(q[0], q[1]))
    circ.append(QutritC1MinusGate().on(q[2], q[3]))
    matrix = get_circuit_unitaries(circ)
    circ2 =  cirq.Circuit()
    circ2.append(cirq.MatrixGate(matrix, qid_shape=[dimension]*size).on(* q))
    if return_base:
        return circ, q
    return circ2, q

def return_toffoli_6(size=6, dimension=3, return_base=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    circ.append(QutritC1PlusGate().on(q[0], q[1]))
    
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritC1PlusGate().on(q[1], q[2]))
    circ.append(QutritX12Gate().on(q[1]))
    
    circ.append(QutritC1MinusGate().on(q[0], q[1]))
    # circ.append(QutritCC1PlusGate().on(q[0], q[1], q[2]))
    circ.append(QutritC1PlusGate().on(q[3], q[4]))

    circ.append(QutritX12Gate().on(q[2]))
    circ.append(QutritX12Gate().on(q[4]))
    circ.append(QutritX01Gate().on(q[5]).controlled_by(q[2], q[4]))
    circ.append(QutritX12Gate().on(q[2]))
    circ.append(QutritX12Gate().on(q[4]))

    circ.append(QutritC1MinusGate().on(q[3], q[4]))
    circ.append(QutritC1PlusGate().on(q[0], q[1]))
    
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritC1MinusGate().on(q[1], q[2]))
    circ.append(QutritX12Gate().on(q[1]))
    
    circ.append(QutritC1MinusGate().on(q[0], q[1]))
    
    print(circ)
    # circ.append(cirq.MatrixGate(np.eye(3), qid_shape=[3]).on(q[-1]))
    # circ.append(QutritCC1MinusGate().on(q[0], q[1], q[2]))
    matrix = get_circuit_unitaries(circ)
    circ2 =  cirq.Circuit()
    circ2.append(cirq.MatrixGate(matrix, qid_shape=[dimension]*(size)).on(* q))

    if return_base:
        return circ, q
    return circ2, q


def return_incrementer_4(size=4, dimension=3, return_base=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    circ2 = cirq.Circuit()
    circ.append(QutritPlusGate().on(q[0]))
    circ.append(QutritX12Gate().on(q[0]))
    circ.append(QutritCC1PlusGate().on(q[0], q[1], q[2]))
    circ.append(QutritX01Gate().on(q[1]).controlled_by(q[0]))
    circ.append(QutritX12Gate().on(q[2]))
    circ.append(QutritX01Gate().on(q[3]).controlled_by(q[2]))
    circ.append(QutritX12Gate().on(q[2]))
    circ.append(QutritX01Gate().on(q[1]))
    circ.append(QutritX02Gate().on(q[2]).controlled_by(* [q[0], q[1]]))
    circ.append(QutritX01Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[0]))
    circ.append(QutritX02Gate().on(q[0]))
    # get_and_print_dag(circ, q, "orig_incrementer_4")

    circ.append(cirq.MatrixGate(np.eye(3), qid_shape=[dimension]).on(q[-1]))

    print(circ)
    matrix = get_circuit_unitaries(circ)
    circ2.append(cirq.MatrixGate(matrix, qid_shape=[dimension]*size).on(* q))
    if return_base:
        return circ, q
    return circ2, q

def return_incrementer_4_(size=4, dimension=3, return_base=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    circ2 = cirq.Circuit()
    circ.append(QutritPlusGate().on(q[0]))
    circ.append(QutritX12Gate().on(q[0]))
    circ.append(QutritCC1PlusGate().on(q[0], q[1], q[2]))
    circ.append(QutritX01Gate().on(q[1]).controlled_by(q[0]))
    circ.append(QutritX12Gate().on(q[2]))
    circ.append(QutritX01Gate().on(q[3]).controlled_by(q[2]))
    circ.append(QutritX12Gate().on(q[2]))
    circ.append(QutritX01Gate().on(q[1]))
    circ.append(QutritX02Gate().on(q[2]).controlled_by(* [q[0], q[1]]))
    circ.append(QutritX01Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[0]))
    circ.append(QutritX02Gate().on(q[0]))
    print(circ)
    matrix = get_circuit_unitaries(circ)
    circ2.append(cirq.MatrixGate(matrix, qid_shape=[dimension]*size).on(* q))
    if return_base:
        return circ, q
    return circ2, q

def return_incrementer_5(size=5, dimension=3, return_base=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    circ2 = cirq.Circuit()
    # circ.append(cirq.MatrixGate(np.eye(3), qid_shape=[dimension]).on(q[5]))
    
    circ.append(QutritPlusGate().on(q[0]))
    
    circ.append(QutritX12Gate().on(q[0]))
    
    circ.append(QutritX12Gate().on(q[1]).controlled_by(q[0]))
    circ.append(QutritX12Gate().on(q[3]).controlled_by(q[2]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[3]))
    circ.append(QutritPlusGate().on(q[4]).controlled_by(q[3], q[1])) #.controlled_by(* [q[0], q[1], q[2]]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[3]))
    circ.append(QutritX12Gate().on(q[1]).controlled_by(q[0]))
    circ.append(QutritX12Gate().on(q[3]).controlled_by(q[2]))
    
    circ.append(QutritX12Gate().on(q[1]).controlled_by(q[0]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritPlusGate().on(q[2]).controlled_by(q[1]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[1]).controlled_by(q[0]))

    circ.append(QutritX12Gate().on(q[2]))
    circ.append(QutritX01Gate().on(q[3]).controlled_by(q[2]))
    circ.append(QutritX12Gate().on(q[2]))

    circ.append(QutritX01Gate().on(q[1]).controlled_by(q[0]))
    
    circ.append(QutritX01Gate().on(q[1]))
    circ.append(QutritCX112Gate().on(q[0], q[1]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX02Gate().on(q[2]).controlled_by(q[1]))
    circ.append(QutritX12Gate().on(q[1]))    
    circ.append(QutritCX112Gate().on(q[0], q[1]))
    circ.append(QutritX01Gate().on(q[1]))
    
    circ.append(QutritX01Gate().on(q[1]))
    circ.append(QutritX01Gate().on(q[2]))
    circ.append(QutritX01Gate().on(q[3]))
    circ.append(QutritX12Gate().on(q[1]).controlled_by(q[0]))
    circ.append(QutritX12Gate().on(q[3]).controlled_by(q[2]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[3]))
    circ.append(QutritX02Gate().on(q[4]).controlled_by(q[1], q[3])) #.controlled_by(* [q[0], q[1], q[2]]))
    circ.append(QutritX12Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[3]))
    circ.append(QutritX12Gate().on(q[1]).controlled_by(q[0]))
    circ.append(QutritX12Gate().on(q[3]).controlled_by(q[2]))
    circ.append(QutritX01Gate().on(q[1]))
    circ.append(QutritX01Gate().on(q[2]))
    circ.append(QutritX01Gate().on(q[3]))
    circ.append(QutritX12Gate().on(q[0]))
    circ.append(QutritX02Gate().on(q[0]))
    matrix = get_circuit_unitaries(circ)
    circ2.append(cirq.MatrixGate(matrix, qid_shape=[dimension]*size).on(* q))
    if return_base:
        return circ, q
    return circ2, q

def return_incrementer_3(size=3, dimension=3):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    circ2 = cirq.Circuit()
    # circ.append(cirq.MatrixGate(np.eye(3), qid_shape=[dimension]).on(q[3]))
    circ.append(QutritPlusGate().on(q[0]))
    circ.append(QutritX12Gate().on(q[0]))
    circ.append(QutritCC1PlusGate().on(q[0], q[1], q[2]))
    circ.append(QutritX01Gate().on(q[1]).controlled_by(q[0]))
    circ.append(QutritX01Gate().on(q[1]))
    circ.append(QutritX02Gate().on(q[2]).controlled_by(* [q[0], q[1]]))
    circ.append(QutritX01Gate().on(q[1]))
    circ.append(QutritX12Gate().on(q[0]))
    circ.append(QutritX02Gate().on(q[0]))
    print(circ)
    matrix = get_circuit_unitaries(circ)
    print(matrix)
    for i in range(matrix.shape[0]):
        print(i, matrix[i, :])
    circ2.append(cirq.MatrixGate(matrix, qid_shape=[dimension]*size).on(* q))
    return circ2, q

def return_takahashi_4(size=4, dimension=3, return_base=False, return_qubit=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    dummy_circ = generate_takahashi_adder(size, to_toffoli=True)
    mat = get_circuit_unitaries(dummy_circ, dim=2)
    mat2 = get_circuit_unitaries_2(dummy_circ, dim=2)
    
    if return_base:
        res_circuit, res_qutrits = create_qutrit_from_qubit(dummy_circ, len(q))
        return res_circuit, res_qutrits
    
    if return_qubit:
        # get_and_print_dag(res_circuit, res_qutrits, "orig_takahashi_4")
        q = cirq.LineQid.range(size, dimension=dimension)
        dummy_circ2 = cirq.Circuit()
        dummy_circ2.append(replace_qubit_circuit(dummy_circ, q))
        single_gate_depth, double_gate_depth = depth_counter(dummy_circ2, q)
        mat = get_circuit_unitaries_2(dummy_circ, dim=2)
        return dummy_circ2, q
    
    mat_printing = return_qubit_circuit(mat)
    assert np.allclose(mat_printing, mat2)
    
    for i in range(mat2.shape[0]):
        print(mat2[i, :])
    circ_dummy = cirq.Circuit()
    circ_dummy.append(QutritCX101Gate().on(q[0], q[2]))
    circ_dummy.append(QutritX02Gate().on(q[2]).controlled_by(q[0]))
    circ_dummy.append(QutritCX201Gate().on(q[2], q[3]))
    circ_dummy.append(QutritX02Gate().on(q[2]).controlled_by(q[0]))
    circ_dummy.append(QutritX01Gate().on(q[3]).controlled_by(q[1]))

    circ_dummy.append(cirq.MatrixGate(np.eye(3, dtype=np.complex128), qid_shape=[3]*1).on(q[1]))
    mat3 = get_circuit_unitaries(circ_dummy)
    print(mat)
    new_matrix = mat.copy() # attain_symmetry(mat.copy())
    mat_printing3 = return_qubit_circuit(new_matrix)
    assert np.allclose(mat_printing3, mat_printing)
    print(dummy_circ)
    
    circ.append(cirq.MatrixGate(new_matrix, qid_shape=[dimension]*(size)).on(* q))
    return circ, q

def return_takahashi_8(size=8, dimension=3, return_base=False, return_qubit=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size , dimension=dimension)
    dummy_circ = generate_takahashi_adder(size, to_toffoli=True)

    if return_base:
        res_circuit, res_qutrits = create_qutrit_from_qubit(dummy_circ, len(q))
        return res_circuit, res_qutrits
        # get_and_print_dag(res_circuit, res_qutrits, "orig_takahashi_8")
    
    if return_qubit:
        dummy_circ2 = cirq.Circuit()
        dummy_circ2.append(replace_qubit_circuit(dummy_circ, q))
        single_gate_depth, double_gate_depth = depth_counter(dummy_circ2, q)
        return dummy_circ2, cirq.LineQid.range(size, dimension=2)
 
    mat = get_circuit_unitaries(dummy_circ, dim=2)
    mat2 = get_circuit_unitaries_2(dummy_circ, dim=2)
    
    mat_printing = return_qubit_circuit(mat)
    print("early_check")
    print(mat_printing)
    print("MATRIX IS THE FOLLOWING FOR TAKAHASHI 4")
    assert np.allclose(mat_printing, mat2)
    
    circ_dummy = cirq.Circuit()
    circ_dummy.append(QutritCX101Gate().on(q[0], q[2]))
    circ_dummy.append(QutritX02Gate().on(q[2]).controlled_by(q[0]))
    circ_dummy.append(QutritCX201Gate().on(q[2], q[3]))
    circ_dummy.append(QutritX02Gate().on(q[2]).controlled_by(q[0]))
    circ_dummy.append(QutritX01Gate().on(q[3]).controlled_by(q[1]))

    circ_dummy.append(cirq.MatrixGate(np.eye(3, dtype=np.complex128), qid_shape=[3]*1).on(q[1]))
    new_matrix = attain_symmetry(mat.copy())
    mat_printing3 = return_qubit_circuit(new_matrix)
    assert np.allclose(mat_printing3, mat_printing)
    print(dummy_circ)
    circ.append(cirq.MatrixGate(new_matrix, qid_shape=[dimension]*size).on(* q))
    return circ, q

def return_takahashi_6(size=6, dimension=3, return_base=False, return_qubit=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    dummy_circ = generate_takahashi_adder(size, to_toffoli=True)

    if return_base:
        res_circuit, res_qutrits = create_qutrit_from_qubit(dummy_circ, len(q))
        return res_circuit, res_qutrits
    
    if return_qubit:
        # get_and_print_dag(res_circuit, res_qutrits, "orig_takahashi_6")
        dummy_circ2 = cirq.Circuit()
        dummy_circ2.append(replace_qubit_circuit(dummy_circ, q))
        single_gate_depth, double_gate_depth = depth_counter(dummy_circ2, q)
        return dummy_circ2, cirq.LineQid.range(size, dimension=2)
        
    mat = get_circuit_unitaries(dummy_circ, dim=2)    
    print("MATRIX IS THE FOLLOWING FOR TAKAHASHI 6")
    print(dummy_circ)
    new_matrix = attain_symmetry(mat.copy())
    
    m1 = return_qubit_circuit(new_matrix)
    m2 = return_qubit_circuit(mat)
    assert np.allclose(m1, m2)
    circ.append(cirq.MatrixGate(new_matrix, qid_shape=[dimension]*size).on(* q))
    return circ, q

def return_cuccaro_6(size=6, dimension=3, return_base=False, return_qubit=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    dummy_circ = generate_cuccaro_adder(size, to_toffoli=True)
    # Curcarro 6 
    if return_base:
        res_circuit, res_qutrits = create_qutrit_from_qubit(dummy_circ, len(q))
        return res_circuit, res_qutrits
    
    if return_qubit:
        dummy_circ2 = cirq.Circuit()
        dummy_circ2.append(replace_qubit_circuit(dummy_circ, q))
        single_gate_depth, double_gate_depth = depth_counter(dummy_circ2, q)
        return dummy_circ2, cirq.LineQid.range(size, dimension=2)
    
    mat2 = get_circuit_unitaries_2(dummy_circ, dim=2)
    mat = get_circuit_unitaries(dummy_circ, dim=2)
    new_matrix = attain_symmetry(mat.copy())
    m1 = return_qubit_circuit(new_matrix)
    m2 = return_qubit_circuit(mat)
    assert np.allclose(m1, m2)
    circ.append(cirq.MatrixGate(new_matrix, qid_shape=[dimension]*size).on(* q))
    return circ, q

def return_cuccaro_4(size=4, dimension=3, return_base=False, return_qubit=False):
    circ = cirq.Circuit()
    q = cirq.LineQid.range(size, dimension=dimension)
    dummy_circ = generate_cuccaro_adder(size, to_toffoli=True)
    if return_base:
        res_circuit, res_qutrits = create_qutrit_from_qubit(dummy_circ, len(q))
        return res_circuit, res_qutrits
    
    if return_qubit:
        dummy_circ2 = cirq.Circuit()
        dummy_circ2.append(replace_qubit_circuit(dummy_circ, q)) 
        single_gate_depth, double_gate_depth = depth_counter(dummy_circ2, q)
        return dummy_circ2, cirq.LineQid.range(size, dimension=2)
            
    mat2 = get_circuit_unitaries_2(dummy_circ, dim=2)
    print("MATRIX IS THE FOLLOWING FOR Cucarro 4")
    for i in range(mat2.shape[0]):
        print(mat2[i, :])
    # plt.imshow(np.asarray(mat2, dtype=np.float), cmap='viridis')
    # plt.savefig("ccc4.png")
    mat = get_circuit_unitaries(dummy_circ, dim=2)
    print(dummy_circ)
    new_matrix = attain_symmetry(mat.copy()) # mat.copy() # 
    circ.append(cirq.MatrixGate(new_matrix, qid_shape=[dimension]*size).on(* q))
    return circ, q


class QubitXGate(cirq.Gate):
    """A gate that adds one in the computational basis of a qutrit.
    This gate acts on two-level systems. In the computational basis of
    this system it enacts the transformation U|0〉= |0 + 1 / √2〉 and 
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        return (2,)

    def _num_qubits_(self):
        return 1

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        return np.array([[0, 1],
                         [1, 0]], dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return "X"


def return_random_circuit(size=6, dimension=3, depth=6, return_base=False):
    for j in range(100):
        circuit = cirq.Circuit()
        qubits = cirq.LineQid.range(size, dimension=2)
        Xgate = np.eye(2, dtype=np.complex128)
        Xgate[0, 0] = 0
        Xgate[1, 0] = 1
        Xgate[1, 1] = 0
        Xgate[0, 1] = 1
        seed = 0
        random.seed(j)
        arr_collection = []

        for i in range(depth):
            arr = []
            while len(arr) < 3:
                num1 = random.random()
                num1 = int(num1 * size)
                if num1 not in arr:
                    arr.append(num1)
            arr_collection.append(arr)

        blocked_arr = list(range(size))            
        for k in range(2):
            if k == 0:
                for elem in arr_collection[:-1]:
                    arr = elem
                    arr.sort()
                    qubits_control = [qubits[idx] for idx in arr[0:2]]
                    print(qubits[arr[-1]])
                    print(qubits_control)
                    
                    circuit.append(QubitXGate().on(qubits[arr[-1]]).controlled_by(* qubits_control))
                    for x in arr:
                        if x in blocked_arr:
                            blocked_arr.remove(x)
                if len(blocked_arr) > 2:
                    assert False
                if size - 1 in blocked_arr:
                    assert False
                if len(blocked_arr) == 2:
                    circuit.append(QubitXGate().on(qubits[size - 1]).controlled_by(* [qubits[blocked_arr[0]], qubits[blocked_arr[1]]]))
                elif len(blocked_arr) == 1 and blocked_arr[0] < size - 2:
                    circuit.append(QubitXGate().on(qubits[size - 1]).controlled_by(* [qubits[blocked_arr[0]], qubits[blocked_arr[0] + 1]]))
                elif len(blocked_arr) == 1 and blocked_arr[0] == size - 2:
                    circuit.append(QubitXGate().on(qubits[size - 1]).controlled_by(* [qubits[blocked_arr[0] - 1], qubits[blocked_arr[0]]]))    
                else:
                    circuit.append(QubitXGate().on(qubits[size - 1]).controlled_by(* [qubits[0], qubits[1]]))
            if k == 2:
                copy = arr_collection[:-1].copy()
                # copy.reverse()
                for elem in copy:
                    arr = elem
                    arr.sort()
                    qubits_control = [qubits[idx] for idx in arr[0:2]]
                    circuit.append(cirq.MatrixGate(Xgate).on(qubits[arr[-1]]).controlled_by(* qubits_control))    
        
        circuit.append(cirq.MatrixGate(np.eye(2)).on(qubits[-1]))
        mat2 = get_circuit_unitaries_2(circuit, dim=2)

        if return_base:
            res_circuit, res_qutrits = create_qutrit_from_qubit(circuit, len(qubits))
            # get_and_print_dag(res_circuit, res_qutrits, "qutrit_replace_random_" + str(depth + 1))
            return res_circuit, res_qutrits

        if return_base:
            dummy_circ2 = cirq.Circuit()
            dummy_circ2.append(replace_qubit_circuit(circuit, qubits))
            # run_qubit_sim(dummy_circ2, size)
            single_gate_depth, double_gate_depth = depth_counter(dummy_circ2, qubits)
            # get_and_print_dag(dummy_circ2, qubits, "random_circ_qubit_" + str(num))
            assert False

        matx = get_circuit_unitaries(circuit, dim=2)
        res, res_arr = check_if_block_diagonal(matx, D=3)
        if res:
            new_matrix = attain_symmetry(matx)
            q = cirq.LineQid.range(size, dimension=dimension)
            final_circ = cirq.Circuit()
            final_circ.append(cirq.MatrixGate(new_matrix, qid_shape=[dimension]*size).on(* q))
            return final_circ, q
        assert False

def return_simple4(size=5, dimension=3, return_base=False):
    matrix = get_simple_circuits()
    c = cirq.Circuit()
    all_qubits = cirq.LineQid.range(size, dimension=3)
    print("all qubit s ", all_qubits)
    matrix = attain_symmetry(matrix)
    c.append(cirq.MatrixGate(matrix, qid_shape=[3]*size).on(* all_qubits))
    return c, all_qubits

def return_twooffive(size=6, dimension=3, return_base=True):
    circuit, qubits_info, matrix = get_twooffive_circuits(return_base)
    print("circuit check")
    print(circuit)
    print("BASE INFO ", return_base)
    if return_base:
        '''Save the circuit into a dot file, which qubitsn are ancilla'''
        dotConvertor(circuit, filename="dot_files/twooffive/twooffive_graph",
                ancillas=list(range(len(circuit.all_qubits())))[-qubits_info:])
        return circuit, cirq.LineQid.range(len(circuit.all_qubits()), dimension=2)
    circuit.all_qubits()
    c = cirq.Circuit()
    all_qubits = cirq.LineQid.range(size, dimension=3)
    matrix = attain_symmetry(matrix)
    c.append(cirq.MatrixGate(matrix, qid_shape=[3]*size).on(* all_qubits))
    return c, all_qubits

def return_errorcorr(size=3, dimension=3, return_base=False):
    circuit, qubits_info, matrix = get_errorcorr_circuits()
    if return_base:
        '''Save the circuit into a dot file, which qubitsn are ancilla'''
        dotConvertor(circuit, filename="dot_files/errorcorr/errorcorr_graph",
                        ancillas=list(range(len(circuit.all_qubits())))[-qubits_info:])
        return circuit, cirq.LineQid.range(len(circuit.all_qubits()), dimension=2)
    c = cirq.Circuit()
    all_qubits = cirq.LineQid.range(size, dimension=3)
    print("all qubit s ", all_qubits)
    matrix = attain_symmetry(matrix)
    c.append(cirq.MatrixGate(matrix, qid_shape=[3]*size).on(* all_qubits))
    return c, all_qubits

def return_sixsim(size=7, dimension=3, return_base=False):
    circuit, qubits_info, matrix = get_sixsim_circuits(return_base)
    if return_base:
        '''Save the circuit into a dot file, which qubitsn are ancilla'''
        dotConvertor(circuit, filename="dot_files/sixsim/sixsim_graph",
                ancillas=list(range(len(circuit.all_qubits())))[-qubits_info:])
        return circuit, cirq.LineQid.range(len(circuit.all_qubits()), dimension=2)

    c = cirq.Circuit()
    all_qubits = cirq.LineQid.range(size, dimension=3)
    matrix = attain_symmetry(matrix)
    c.append(cirq.MatrixGate(matrix, qid_shape=[3]*size).on(* all_qubits))
    return c, all_qubits

def return_rd53(size=8, dimension=3, return_base=False):
    circuit, qubits_info, matrix = get_rd53_circuits(return_base)
    if return_base:
        dotConvertor(circuit, filename="dot_files/rd53/rd53_graph",
                ancillas=list(range(len(circuit.all_qubits())))[-qubits_info:])
        return circuit, cirq.LineQid.range(len(circuit.all_qubits()), dimension=2)

    c = cirq.Circuit()
    all_qubits = cirq.LineQid.range(size, dimension=3)
    matrix = attain_symmetry(matrix)
    c.append(cirq.MatrixGate(matrix, qid_shape=[3]*size).on(* all_qubits))
    return c, all_qubits

def return_adder(size=6, dimension=3, return_base=False):
    print("ADDER CHWECK")
    matrix = get_adder_circuits()
    c = cirq.Circuit()
    all_qubits = cirq.LineQid.range(size, dimension=3)
    matrix = attain_symmetry(matrix)
    c.append(cirq.MatrixGate(matrix, qid_shape=[3]*size).on(* all_qubits))
    return c, all_qubits

def return_decomp(size=8, dimension=3, return_base=False):
    c = cirq.Circuit()
    qubits = cirq.LineQid.range(size, dimension=3)
    # c.append(QutritPlusGate().on(qubits[0]).controlled_by(qubits[1]))
    # c.append(QutritX12Gate().on(qubits[0]))
    # c.append(QutritX01Gate().on(qubits[1]).controlled_by(qubits[0]))
    # c.append(QutritX12Gate().on(qubits[0]))
    c.append(QutritX01Gate().on(qubits[3]).controlled_by(qubits[2]))
    c.append(QutritX01Gate().on(qubits[5]).controlled_by(qubits[4]))
    c.append(QutritX01Gate().on(qubits[7]).controlled_by(qubits[6]))
    c.append(QutritX01Gate().on(qubits[2]).controlled_by(* [qubits[0], qubits[1]]))
    c.append(QutritX01Gate().on(qubits[4]).controlled_by(qubits[3], qubits[2]))
    c.append(QutritX01Gate().on(qubits[6]).controlled_by(qubits[5], qubits[4]))
    c.append(QutritX01Gate().on(qubits[7]).controlled_by(qubits[6]))
    c.append(QutritX01Gate().on(qubits[6]).controlled_by(qubits[5], qubits[4]))
    
    c.append(QutritX01Gate().on(qubits[5]).controlled_by(qubits[4]))
    c.append(QutritX01Gate().on(qubits[4]).controlled_by(qubits[3], qubits[2]))
    c.append(QutritX01Gate().on(qubits[3]).controlled_by(qubits[2]))
    c.append(QutritX01Gate().on(qubits[2]).controlled_by(qubits[0], qubits[1]))
    c.append(QutritX01Gate().on(qubits[1]).controlled_by(qubits[0]))
    c.append(QutritX01Gate().on(qubits[4]).controlled_by(qubits[2]))

    c.append(QutritX01Gate().on(qubits[3]).controlled_by(qubits[2]))
    c.append(QutritX01Gate().on(qubits[6]).controlled_by(qubits[4]))

    c.append(QutritX01Gate().on(qubits[5]).controlled_by(qubits[4]))
    c.append(QutritX01Gate().on(qubits[7]).controlled_by(qubits[6]))
    print(c)
    matrix = cirq.unitary(c)
    print(matrix)
    matrix_new = attain_symmetry(matrix.copy())
    c2 = cirq.Circuit()
    c2.append(cirq.MatrixGate(matrix_new, qid_shape=[3]*size).on(* qubits))
    return c2, qubits

# def __name__ == "__main__":

