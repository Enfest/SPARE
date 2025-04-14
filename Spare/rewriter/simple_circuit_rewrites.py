import cirq
import numpy as np
from collections import defaultdict
from src.toffoli import *
from src.unitary_operations import get_circuit_unitaries
import re

def merge_same_qubit_gates(circuit, qutrits, d=3, end=None, iffinal=False, dtype=np.complex128, checking_blocked=False):
    # print("initial circuit ")
    # print(circuit)
    if_merge = False
    directory_if_index_processed = {}
    new_circuit = cirq.Circuit()
    other_circuit = cirq.Circuit()
    if iffinal:
        pass
    for m_index, moment in enumerate(circuit):
        directory_if_index_processed[m_index] = False
        if end is not None and m_index == end:
            break
    for m_index, moment in enumerate(circuit):
        if end is not None and m_index >= end:
            break
        other_circuit.append(moment)
        if not directory_if_index_processed[m_index]:
            qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
            for q in qd_map:
                if isinstance(q, cirq.LineQid):
                    qd_map[q] = q.x
                elif isinstance(q, cirq.NamedQubit):
                    qd_map[q] = int(re.findall(r'\d+', q.name)[-1])

            qubits_used_by_same_gate = [None] * len(qutrits)
            gate_matrix_arr = {}
            operation_states = {}
            # For curr moment
            moment_arr = {}
            for gates in moment:
                for q in gates.qubits:
                    qubits_used_by_same_gate[qd_map[q]] = gates.qubits
                    if q == gates.qubits[-1]:
                        moment_arr[qd_map[q]] = cirq.unitary(gates)
                        gate_matrix_arr[qd_map[q]] = cirq.unitary(gates)
                        operation_states[qd_map[q]] = gates.qubits 
            
            process_next = True
            k = 0
            moment_arr2 = moment_arr.copy()
            while process_next and m_index + k + 1 < len(directory_if_index_processed.keys()):
                k += 1
                moment2 = circuit[m_index + k]
                for gates in moment2:
                    # print(qd_map[gates.qubits[-1]], qubits_used_by_same_gate)\
                    # assert qd_map[gates.qubits[-1]] in qubits_used_by_same_gate 
                    
                    if qd_map[gates.qubits[-1]] in qubits_used_by_same_gate and qubits_used_by_same_gate[qd_map[gates.qubits[-1]]] is not None and \
                            len(gates.qubits) == len(qubits_used_by_same_gate[qd_map[gates.qubits[-1]]]):
                        # Check if gates in kast moment is same as curr moment
                        for q in gates.qubits:
                            if qubits_used_by_same_gate[qd_map[q]] is None or q not in qubits_used_by_same_gate[qd_map[gates.qubits[-1]]]:
                                process_next = False
                        # Only if we can proceed
                        if process_next:
                            if gates.qubits[-1] != qubits_used_by_same_gate[qd_map[gates.qubits[-1]]][-1]:
                                assert False
                                process_next = False
                        if process_next:
                            q = gates.qubits[-1]
                            g1 = moment_arr2[qd_map[q]]
                            g2 = cirq.unitary(gates)
                            moment_arr2[qd_map[q]] = g2 @ g1
                            if_merge = True
                    else:
                        process_next = False
                if process_next:
                    directory_if_index_processed[m_index + k] = True
                    # New moment post merging is
                    moment_arr = moment_arr2
            for elems in moment_arr.keys():
                if moment_arr[elems].shape[0] != d:
                    new_gate = moment_arr[elems][(moment_arr[elems].shape[0] - d) // (d-1): (moment_arr[elems].shape[0] + d) // (d-1),
                                                 (moment_arr[elems].shape[0] - d) // (d-1): (moment_arr[elems].shape[0] + d) // (d-1)]
                    if len(operation_states[elems]) > 1 and not np.allclose(np.eye(moment_arr[elems].shape[0]), moment_arr[elems]):
                        
                        if np.allclose(new_gate, cirq.unitary(QutritX01Gate())[:new_gate.shape[0], :new_gate.shape[0]]): 
                            if d ==2:
                                new_circuit.append(cirq.X(operation_states[elems][-1]).controlled_by(* operation_states[elems][:-1]))            
                            elif d == 3:
                                new_circuit.append(QutritX01Gate().on(operation_states[elems][-1]).controlled_by(* operation_states[elems][:-1]))            
                                
                        else:
                            new_circuit.append(cirq.MatrixGate(new_gate,
                                                               qid_shape=[d]).on(operation_states[elems][-1]).controlled_by(* operation_states[elems][:-1]))
                else:
                    if not np.allclose(np.eye(moment_arr[elems].shape[0]), moment_arr[elems]):
                        new_circuit.append(cirq.MatrixGate(moment_arr[elems], qid_shape=[d]).on(operation_states[elems][-1]))
    # print(new_circuit)
    if not checking_blocked:
        unitary1 = get_circuit_unitaries(circuit, dim=d)
        unitary2 = get_circuit_unitaries(new_circuit, dim=d)
        assert np.allclose(unitary1, unitary2)
    return new_circuit, if_merge

def merge_controlled_and_single_gates(circuit, qutrits, d=3, dtype=np.complex128, checking_blocked=False):
    merged = False
    qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
    min_val = float('inf')
    for a in qd_map.keys():
        if isinstance(a, cirq.LineQid):
            qd_map[a] = a.x
            min_val = min(min_val, a.x)
        elif isinstance(a, cirq.NamedQubit):
            qd_map[a] = int(re.findall(r'\d+', a.name)[-1])
            min_val = min(min_val, qd_map[a])
    
    new_circuit = cirq.Circuit()
    gate_to_qubits = defaultdict(dict)
    gate_to_matrix = defaultdict(dict)
    default_matrix = np.eye(d, dtype=np.complex128)
    for m_index, moment in enumerate(circuit):
        for op in moment:
            gate_to_matrix[m_index][op.qubits[-1]] = op
            gate_to_qubits[m_index][op.qubits[-1]] = op.qubits

    for m_index, moment in enumerate(circuit):
        new_moment = {}
        for op in moment:
            new_moment[op.qubits[-1]] = op
            if len(op.qubits) > 1 and m_index > 0:
                # Already seen a gate and is a controlled gate
                if op.qubits[-1] in gate_to_matrix[m_index - 1].keys() and\
                        gate_to_matrix[m_index - 1][op.qubits[-1]] is not None and\
                        len(gate_to_qubits[m_index - 1][op.qubits[-1]]) == 1:
                    mat1 = cirq.unitary(gate_to_matrix[m_index - 1][op.qubits[-1]])
                else:
                    mat1 = default_matrix.copy()
                if op.qubits[-1] in gate_to_matrix[m_index + 1].keys() and\
                        gate_to_matrix[m_index + 1][op.qubits[-1]] is not None and\
                        len(gate_to_qubits[m_index + 1][op.qubits[-1]]) == 1:
                    mat2 = cirq.unitary(gate_to_matrix[m_index + 1][op.qubits[-1]])
                else:
                    mat2 = default_matrix.copy()
                # If both are not identity
                if not np.allclose(mat1, np.eye(d)) and np.allclose(mat1 @ mat2, np.eye(d)):
                    mat = cirq.unitary(new_moment[op.qubits[-1]])
                    gate_to_matrix[m_index][op.qubits[-1]] = cirq.MatrixGate(mat2 @ mat[((mat.shape[0] - d) // (d - 1)) : ((mat.shape[0] + d) // (d - 1)),
                                                                                        ((mat.shape[0] - d) // (d - 1)) : ((mat.shape[0] + d) // (d - 1))] @
                                                                            mat1, qid_shape=[d]).on(op.qubits[-1]).controlled_by(* op.qubits[:-1])
                    gate_to_matrix[m_index - 1][op.qubits[-1]] = None
                    gate_to_matrix[m_index + 1][op.qubits[-1]] = None
                    merged = True
    
    for m_index in gate_to_matrix.keys():
        for q in gate_to_matrix[m_index]:
            if gate_to_matrix[m_index][q] is not None:
                new_circuit.append(gate_to_matrix[m_index][q])
    if not checking_blocked:
        unitary1 = get_circuit_unitaries(circuit, dim=d)
        unitary2 = get_circuit_unitaries(new_circuit, dim=d)
        assert np.allclose(unitary1, unitary2)
    final_circuit = cirq.Circuit()
    final_circuit.append(merge_same_qubit_gates(new_circuit,qutrits,checking_blocked=checking_blocked, d=d)[0])
    return final_circuit, merged

def change_2q_gates_to_controlled(circuit, qutrits, checking_blocked, d=3):
    new_circuit = cirq.Circuit()
    for m_index, moments in enumerate(circuit):
        for op in moments:
            if len(op.qubits) == 2:
                mat = cirq.unitary(op)
                if np.allclose(np.eye(d), mat[d:d*2, d:d*2]):
                    if not np.allclose(np.eye(d), mat[-d:, -d:]):
                        new_circuit.append(QutritX12Gate().on(op.qubits[0]))
                        if np.allclose(mat[-3:, -3:], cirq.unitary(QutritX01Gate())[-3:, -3:]):
                            new_circuit.append(QutritX01Gate().on(op.qubits[-1]).controlled_by(* op.qubits[:-1]))
                        else:
                            new_circuit.append(cirq.MatrixGate(mat[-3:, -3:], qid_shape=[d]).on(op.qubits[-1]).controlled_by(* op.qubits[:-1]))
                        new_circuit.append(QutritX12Gate().on(op.qubits[0]))
                    if not np.allclose(np.eye(d), mat[:d, :d]):
                        new_circuit.append(QutritX01Gate().on(op.qubits[0]))
                        new_circuit.append(cirq.MatrixGate(mat[:d, :d], qid_shape=[d]).on(op.qubits[-1]).controlled_by(* op.qubits[:-1]))
                        new_circuit.append(QutritX01Gate().on(op.qubits[0]))
                else:
                    new_circuit.append(cirq.MatrixGate(mat[d:2*d, d:2*d], qid_shape=[d]).on(op.qubits[-1]).controlled_by(* op.qubits[:-1]))
            else:
                new_circuit.append(op)
    if not checking_blocked:
        unitary1 = get_circuit_unitaries(circuit, dim=d)
        unitary2 = get_circuit_unitaries(new_circuit, dim=d)
        assert np.allclose(unitary1, unitary2)
    return new_circuit

def loop_through_rewrites(circuit, qutrits, dtype=np.complex128, checking_blocked=False, d=3):
    if qutrits is None:
        qutrits = circuit.all_qubits()
    circuit = change_2q_gates_to_controlled(circuit, qutrits, checking_blocked, d=d)
    old_circ = circuit
    for i in range(4):
        circ = cirq.Circuit()
        c, res = merge_controlled_and_single_gates(old_circ, qutrits, checking_blocked=checking_blocked, d=d)
        circ.append(c)
        old_circ = circ
        if not res:
            break
    old_finl_circ = circ # finl_circ
    for i in range(4):
        circ = cirq.Circuit()
        c, res = merge_controlled_and_single_gates(old_finl_circ, qutrits, checking_blocked=checking_blocked, d=d)
        circ.append(c)
        # circ.append(merge_controlled_and_single_gates(old_finl_circ, qutrits))
        old_finl_circ = circ
        if not res:
            break
    for i in range(4):
        circ2 = cirq.Circuit()
        c, res = merge_same_qubit_gates(circ, qutrits, checking_blocked=checking_blocked, d=d)
        circ2.append(c)
        circ = circ2
        if not res:
            break
    return circ2
