import cirq 
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit import Aer, assemble, execute
from qiskit.quantum_info.synthesis import two_qubit_decompose
# from qiskit.circuit.library import RCCXGate
from rewriter.utils.node_gate_types import GateOpType
from rewriter.utils.graph_utils import if_mat_preserves_binary
import cmath
import src.qutrit_decompositions as gateimpl

def matrix_convert(gate, gatetypeinfo, base=False, debug=True):
    new_matrix = np.eye(4, dtype=np.complex128)
    if base:
        new_matrix[2, 2] = 0
        new_matrix[2, 3] = 1
        new_matrix[3, 2] = 1
        new_matrix[3, 3] = 0
        if debug:
            print(new_matrix)
        return new_matrix
    if if_mat_preserves_binary(gate @ cirq.unitary(gateimpl.QutritX12Gate())):
        new_matrix[:2, :2] = gate[:2, :2]
        new_matrix[3, :2] = gate[2, :2].reshape((2,))
        new_matrix[3, 3] = gate[2, 2]
        new_matrix[:2, 3] = gate[:2, 2].reshape((2,))
        # print("for x12 gate we get ", new_matrix)
    elif if_mat_preserves_binary(gate @ cirq.unitary(gateimpl.QutritMinusGate())):
        for i in range(new_matrix.shape[0]):
            new_matrix[i, i] = 0
            new_matrix[2, 1] = 1
        new_matrix[:2, 2:] = gate[:2, 1:]
        new_matrix[:2, 0] =  gate[:2, 0]
        new_matrix[3, 0] = gate[2, 0]
        new_matrix[3, 2:] = gate[2, 1:].reshape((2,))
    elif if_mat_preserves_binary(gate @ cirq.unitary(gateimpl.QutritPlusGate())):
        for i in range(new_matrix.shape[0]):
            new_matrix[i, i] = 0
        new_matrix[2, 2] = 1
        new_matrix[:2, :2] = gate[:2, :2]
        new_matrix[3, :2] =  gate[2, :2]
        new_matrix[:2, 3] = gate[:2, 2]
        new_matrix[3, 3] = gate[2, 2]
    elif if_mat_preserves_binary(gate @ cirq.unitary(gateimpl.QutritX02Gate())):
        for i in range(new_matrix.shape[0]):
            new_matrix[i, i] = 0
        new_matrix[1:, :3] = gate[:, :]
        new_matrix[3, 0] =  1
    else:
        print(gate)
        raise NotImplementedError()     
    if debug:
        print(new_matrix)
    return new_matrix

def circuit_gen(matrix, original):
    qc = QuantumCircuit(2)
    is_rccx = False
    if if_mat_preserves_binary(original @ cirq.unitary(gateimpl.QutritX12Gate())):
        # qc.append(RCCXGate(), [0, 1])  
        qc.cx(0, 1)
        return qc, True
    elif if_mat_preserves_binary(original @ cirq.unitary(gateimpl.QutritMinusGate())):
        qc.cx(0, 1)
        qc.x(0)
        return qc, is_rccx
    elif if_mat_preserves_binary(original @ cirq.unitary(gateimpl.QutritPlusGate())):
        qc.x(0)
        qc.cx(0, 1)
        return qc, is_rccx
    elif if_mat_preserves_binary(original @ cirq.unitary(gateimpl.QutritX02Gate())):
        qc.x(0)
        qc.x(1)
        return qc, is_rccx
    elif if_mat_preserves_binary(original):
        return qc, is_rccx
    qc.unitary(matrix, [0, 1])
    trans_qc = transpile(qc, optimization_level=3, basis_gates=['cx', 'x', 'y', 'z', 'h', 's', 't', 'u']) # , decomposition=[two_qubit_decompose])
    return trans_qc

# Function to get the unitary matrix of a gate
def get_gate_unitary(gate, qubits):
    gatecopy = gate.copy()
    qc = QuantumCircuit(qubits)
    gatecopy.qubits = qc.qubits
    qc.append(gatecopy, qubits)
    simulator = Aer.get_backend('unitary_simulator')
    compiled_circuit = transpile(qc, simulator)
    result = execute(compiled_circuit, simulator).result()
    return np.array(result.get_unitary().data)

def unroll_controls(gate, new_qubits):
    for i, c in enumerate(gate.get_control()):
        if gate.get_control_value_for_qubit(c) == 2:
            gate.replace_qubit(c, new_qubits[c])
            gate.set_control_value_for_qubit(new_qubits[c], 1)
    # print("check the status of the conversion ", gate)
    return gate

def gate_convertor(orig_gate, new_qubits, gatetypeinfo, ancilla_data):
    if orig_gate.get_gate_type().is_binary():
        new_matrix = orig_gate.get_matrix()[:2, :2]
        # print("Check status ", orig_gate.get_target(), ancilla_data[orig_gate.get_target()])
        if ancilla_data[orig_gate.get_target()] and len(orig_gate.get_control()) > 1:
            if np.allclose(orig_gate.get_matrix()[:2, :2], np.asarray([[0, 1], [1, 0]])):
                new_matrix = np.asarray([[0, -1], [1, 0]])
        orig_gate_copy = orig_gate.deepcopy()
        if orig_gate_copy.get_matrix().shape[0] == 3:
            assert np.allclose(orig_gate_copy.get_matrix()[2, 2] * orig_gate_copy.get_matrix()[2, 2].conjugate(), 1)
            orig_gate_copy.set_matrix(new_matrix, assert_disabled=True)  # rig_gate_copy.get_matrix()[:2, :2]
        else:
            orig_gate_copy.set_matrix(new_matrix, assert_disabled=True)  # orig_gate_copy.get_matrix()
        result = unroll_controls(orig_gate_copy, new_qubits)
        return result
    # print("target qubits could be ", new_qubits, orig_gate.get_target(), new_qubits[orig_gate.get_target()])
    new_target = new_qubits[orig_gate.get_target()]
    targets = [orig_gate.get_target(), new_target]
    circuit, is_rccx = circuit_gen(matrix_convert(gatetypeinfo.normalize_matrix(orig_gate.get_matrix()), gatetypeinfo, debug=False), original=orig_gate.get_matrix())
    
    # print("Check status ", new_target, orig_gate.get_target(), ancilla_data[orig_gate.get_target()], new_target < len(ancilla_data) and ancilla_data[new_target])
    
    is_rccx = is_rccx or (new_target < len(ancilla_data) and ancilla_data[new_target])
    gates = []
    new_gate = orig_gate.deepcopy()
    phase_angle = circuit.global_phase
    new_gate = unroll_controls(new_gate, new_qubits)
    new_gate.set_matrix(np.asarray([[cmath.exp(complex(0, phase_angle)), 0], [0, cmath.exp(complex(0, phase_angle))]]), assert_disabled=True)
    new_gate.set_target(new_target)
    new_gate.delete_last_qubit(new_target)
    if not np.allclose(new_gate.get_matrix(), np.eye(2)):
        gates.append(new_gate)
    for i, gate in enumerate(circuit.data):
        new_gate = orig_gate.deepcopy()
        if len(gate.qubits) == 1:
            if np.allclose(get_gate_unitary(gate, 1), np.eye(2)):
                continue
            gate.qubits = (gate.qubits[0], 0)
            new_gate.set_matrix(get_gate_unitary(gate, 1), assert_disabled=True)
            # Ancilla target
            new_gate.set_target(targets[gate.qubits[0].index])
        elif gate.operation.name == "cx" and (not is_rccx):
            new_gate.set_matrix(np.asarray([[0, 1],[1, 0]], dtype=np.complex128), assert_disabled=True)
            new_gate.update_qubit_targets(control=targets[gate.qubits[0].index], target=targets[gate.qubits[1].index])
            # print("1 ", new_gate)
        elif gate.operation.name == "cx" and (is_rccx):          
            new_gate.set_matrix(np.asarray([[0, -1],[1, 0]], dtype=np.complex128), assert_disabled=True)
            new_gate.update_qubit_targets(control=targets[gate.qubits[0].index], target=targets[gate.qubits[1].index])     
            # print("2 ", new_gate)
        else:
            assert False
        gates.append(new_gate)
    result = [unroll_controls(g, new_qubits) for g in gates]
    return result
