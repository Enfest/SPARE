import cirq
import numpy as np
from cmath import exp
from math import pi
import math
# import cirq.setup
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate as UGate
from qiskit.circuit.library import CCXGate, RYGate, RXGate, CXGate, RCCXGate, CCZGate
from qiskit.circuit import Gate
# from qiskit.transpiler import PassManager
from qiskit.compiler import transpile
from src.cirq_depth_counter import depth_counter

class QubitHGate(cirq.Gate):
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
        ones = np.ones((2, 2), dtype=np.complex128)
        ones[1, 1] = -1
        ones = ones / np.sqrt(np.linalg.norm(ones))
        # print("Matrix is")
        # print(ones)
        return ones

    def _circuit_diagram_info_(self, args):
        return "H"

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

class QubitTGate(cirq.Gate):
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
        return np.array([[1, 0],
                         [0, exp(complex(0, 1)  * pi / 4)]], dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return "T"

class QubitT_Gate(cirq.Gate):
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

        return np.array([[1, 0],
                         [0, exp(-1 * complex(0, 1) * pi / 4)]], dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return "T+"
    
class QubitSGate(cirq.Gate):
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
        return np.array([[1, 0],
                         [0, 1*j]], dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return "S"
    
class QubitCNOTGate(cirq.Gate):
    """A gate that adds one in the computational basis of a qutrit.
    This gate acts on two-level systems. In the computational basis of
    this system it enacts the transformation U|0〉= |0 + 1 / √2〉 and 
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        return (2, 2)

    def _num_qubits_(self):
        return 2

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        return np.array([[1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 0, 1],
                         [0, 0, 1, 0]], dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return '(1)', 'X'

class QubitCNOTGate(cirq.Gate):
    """A gate that adds one in the computational basis of a qutrit.
    This gate acts on two-level systems. In the computational basis of
    this system it enacts the transformation U|0〉= |0 + 1 / √2〉 and 
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        return (2, 2)

    def _num_qubits_(self):
        return 2

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        return np.array([[1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, -1]], dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return '(1)', 'Z'


def toffoli_circ(qubits):
    qubits = list(qubits)
    circuit = cirq.Circuit()
    Hgate =  np.ones((2, 2), dtype=np.complex128)
    Hgate[1, 1] = -1
    circuit.append(QubitHGate().on(qubits[2]))
    circuit.append(QubitCNOTGate().on(qubits[1], qubits[2]))
    circuit.append(QubitT_Gate().on(qubits[2]))
    circuit.append(QubitCNOTGate().on(qubits[0], qubits[2]))
    circuit.append(QubitTGate().on(qubits[2]))
    circuit.append(QubitCNOTGate().on(qubits[1], qubits[2]))
    circuit.append(QubitT_Gate().on(qubits[2]))
    circuit.append(QubitCNOTGate().on(qubits[0], qubits[2]))
    circuit.append(QubitTGate().on(qubits[2]))
    circuit.append(QubitHGate().on(qubits[2]))
    circuit.append(QubitTGate().on(qubits[1]))
    circuit.append(QubitCNOTGate().on(qubits[0], qubits[1]))
    circuit.append(QubitTGate().on(qubits[0]))
    circuit.append(QubitT_Gate().on(qubits[1]))
    circuit.append(QubitCNOTGate().on(qubits[0], qubits[1]))
    return circuit 

# def replace_qubit_circuit(circuit, qubits=None):
#     final_circ = cirq.Circuit()
#     xgate = np.asarray([[0, 1],
#                         [1, 0]], dtype=np.complex128)
#     # clifford_t_circuit = cirq.optimizers.decompose_clifford_t(circuit)
#     # return clifford_t_circuit
#     for m_index, moment in enumerate(circuit):
#         for op in moment:
#             print(op)
#             if cirq.unitary(op).shape[0] == 2:
#                 # final_circ.append(op)
#                 U = cirq.unitary(op)
#                 q = op.qubits[0]
#                 if not np.allclose(U, xgate):
#                     print("MISMATCH")
#                 if np.allclose(U, xgate):  # Identity
#                     final_circ.append(cirq.X(q))
#                 else:
#                     assert False
#             elif cirq.unitary(op).shape[0] == 4:
#                 control = op.qubits[0]
#                 target = op.qubits[1]
#                 matrix =  cirq.unitary(op)[-2:, -2:]
#                 final_circ.append(cirq.MatrixGate(matrix).on(target).controlled_by(control))
#             elif cirq.unitary(op).shape[0] == 8:
#                 matrix = cirq.unitary(op)[-2:, -2:]                
#                 if np.allclose(matrix, xgate):
#                     final_circ.append(toffoli_circ(op.qubits))
#                 else:
#                     raise Exception("Not supported for decomposing yet 3 qubits")
#             else:
#                 print("last operation was ")
#                 raise Exception("Not supported for decomposing yet")
#     return final_circ

# Function to convert Cirq to Qiskit
def cirq_to_qiskit(cirq_circuit):
    """Convert Cirq circuit to Qiskit circuit."""
    qiskit_circuit = QuantumCircuit(1 + max([x.x for x in cirq_circuit.all_qubits()]))
    for moment in cirq_circuit:
        for gate in moment:
            qubits = [qiskit_circuit.qubits[qubit.x] for qubit in gate.qubits]
            qiskit_gate = _convert_cirq_gate_to_qiskit(gate)
            if isinstance(qiskit_gate, list):
                for g, qs in qiskit_gate:
                    # print([x.index for x in qs])
                    g_qs = [qubits[x.index] for x in qs]
                    # print(g_qs)
                    qiskit_circuit.append(g, g_qs)
            else:
                qiskit_circuit.append(qiskit_gate, qubits)
    # print(qiskit_circuit)
    return qiskit_circuit

def is_ry_gate(matrix, theta_tolerance=1e-12):
    if matrix.shape != (2, 2) or np.imag(matrix[0, 0]) != 0 or np.imag(matrix[0, 1]) != 0:
        return (False, None)
    theta = 2*np.arccos(matrix[0, 0])
    cos_theta = np.cos(theta/2)
    sin_theta = np.sin(theta/2)
    # print(matrix[0, 0], theta, cos_theta)
    assert np.allclose(matrix[0, 0], cos_theta)
    expected_matrix = np.array([[cos_theta, -sin_theta], [sin_theta, cos_theta]])
    return (np.allclose(matrix, expected_matrix, atol=theta_tolerance), np.real(theta).astype(np.float64))

def is_rx_gate(matrix, theta_tolerance=1e-12):
    if matrix.shape != (2, 2) or np.imag(matrix[0, 0]) != 0 or np.real(matrix[0, 1]) != 0:
        return (False, None)
    theta = 2*np.arccos(matrix[0, 0])
    cos_theta = np.cos(theta/2)
    sin_theta = np.sin(theta/2)
    assert np.allclose(matrix[0, 0], cos_theta)
    expected_matrix = np.array([[cos_theta, -1j * sin_theta], [-1j * sin_theta, cos_theta]])
    return (np.allclose(matrix, expected_matrix, atol=theta_tolerance), np.real(theta).astype(np.float64))

def is_rccx_gate(matrix, theta_tolerance=1e-12):
    if matrix.shape != (2, 2) or not np.allclose(matrix[0, 0], 0): 
        return False
    # print(matrix)
    return np.allclose(matrix, np.asarray([[0, -1], [1, 0]])) or np.allclose(matrix, np.asarray([[0, 1], [-1, 0]]))

def is_zgate(matrix, theta_tolerance=1e-12):
    if matrix.shape != (2, 2) or not np.allclose(matrix[0, 1], 0): 
        return False
    return np.allclose(matrix, np.asarray([[-1, 0], [0, 1]])) or np.allclose(matrix, np.asarray([[1, 0], [0, -1]]))

# Helper function to convert a Cirq gate to a Qiskit gate
def _convert_cirq_gate_to_qiskit(cirq_gate):
    """Convert Cirq gate to Qiskit gate."""
    # print(cirq_gate)
    if cirq.unitary(cirq_gate).shape[0] == 2:
        if np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.X)):
            return Gate(name="x", num_qubits=1, params=[])
        elif np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.Y)):
            return Gate(name="y", num_qubits=1, params=[])
        elif np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.Z)):
            return Gate(name="z", num_qubits=1, params=[])
        elif np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.H)):
            return Gate(name="h", num_qubits=1, params=[])
        else:
            return UGate(cirq.unitary(cirq_gate))
    elif cirq.unitary(cirq_gate).shape[0] == 4:
        if np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.CX)):
            return Gate(name="cx", num_qubits=2, params=[])
        elif np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.CZ)):
            return Gate(name="cz", num_qubits=2, params=[])
        else:       
            return UGate(cirq.unitary(cirq_gate))
    elif cirq.unitary(cirq_gate).shape[0] == 8:
        # You can add more gate conversions as needed
        # print(cirq_gate)
        if np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.CCX)):
            #print("CCX gate")
            return Gate(name="ccx", num_qubits=3, params=[])
        elif is_zgate(cirq.unitary(cirq_gate)[-2:, -2:]):
            qc = QuantumCircuit(3)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
            qc.append(CCZGate(), [0, 1, 2]) # Do a CCX gate if we do not want improvements 
            # return Gate(name="ccz", num_qubits=3, params=[])
            check = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=3)
            # print(check)
            return [(x.operation, x.qubits) for x in check.data]
            # return UGate(cirq.unitary(cirq_gate))    
        elif is_rccx_gate(cirq.unitary(cirq_gate)[-2:, -2:]):
            qc = QuantumCircuit(3)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
            qc.append(RCCXGate(), [0, 1, 2]) # Do a CCX gate if we do not want improvements 
            check = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=3)
            return [(x.operation, x.qubits) for x in check.data]
            return UGate(cirq.unitary(cirq_gate))
        elif is_ry_gate(cirq.unitary(cirq_gate)[-2:, -2:])[0]:
            y_gate_theta = is_ry_gate(cirq.unitary(cirq_gate)[-2:, -2:])[1]
            # gate = RYGate(theta=).control(num_ctrl_qubits=2)
            qc = QuantumCircuit(3)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
            g1 = RYGate(theta=y_gate_theta/2).control(num_ctrl_qubits=1)
            g2 = RYGate(theta=-1*y_gate_theta/2).control(num_ctrl_qubits=1)
            g3 = RYGate(theta=1*y_gate_theta/2).control(num_ctrl_qubits=1)
            qc.append(g1, [1, 2])
            qc.append(CXGate(), [0, 1])
            qc.append(g2, [1, 2])
            qc.append(CXGate(), [0, 1])
            qc.append(g3, [0, 2])
            check = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=3)
            # print("RCCY gate ", len(check.data))
            return [(x.operation, x.qubits) for x in check.data]
            # return UGate(cirq.unitary(cirq_gate))
        elif is_rx_gate(cirq.unitary(cirq_gate)[-2:, -2:])[0]:
            gate = RXGate(theta=is_rx_gate(cirq.unitary(cirq_gate)[-2:, -2:])[1]).control(num_ctrl_qubits=2)
            qc = QuantumCircuit(3)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
            qc.append(gate, [0, 1, 2])
            check = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=3)
            return UGate(cirq.unitary(cirq_gate))
        else:
            qc = QuantumCircuit(3)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
            gate = UGate(cirq.unitary(cirq_gate)[-2:, -2:]).control(num_ctrl_qubits=2)
            qc.append(gate, [0, 1, 2])
            check = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=3)
            return UGate(cirq.unitary(cirq_gate))
    else:
        raise ValueError(f"Unsupported Cirq gate type: {type(cirq_gate)}")

# Function to decompose Qiskit circuit to Clifford+T
def decompose_to_clifford_t(qiskit_circuit):
    """Decompose the Qiskit circuit to Clifford+T gateset.""" # 
    transpiled = transpile(qiskit_circuit, basis_gates=['h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3', 'cx'], optimization_level=0)
    return transpiled

# Function to convert Qiskit circuit back to Cirq
def qiskit_to_cirq(qiskit_circuit):
    """Convert Qiskit circuit back to Cirq."""
    cirq_circuit = cirq.Circuit()
    cirq_qubits = [cirq.LineQid(i, dimension=2) for i in range(len(qiskit_circuit.qubits))]
    for gate in qiskit_circuit.data:
        qubits = gate.qubits
        cirq_gate = _convert_qiskit_gate_to_cirq(gate[0])
        if isinstance(qubits, list) or isinstance(qubits, tuple):
            cirq_circuit.append(cirq_gate.on(*[cirq_qubits[q.index] for q in qubits]))
        else:
            cirq_circuit.append(cirq_gate.on(cirq_qubits[qubits.index]))
    return cirq_circuit

# Helper function to convert a Qiskit gate to a Cirq gate
def _convert_qiskit_gate_to_cirq(qiskit_gate):
    """Convert Qiskit gate to Cirq gate."""
    if qiskit_gate.name == 'x':
        return cirq.X
    elif qiskit_gate.name == 'y':
        return cirq.Y
    elif qiskit_gate.name == 'z':
        return cirq.Z
    elif qiskit_gate.name == 'h':
        return cirq.H
    elif qiskit_gate.name == 'cx':
        return cirq.CX
    elif qiskit_gate.name == 'cz':
        return cirq.CZ
    elif qiskit_gate.name == 't':
        return cirq.T
    elif qiskit_gate.name == 'tdag':
        return cirq.Tdag
    elif qiskit_gate.name == 's':
        return cirq.S
    elif qiskit_gate.name == 'sdag':
        return cirq.Sdag
    elif qiskit_gate.name == 'u3' or qiskit_gate.name == 'u' or qiskit_gate.name == 'u1' or qiskit_gate.name == 'u2':
        # print(qiskit_gate)
        return cirq.MatrixGate(qiskit_gate.to_matrix())
        # theta=qiskit_gate.data[0][0].params[0], phi=qiskit_gate.data[0][0].params[1], lambda_=qiskit_gate.data[0][0].params[2]) 
        # You can add more gate conversions as needed
    else:
        print(qiskit_gate)
        raise ValueError(f"Unsupported Qiskit gate type: {qiskit_gate.name}")

# Main function
def replace_qubit_circuit(cirq_circuit):
    """Convert Cirq circuit to Qiskit, decompose it to Clifford+T, and convert it back to Cirq."""
    # Step 1: Convert Cirq to Qiskit
    qiskit_circuit = cirq_to_qiskit(cirq_circuit)
    print("1 STEPS DONE")

    # Step 2: Decompose the Qiskit circuit to Clifford + T
    decomposed_qiskit_circuit = decompose_to_clifford_t(qiskit_circuit)
    print("2 STEPS DONE")

    # Step 3: Convert the decomposed Qiskit circuit back to Cirq
    final_cirq_circuit = qiskit_to_cirq(decomposed_qiskit_circuit)
    print("ALL STEPS DONE")
    return final_cirq_circuit

# def replace_qubit_circuit(circuit, qubits=None):
#     final_circ = cirq.Circuit()
#     xgate = np.asarray([[0, 1],
#                         [1, 0]], dtype=np.complex128)
#     # Define known Clifford+T matrices
#     T_matrix = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]])
#     T_dag_matrix = np.array([[1, 0], [0, np.exp(-1j * np.pi / 4)]])
#     H_matrix = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
#     S_matrix = np.array([[1, 0], [0, 1j]])
#     S_dag_matrix = np.array([[1, 0], [0, -1j]])
#     for m_index, moment in enumerate(circuit):
#         for op in moment:
#             print(op)
#             if cirq.unitary(op).shape[0] == 2:
#                 q = op.qubits[0]
#                 try:
#                     U = cirq.unitary(op.gate)
#                 except:
#                     raise ValueError(f"Cannot extract unitary for {op.gate}")
#                 # Check if the unitary is already in the Clifford+T set
#                 if np.allclose(U, np.eye(2)):  # Identity
#                     continue 
#                 elif np.allclose(U, xgate):  # Identity
#                     final_circ.append(cirq.X(q))    
#                 elif np.allclose(U, T_matrix):  # T gate
#                     final_circ.append(cirq.T(q))
#                 elif np.allclose(U, T_dag_matrix):  # T† gate
#                     final_circ.append(cirq.T(q)**-1)
#                 elif np.allclose(U, H_matrix):  # H gate (Clifford)
#                     final_circ.append(cirq.H(q))
#                 elif np.allclose(U, S_matrix):  # S gate (Clifford)
#                     final_circ.append(cirq.S(q))
#                 elif np.allclose(U, S_dag_matrix):  # S† gate
#                     final_circ.append(cirq.S(q)**-1)
#                 # Now check if the unitary matches Rx, Ry, Rz, or U3
#                 elif np.allclose(U, H_matrix @ T_matrix @ H_matrix):  # Rx(π/4) = H T H
#                     final_circ.extend([cirq.H(q), cirq.T(q), cirq.H(q)])
#                 elif np.allclose(U, S_matrix @ H_matrix @ T_matrix @ H_matrix @ S_dag_matrix):  # Ry(π/4)
#                     final_circ.extend([cirq.S(q), cirq.H(q), cirq.T(q), cirq.H(q), cirq.S(q)**-1])
#                 elif np.allclose(U, [[1, 0], [0, np.exp(1j * np.pi / 2)]]):  # Rz(π/2) = S
#                     final_circ.append(cirq.S(q))
#                 elif np.allclose(U, [[1, 0], [0, np.exp(-1j * np.pi / 2)]]):  # Rz(-π/2) = S†
#                     final_circ.append(cirq.S(q)**-1)
#                 elif np.allclose(U, [[1, 0], [0, np.exp(1j * np.pi)]]):  # Rz(π) = Z
#                     final_circ.append(cirq.Z(q))
#                 elif np.allclose(U, [[1, 0], [0, np.exp(-1j * np.pi)]]):  # Rz(-π) = Z†
#                     final_circ.append(cirq.Z(q)**-1)
#                 else:
#                     check_u3 = is_u3_matrix(U)
#                     if check_u3[0]:
#                         # theta, phi, lamb = check_u3[1], check_u3[2], check_u3[3]
#                         final_circ.append(op)
#                         # ([cirq.Rz(lamb)(q), cirq.H(q), cirq.Rz(theta)(q), cirq.H(q), cirq.Rz(phi)(q)])
#                     else:
#                         raise ValueError(f"Unsupported gate (not in Clifford+T or Rx, Ry, Rz, U3): {op.gate}, matrix:\n{U}")
#                 # final_circ.append(op)
#             elif cirq.unitary(op).shape[0] == 4:
#                 control = op.qubits[0]
#                 target = op.qubits[1]
#                 matrix =  cirq.unitary(op)[-2:, -2:]
#                 if not np.allclose(matrix, xgate):
#                     print("not a cnot")
#                 final_circ.append(cirq.MatrixGate(matrix).on(target).controlled_by(control))
#             elif cirq.unitary(op).shape[0] == 8:
#                 matrix = cirq.unitary(op)[-2:, -2:]                
#                 if np.allclose(matrix, xgate):
#                     final_circ.append(toffoli_circ(op.qubits))
#                 else:
#                     print(cirq.unitary(op))
#                     raise Exception("Not supported for decomposing yet 3 qubits")
#             else:
#                 print("last operation was ")
#                 raise Exception("Not supported for decomposing yet")
#     return final_circ

if __name__ == "__main__":
    circ1 = cirq.Circuit()
    all_qubits = cirq.LineQid.range(8, dimension=2)
    circ1.append(QubitCNOTGate().on(all_qubits[0], all_qubits[1]))
    circ1.append(QubitCNOTGate().on(all_qubits[2], all_qubits[3]))
    circ1.append(toffoli_circ([all_qubits[1], all_qubits[3], all_qubits[6]]))
    circ1.append(QubitXGate().on(all_qubits[4]))
    circ1.append(toffoli_circ([all_qubits[6], all_qubits[4], all_qubits[5]]))
    circ1.append(QubitXGate().on(all_qubits[4]))
    circ1.append(toffoli_circ([all_qubits[1], all_qubits[3], all_qubits[6]]))
    circ1.append(QubitCNOTGate().on(all_qubits[0], all_qubits[1]))
    circ1.append(QubitCNOTGate().on(all_qubits[2], all_qubits[3]))

    circ1.append(QubitXGate().on(all_qubits[0]))
    circ1.append(QubitXGate().on(all_qubits[1]))
    circ1.append(QubitCNOTGate().on(all_qubits[3], all_qubits[4]))
    circ1.append(toffoli_circ([all_qubits[0], all_qubits[1], all_qubits[6]]))
    circ1.append(toffoli_circ([all_qubits[2], all_qubits[4], all_qubits[7]]))
    circ1.append(toffoli_circ([all_qubits[6], all_qubits[7], all_qubits[5]]))
    circ1.append(toffoli_circ([all_qubits[0], all_qubits[1], all_qubits[6]]))
    circ1.append(toffoli_circ([all_qubits[2], all_qubits[4], all_qubits[7]]))
    circ1.append(QubitCNOTGate().on(all_qubits[3], all_qubits[4]))
    circ1.append(QubitXGate().on(all_qubits[0]))
    circ1.append(QubitXGate().on(all_qubits[1]))


    circ1.append(QubitXGate().on(all_qubits[2]))
    circ1.append(QubitXGate().on(all_qubits[3]))
    circ1.append(QubitCNOTGate().on(all_qubits[1], all_qubits[4]))
    circ1.append(toffoli_circ([all_qubits[0], all_qubits[2], all_qubits[6]]))
    circ1.append(toffoli_circ([all_qubits[3], all_qubits[4], all_qubits[7]]))
    circ1.append(toffoli_circ([all_qubits[6], all_qubits[7], all_qubits[5]]))
    circ1.append(toffoli_circ([all_qubits[0], all_qubits[2], all_qubits[6]]))
    circ1.append(toffoli_circ([all_qubits[3], all_qubits[4], all_qubits[7]]))
    circ1.append(QubitCNOTGate().on(all_qubits[1], all_qubits[4]))
    circ1.append(QubitXGate().on(all_qubits[2]))
    circ1.append(QubitXGate().on(all_qubits[3]))

    circ1.append(QubitXGate().on(all_qubits[0]))
    circ1.append(QubitXGate().on(all_qubits[2]))
    circ1.append(QubitCNOTGate().on(all_qubits[1], all_qubits[3]))
    circ1.append(toffoli_circ([all_qubits[0], all_qubits[2], all_qubits[6]]))
    circ1.append(toffoli_circ([all_qubits[3], all_qubits[4], all_qubits[7]]))
    circ1.append(toffoli_circ([all_qubits[6], all_qubits[7], all_qubits[5]]))
    circ1.append(toffoli_circ([all_qubits[0], all_qubits[2], all_qubits[6]]))
    circ1.append(toffoli_circ([all_qubits[3], all_qubits[4], all_qubits[7]]))
    circ1.append(QubitCNOTGate().on(all_qubits[1], all_qubits[3]))
    circ1.append(QubitXGate().on(all_qubits[0]))
    circ1.append(QubitXGate().on(all_qubits[2]))
    # print(circ1)


    circ2 = cirq.Circuit()
    all_qubits = cirq.LineQid.range(8, dimension=2)
    circ2.append(toffoli_circ([all_qubits[0], all_qubits[1], all_qubits[5]]))
    circ2.append(QubitCNOTGate().on(all_qubits[0], all_qubits[1]))
    circ2.append(toffoli_circ([all_qubits[2], all_qubits[5], all_qubits[6]]))
    
    circ2.append(toffoli_circ([all_qubits[1], all_qubits[2], all_qubits[5]]))
    circ2.append(QubitCNOTGate().on(all_qubits[1], all_qubits[2]))
    circ2.append(toffoli_circ([all_qubits[3], all_qubits[5], all_qubits[6]]))
    
    circ2.append(toffoli_circ([all_qubits[2], all_qubits[3], all_qubits[5]]))
    circ2.append(QubitCNOTGate().on(all_qubits[2], all_qubits[3]))
    circ2.append(toffoli_circ([all_qubits[4], all_qubits[5], all_qubits[6]]))
    
    circ2.append(toffoli_circ([all_qubits[3], all_qubits[4], all_qubits[5]]))
    circ2.append(QubitCNOTGate().on(all_qubits[3], all_qubits[4]))
    circ2.append(QubitCNOTGate().on(all_qubits[6], all_qubits[5]))

    circ2.append(QubitCNOTGate().on(all_qubits[5], all_qubits[7]))

    circ2.append(QubitCNOTGate().on(all_qubits[6], all_qubits[5]))
    circ2.append(QubitCNOTGate().on(all_qubits[3], all_qubits[4]))
    circ2.append(toffoli_circ([all_qubits[3], all_qubits[4], all_qubits[5]]))
    
    circ2.append(toffoli_circ([all_qubits[3], all_qubits[5], all_qubits[6]]))
    circ2.append(QubitCNOTGate().on(all_qubits[1], all_qubits[2]))
    circ2.append(toffoli_circ([all_qubits[1], all_qubits[2], all_qubits[5]]))

    circ2.append(toffoli_circ([all_qubits[4], all_qubits[5], all_qubits[6]]))    
    circ2.append(QubitCNOTGate().on(all_qubits[2], all_qubits[3]))
    circ2.append(toffoli_circ([all_qubits[2], all_qubits[3], all_qubits[5]]))

    circ2.append(toffoli_circ([all_qubits[2], all_qubits[5], all_qubits[6]]))
    circ2.append(QubitCNOTGate().on(all_qubits[0], all_qubits[1]))
    circ2.append(toffoli_circ([all_qubits[0], all_qubits[1], all_qubits[5]]))
    print(circ2)
    depth_counter(circ1, all_qubits)

    depth_counter(circ2, all_qubits)