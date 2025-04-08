# Helper function to convert a Cirq gate to a Qiskit gate
def _convert_cirq_gate_to_qiskit(cirq_gate):
    """Convert Cirq gate to Qiskit gate."""
    if cirq.unitary(cirq_gate).shape[0] == 2:
            # if np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.X)):
            #     return Gate(name="x", num_qubits=1, params=[])
            # elif np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.Y)):
            #     return Gate(name="y", num_qubits=1, params=[])
            # elif np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.Z)):
            #     return Gate(name="z", num_qubits=1, params=[])
            # elif np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.H)):
            #     return Gate(name="h", num_qubits=1, params=[])
            # else:
        return UGate(cirq.unitary(cirq_gate))
    elif cirq.unitary(cirq_gate).shape[0] == 4:
        if np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.CX)):
            return Gate(name="cx", num_qubits=2, params=[])
        elif np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.CZ)):
            return Gate(name="cz", num_qubits=2, params=[])
        else:       
            qc = QuantumCircuit(2)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
            gate = UGate(cirq.unitary(cirq_gate)[-2:, -2:]).control(num_ctrl_qubits=1)
            qc.append(gate, [0, 1])
            check1 = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=0)
            qc = QuantumCircuit(2)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
            gate = UGate(cirq.unitary(cirq_gate))
            qc.append(gate, [0, 1])
            check2 = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=3)
            # if len(check1.data) < len(check2.data):
            #     return [(x[0], x.qubits) for x in check1.data]
            # else:
            #     return [(x[0], x.qubits) for x in check2.data] 
            return UGate(cirq.unitary(cirq_gate))
    elif cirq.unitary(cirq_gate).shape[0] == 8:
        # You can add more gate conversions as needed
        # if np.allclose(cirq.unitary(cirq_gate), cirq.unitary(cirq.CCX)):
        #     return Gate(name="ccx", num_qubits=3, params=[])
        # elif is_ry_gate(cirq.unitary(cirq_gate)[-2:, -2:])[0]:
        #     gate = RYGate(theta=is_ry_gate(cirq.unitary(cirq_gate)[-2:, -2:])[1]).control(num_ctrl_qubits=2)
        #     qc = QuantumCircuit(3)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
        #     qc.append(gate, [0, 1, 2])
        #     check = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=3)
        #     return UGate(cirq.unitary(cirq_gate))
        # elif is_rx_gate(cirq.unitary(cirq_gate)[-2:, -2:])[0]:
        #     gate = RXGate(theta=is_rx_gate(cirq.unitary(cirq_gate)[-2:, -2:])[1]).control(num_ctrl_qubits=2)
        #     qc = QuantumCircuit(3)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
        #     qc.append(gate, [0, 1, 2])
        #     check = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=3)
        #     return UGate(cirq.unitary(cirq_gate))
        # else:
        qc = QuantumCircuit(3)  # Assume it's a 1-qubit UGate, adjust for multi-qubit
        gate = UGate(cirq.unitary(cirq_gate)[-2:, -2:]).control(num_ctrl_qubits=2)
        qc.append(gate, [0, 1, 2])
        check = transpile(qc, basis_gates=['cx', 'h', 't', 'tdag', 'x', 'y', 'x', 's', 'sdg', 'u3'], optimization_level=3)
        return UGate(cirq.unitary(cirq_gate))
    else:
        raise ValueError(f"Unsupported Cirq gate type: {type(cirq_gate)}")