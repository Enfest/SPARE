import cirq
import math
import numpy as np
from itertools import product
from toffoli import *
from qutrit_decompositions import *
from controlled_ternary_gate import *
from verifying_gates import *
from gate_decomposition import *
from qutrit_decomposition_kak import *
from cirq_depth_counter import * 
from check_subspace_equality import *
from rewritten_circuits import *
from controlled_ternary_gate import * 
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from check_subspace_equality import return_qubit_circuit
import rustworkx as rx
import pydot
import os

def build_rx_graph(circuit, qubits):
  dag = rx.PyDiGraph()
  in_nodes = dag.add_nodes_from(qubits)
  out_nodes = dag.add_nodes_from(qubits)
  parent_node = {}
  # print(in_nodes)
  qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
  min_val = 100000000000
  for a in qd_map.keys():
    qd_map[a] = a.x
    min_val = min(min_val, a.x)
  for a in qd_map.keys():
    qd_map[a] = a.x - min_val
  for q in circuit.all_qubits():
    parent_node[qd_map[q]] = in_nodes[qd_map[q]]
  for m_index, moment in enumerate(circuit):
    for ops in moment:
      if len(ops.qubits) == 1:
        parent_node[qd_map[ops.qubits[0]]] = dag.add_child(parent_node[qd_map[ops.qubits[0]]], ops, str(qd_map[ops.qubits[0]]))
      else:
        gate = dag.add_child(parent_node[qd_map[ops.qubits[-1]]], ops, str(qd_map[ops.qubits[-1]]))
        parent_node[qd_map[ops.qubits[-1]]] = gate
        for i in range(len(ops.qubits[:-1])):
          dag.add_edge(parent_node[qd_map[ops.qubits[i]]], gate, str(qd_map[ops.qubits[i]]))
          parent_node[qd_map[ops.qubits[i]]] = gate
  for nodes in parent_node.keys():
    dag.add_edge(parent_node[nodes], out_nodes[nodes], str(nodes))
  return dag

def get_and_print_dag(circuit, qubits, name=""):
    dag = build_rx_graph(circuit, qubits)
    dot_str = dag.to_dot(node_attr=lambda node: {"label": str(node)}, edge_attr=lambda edge: {"label": str(edge)})
    dot = pydot.graph_from_dot_data(dot_str)[0]
    # Save the DOT graph into a DOT file
    dot_file_path = name + "_graph.dot"
    with open(dot_file_path, 'w') as dot_file:
        dot_file.write(dot.to_string())
    tmp_path = os.path.join(name + 'dag.png')
    dot.write_png(tmp_path)
    return

def is_cnot_gate(matrix):
    if matrix.shape[0] == 4:
        return True
    else:
        return False

def is_toffoli_gate(matrix):
    if matrix.shape[0] == 8:
        return True
    return False

def is_identity_gate(matrix):
    if matrix.shape[0] == 2:
        return True
    return False

def create_qutrit_from_qubit(circuit, size):
    new_circuit = cirq.Circuit()
    qutrits = cirq.LineQid.range(size, dimension=3)
    qd_map = {}
    for qubits in circuit.all_qubits():
        qd_map[qubits] = qubits.x
    
    qd2 = {}
    for qs in circuit.all_qubits():
        qd2[qs] = qs.x

    print("CHECK")
    print(qd2)
    for m_idx, moment in enumerate(circuit):
        for op in moment:
            print("op is ", is_identity_gate(cirq.unitary(op)), is_toffoli_gate(cirq.unitary(op)), is_cnot_gate(cirq.unitary(op)))
            if is_identity_gate(cirq.unitary(op)):
                pass
            elif is_toffoli_gate(cirq.unitary(op)):
                print(op)
                qs = op.qubits
                print(qs, qs[-1], qd_map)
                print(qs[-1])
                print(qd2[qs[-1]])
                
                t = qutrits[qd2[qs[-1]]]
                ctrl = [qutrits[qd2[x]] for x in qs[:-1]]
                new_circuit.append(QutritX12Gate().on(ctrl[1]).controlled_by(ctrl[0]))
                new_circuit.append(QutritX12Gate().on(ctrl[1]))
                new_circuit.append(QutritX01Gate().on(t).controlled_by(* ctrl))
                new_circuit.append(QutritX12Gate().on(ctrl[1]))
                new_circuit.append(QutritX12Gate().on(ctrl[1]).controlled_by(ctrl[0]))
            elif is_cnot_gate(cirq.unitary(op)):
                # print("1 val ", op)
                qs = op.qubits
                t = qutrits[qd2[qs[-1]]]
                ctrl = [qutrits[qd2[x]] for x in qs[:-1]]
                new_circuit.append(QutritX01Gate().on(t).controlled_by(ctrl[0]))
            elif cirq.matrix(op).shape[0] == 2:
                matrix = np.eye(3, dtype=np.complex128)
                matrix[:2, :2] = cirq.matrix(op)
                new_circuit.append(cirq.MatrixGate(matrix, qid_shape=[3]).on(qutrits[qs[-1]]))
            elif cirq.matrix(op).shape[0] == 4:
                matrix = np.eye(3, dtype=np.complex128)
                if not np.allclose(cirq.matrix(op)[:2, :2], np.eye(2)):
                        raise Exception()
                matrix[:2, :2] = cirq.matrix(op)[2:, 2:]
                new_circuit.append(cirq.MatrixGate(matrix, qid_shape=[3])).on(qutrits[qs[-1]])
            else:
                assert False

    return new_circuit, qutrits

def get_gate(matrix):
    if matrix.shape[0] % 3 == 0:
        d = 3
        matx2 = matrix[(matrix.shape[0] - d) // 2: (matrix.shape[0] + d) // 2,
                       (matrix.shape[0] - d) // 2: (matrix.shape[0] + d) // 2]
    else:
        d = 2
    if matrix.shape == 9:
        return True
    if matrix.shape == 27:
        return False

def create_qutrit_from_qubit_(circuit, size):
    new_circuit = cirq.Circuit()
    qutrits = cirq.LineQid(size, dimension=3)
    qd_map = {}
    for qubits in circuit.all_qubits():
        qd_map[qubits.x] = qubits
    
    for m_idx, moment in enumerate(circuit):
        for op in moment:

            if is_toffoli_gate(cirq.unitary(op)):
                qs = op.qubits()
                t = qutrits[qd_map[qs[-1]]]
                ctrl = [qutrits[qd_map[x]] for x in qs[:-1]]
                new_circuit.append(QutritX12Gate().on(ctrl[1]).controlled_by(ctrl[0]))
                new_circuit.append(QutritX01Gate().on(t).controlled_by(* ctrl))
                new_circuit.append(QutritX12Gate().on(ctrl[1]).controlled_by(ctrl[0]))
            
            elif is_cnot_gate:
                qs = op.qubits()
                t = qutrits[qd_map[qs[-1]]]
                ctrl = [qutrits[qd_map[x]] for x in qs[:-1]]
                new_circuit.append(QutritX01Gate().on(t).controlled_by(ctrl[0]))
