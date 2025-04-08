import numpy as np
import cirq
import rustworkx as rx
import pydot
import os

def is_binary_mapping(matrix):
    assert matrix.shape[0] == 3
    if np.allclose(matrix[2,2] * matrix[2,2].conjugate(), np.eye(1)):
        return True
    return False

def if_mat_preserves_binary(matrix):
  for i in range(matrix.shape[0]):
    if not np.allclose(matrix[i, i] * matrix[i, i].conjugate(), np.eye(1)):
      return False
  return True

def check_matrix_phase(matrix):
   eigenvalues, eigenvectors = np.linalg.eig(matrix)
   common_phase = np.angle(eigenvalues)
   theta = np.sum(common_phase)
   while theta > 3 * np.pi:
      theta = theta - 3 * np.pi
   while theta < 0:
      theta = theta + 3 * np.pi
   if np.allclose(theta, 0):
    return matrix, theta
   matrix_2 =  matrix * np.exp(1j * (theta - np.pi))
   return matrix_2, theta - np.pi

def get_circuit_unitaries(circuit, dtype=np.complex128, dimension=3):
    dim = dimension
    qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
    elems_avail = []
    for a in qd_map.keys():
      qd_map[a] = a.x
      elems_avail.append(a.x)
    elems_avail = sorted(elems_avail)
    for a in qd_map.keys():
      qd_map[a] = elems_avail.index(a.x)
    
    num_qudits = len(qd_map.keys())
    n = num_qudits
    new_state = np.zeros((dim,)*n, dtype=np.complex128)
    unitary_matrix = np.zeros((dim ** n, dim ** n), dtype=np.complex128)
    j = 0
    for i in range(dim ** n) :
        temp_state = np.zeros(dim ** n, dtype=np.complex128)
        temp_state[i] = 1.0
        state = temp_state.reshape(new_state.shape)
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
    assert np.allclose(np.eye(len(unitary_matrix)), unitary_matrix.dot(unitary_matrix.T.conj()))
    return unitary_matrix

def build_rx_graph(circuit, qubits):
  dag = rx.PyDiGraph()
  in_nodes = dag.add_nodes_from(qubits)
  out_nodes = dag.add_nodes_from(qubits)
  parent_node = {}
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

def remove_broadcast_nodes(G):
    g = G.copy()
    counter = 0
    for a in g:
        g0 = g.copy()
        node_i = g.nodes[a]
        if len(node_i.keys()) > 0 and "C" not in node_i["label"]:
            counter += 1
            for preds in g0.predecessors(a):
                for succs in g0.neighbors(a):
                    g0.add_edge(preds, succs, **(g0.get_edge_data(a, succs)[0]))
            g0.remove_node(a)
        g = g0
    for a in list(nx.topological_sort(g)):
        for preds in g.predecessors(a):
            pass
    return g

def is_in(t, l):
  if isinstance(t, list) and isinstance(l, str):
    return ''.join(t) in l
  if isinstance(t, list) and isinstance(l, list):
    return set(t) <= set(l)
  if isinstance(t, set) and isinstance(l, set):
    return t <= l
  if isinstance(t, str) and isinstance(l, str):
    if t in l:
      return True
    else:
      return False
  print("case has not been handled yet for types ", type(t), type(l))
  assert False

def is_equal(t, l):
  if isinstance(t, list) and isinstance(l, list):
    t = ''.join(str(t))
    l = ''.join(str(l))
    return t == l  
  if isinstance(t, str) and isinstance(l, str):
    return t == l
  print("case has not been handled yet for types ", type(t), type(l))
  assert False