import cirq
from collections import defaultdict
import networkx as nx
from src.toffoli import *
from decomp.src.circuit_lowering import *
# from decomp.src.circuit_decomposition import *
from rewriter.simple_circuit_rewrites import *
import numpy as np
import time
import matplotlib.pyplot as plt
import collections
import time
import warnings
import rustworkx as rx
from src.error_models import CurrentSuperconductingQCErrors, ModifiedCurrSuperconductingQCErrors,\
  FutureSuperconductingQCErrors, ModifiedFutureSuperconductingQCErrors, FutureSuperconductingQCErrorsBetterT1AndGates,\
        DressedQutritErrors, BareQutritErrors, FutureSuperconductingQCErrorsBetterT1, FutureSuperconductingQCErrorsBetterGates
from src.unitary_operations import get_random_state, apply_unitaries_to_circuit
from decomp.src.run_qubit_decomposition import run_qubit_decomposition
from decomp.src.gate_decomposition import *
from rewriter.utils.check_subspace_equality import *
from rewriter.src.nodeset import LabelledNode
from rewriter.src.node_directory import NodeDirectory
from rewriter.src.gate_inverses import TrackingDict
from rewriter.utils.graph_utils import *
from replacement_lib.replace_qubit_circ import replace_qubit_circuit


class GraphInfo:
  def __init__(self):
    self.node_info = NodeDirectory()
    self.edge_info = defaultdict(dict)
    self.dag = None

  def check_gate_type(self, node):
    return self.node_info[node].get_gate_type()
  
  def create_node_directory(self, dag):
    self.dag = dag
    for u in nx.topological_sort(dag):
      if "label" in dag.nodes()[u] and dag.nodes()[u]["label"] != "":
        self.node_info[u] = LabelledNode(node_info=dag.nodes()[u])
    s = 0
    for u, v, _ in dag.edges:
      s += 1
      assert isinstance(_, int)
      a = dag.get_edge_data(u, v)[_]
      for k in a.keys():
        self.edge_info[(u, v, _)][k] = str(a[k]).replace('"', '')
  
  def get_all_qubits(self):
    qubit_nodes = []
    for nodes in self.node_info:
      if self.node_info[nodes].is_root():
        qubit_nodes.append(nodes)
    if len(qubit_nodes) == 0:
      raise Exception("This is a manual circuit")
    return qubit_nodes


def post_rewrite_compile(circuit, qutrits, run_sim=False, bench=False, if_verify=True, compile_mode=""):
    if "qubit" in compile_mode:
      single_gate_depth2, double_gate_depth2, one_qubit_gate_count, two_qubit_gate_count = depth_counter(replace_qubit_circuit(circuit), qutrits)
      return single_gate_depth2, two_qubit_gate_count, one_qubit_gate_count, {}

    circ_cs_decomp = cirq.Circuit()
    start = time.time() 
    # Gate simplify was set to True earlier
    single_gate_depth2, double_gate_depth2, _, _ = depth_counter(circuit, qutrits)
    circ_cs_decomp.append(return_decomposed_circuit(circuit, qutrits, dim=3, dtype=np.complex128,
                                                    debug_mode_=False, gate_simplify=False, checking_blocked=True))
    end = time.time()
    single_gate_depth, double_gate_depth, _, _ = depth_counter(circ_cs_decomp, qutrits)
    if not bench:
      for i in range(3):
        if i == 0:
          final_circ = cirq.Circuit()
          final_circ.append(loop_through_rewrites(circ_cs_decomp, qutrits, checking_blocked=True))
        else:
          final_circ2 = cirq.Circuit()
          final_circ2.append(loop_through_rewrites(final_circ, qutrits, checking_blocked=True))
    else:
      final_circ2 = cirq.Circuit()
      final_circ2.append(circ_cs_decomp)
    T2 = time.time()

    print(final_circ2)
    single_gate_depth2, double_gate_depth2, one_qubit_gate_count, two_qubit_gate_count = depth_counter(final_circ2, qutrits)
    print("COMPARISON: ", "original circuit ", single_gate_depth, double_gate_depth,
          "basic rewrites", single_gate_depth2, double_gate_depth2)
    if run_sim:
      example_noise_model1 = CurrentSuperconductingQCErrors()
      example_noise_model2 = DressedQutritErrors()
      example_noise_model3 = FutureSuperconductingQCErrors()
      example_noise_model4 = ModifiedCurrSuperconductingQCErrors()
      example_noise_model5 = ModifiedFutureSuperconductingQCErrors()
      example_noise_model6 = FutureSuperconductingQCErrorsBetterT1()
      example_noise_model7 = FutureSuperconductingQCErrorsBetterGates()
      example_noise_model8 = FutureSuperconductingQCErrorsBetterT1AndGates()

      results1 = []
      results2 = []
      results3 = []
      results4 = []
      results5 = []
      results6 = []
      results7 = []
      results8 = []

      num = len(qutrits)
      for i in range(1000):
        random_state = get_random_state(num)
        ideal_state = apply_unitaries_to_circuit(final_circ2, random_state)
        noisy_state = apply_unitaries_to_circuit(final_circ2, random_state, noise_model=example_noise_model1)
        result = np.linalg.norm(np.vdot(ideal_state, noisy_state)) ** 2
        results1.append(result)
        noisy_state = apply_unitaries_to_circuit(final_circ2, random_state, noise_model=example_noise_model2)
        result = np.linalg.norm(np.vdot(ideal_state, noisy_state)) ** 2                                  
        results2.append(result)
        noisy_state = apply_unitaries_to_circuit(final_circ2, random_state, noise_model=example_noise_model3)
        result = np.linalg.norm(np.vdot(ideal_state, noisy_state)) ** 2
        results3.append(result)
        noisy_state = apply_unitaries_to_circuit(final_circ2, random_state, noise_model=example_noise_model4)
        result = np.linalg.norm(np.vdot(ideal_state, noisy_state)) ** 2
        results4.append(result)
        noisy_state = apply_unitaries_to_circuit(final_circ2, random_state, noise_model=example_noise_model5)
        result = np.linalg.norm(np.vdot(ideal_state, noisy_state)) ** 2
        results5.append(result)
        noisy_state = apply_unitaries_to_circuit(final_circ2, random_state, noise_model=example_noise_model6)
        result = np.linalg.norm(np.vdot(ideal_state, noisy_state)) ** 2
        results6.append(result)
        noisy_state = apply_unitaries_to_circuit(final_circ2, random_state, noise_model=example_noise_model7)
        result = np.linalg.norm(np.vdot(ideal_state, noisy_state)) ** 2
        results7.append(result)
        noisy_state = apply_unitaries_to_circuit(final_circ2, random_state, noise_model=example_noise_model8)
        result = np.linalg.norm(np.vdot(ideal_state, noisy_state)) ** 2
        results8.append(result)        
      results1 = sum(results1)/len(results1)
      results2 = sum(results2)/len(results2)
      results3 = sum(results3)/len(results3)
      results4 = sum(results4)/len(results4)
      results5 = sum(results5)/len(results5)
      results6 = sum(results6)/len(results6)
      results7 = sum(results7)/len(results7)
      results8 = sum(results8)/len(results8)
      print("Rewritten circuit simulation results ", results1, results2, results3, "modified", results4, results5, "future", results6, results7, results8)
      return single_gate_depth2, two_qubit_gate_count, one_qubit_gate_count, \
              {"current_sc":results1, "dressed_ti":results2, "removed_2_scaling":results4, "all_scaled":results6,
                      "others are ": ["Rewritten circuit simulation results ", results1, results2, results3,
                          "modified", results4, results5, "future", results6, results7, results8]} 
    return single_gate_depth2, two_qubit_gate_count, one_qubit_gate_count, {}
