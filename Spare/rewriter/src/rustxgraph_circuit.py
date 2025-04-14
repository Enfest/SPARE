from collections import defaultdict
import networkx as nx
from src.toffoli import *
from rewriter.simple_circuit_rewrites import *
from rewriter.utils.check_subspace_equality import *
from rewriter.src.nodeset import LabelledNode
from rewriter.src.node_directory import NodeDirectory
from rewriter.utils.graph_utils import *
from Spare.src.cirq_depth_counter import depth_counter
from Spare.replacement_lib.replace_qubit_circ import replace_qubit_circuit


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
    single_gate_depth2, double_gate_depth2, one_qubit_gate_count, two_qubit_gate_count = depth_counter(replace_qubit_circuit(circuit), qutrits)
    return single_gate_depth2, two_qubit_gate_count, one_qubit_gate_count, {}