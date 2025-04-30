import random
import networkx as nx
import time
import pydot
import json
import os
import argparse
from src.toffoli import *
from utils.check_subspace_equality import *
import rewriter.src.rustxgraph_circuit as rx_graph
from rewriter.simple_circuit_rewrites import *
from rewriter.src.rustxgraph_circuit import post_rewrite_compile
from rewriter.src.circuit_lifting import *
from rewriter.src.compile import RewriteEngine
from rewriter.src.rewrites import apply_rule_2
from rewriter.src.compile_basic import *
from rewriter.utils.passmanage import rewritePassesManager, rewritePassManage, passPair
import Spare.ancilla_convertor.circuit_convertor
from Spare.src.cirq_depth_counter import toffoli_estimator, depth_estimater


class CLIR(rx_graph.GraphInfo):
    def __init__(self, name, manual_lifting=True, dimension=3):
        super().__init__()
        self.dimension=dimension
        self.dot = pydot.graph_from_dot_file(name + ".dot")[0]
        networkx_graph = nx.nx_pydot.from_pydot(self.dot)
        self.create_node_directory(networkx_graph)
        if not manual_lifting:
          self.lift_dag_circuit()
        self.basic_copy = None
    
    def lift_dag_circuit(self):
      arr = []
      for nodes in self.node_info:
        if len(self.node_info[nodes].get_label()) == 0: 
           arr.append(nodes)
      for a in arr:
         del self.node_info[a]
      for nodes in self.node_info:
        if len(self.node_info[nodes].get_label()) != 0: 
          if not self.node_info[nodes].data_initialized:
            self.node_info[nodes].set_qubits(get_qubits(self.node_info[nodes]))
            self.node_info[nodes].set_matrix(get_matrix(self.node_info[nodes]))
            self.node_info[nodes].set_control(get_controls(self.node_info[nodes]))
          if len(list(self.dag.predecessors(nodes))) == 0 and self.node_info[nodes].get_matrix() is None:
            self.node_info[nodes].set_root()
          elif len(list(self.dag.predecessors(nodes))) == 0:
            print(nodes)
            raise Exception("This gate has no qubits feeding it") 
          else:
            pass
          if len(list(self.dag.successors(nodes))) == 0 and self.node_info[nodes].get_matrix() is None:
            self.node_info[nodes].set_sink()
          elif len(list(self.dag.successors(nodes))) == 0:
            raise Exception("This gate has no qubits outputting to")
          else:
            pass
      return self
    
    def create_basic_copy(self):
       self.old_graph = self.dag.copy()
       self.old_nodes = self.node_info.copy()
       self.old_edges = self.edge_info.copy()
       
    def return_object(self):
      self.dag = self.old_graph.copy()
      self.node_info = self.old_nodes.copy()
      self.edge_info = self.old_edges.copy()
      return self.basic_copy
    
    def create_gate_slices(self):
      dag = self.dag
      gateslices = RewriteEngine(self)
      return gateslices

speculative_history = []
def speculatively_apply_rewrite1(circuit_obj, rewrite_passes):
  global speculative_history
  circuit_copy = circuit_obj.deepcopy()
  c, _ = circuit_obj.build_circuit(breakdown="vvv")
  unitary1 = get_circuit_unitaries(c, dimension=circuit_obj.dimension)
  for re_pass in rewrite_passes:
    if len(re_pass) > 1:
      target_finding_func = getattr(circuit_copy, re_pass[0])
      target_opt_func = getattr(circuit_copy, re_pass[1])
      targets, corr_bookend_gates = target_finding_func()
      if len(targets) == 0:
        return circuit_obj
      selected_target = re_pass[2](targets)
      if targets.index(selected_target) not in speculative_history:
        speculative_history.append(targets.index(selected_target))
        selected_bookends = corr_bookend_gates[speculative_history[-1]]
        anySimplificationMade = target_opt_func(selected_target, selected_bookends)
        if not anySimplificationMade:
          '''Abort'''
          return circuit_obj
      else:
        return circuit_obj
    else:
      target_opt_func = getattr(circuit_copy, re_pass[0])
      target_opt_func(debug=True)
      verify_gateslices(circuit_copy, True, unitary1)
  '''Check if the circuit complexity has reduced with this'''
  circuit_copy = gateslices_simplify(circuit_copy)
  orig_estimate_value = depth_estimater(circuit_obj.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
  new_estimate_value = depth_estimater(circuit_copy.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
  old2gates = toffoli_estimator(circuit_obj.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
  new2gates = toffoli_estimator(circuit_copy.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
  if new_estimate_value < orig_estimate_value or new2gates < old2gates:
    speculative_history = [targets.index(selected_target)]
    return circuit_copy
  else:
    return circuit_obj

def compile_basic(name="temp", manual_lifting=True, if_base=False, json_data=None, if_verify=True, dimension=3):
    # Currently the starting point is the old CLIR object that would read in dot files
    info = CLIR(name, manual_lifting=manual_lifting, dimension=dimension)
    gateslices = info.create_gate_slices()
    gateslices.merge_greedily(dont_create_mappings=True)
    unitary1 = None
    if if_verify:
      unitary1 = get_circuit_unitaries(gateslices.build_circuit(breakdown="vvv")[0], dimension=gateslices.dimension)
    if if_base:
      return gateslices
    gateslices.control_rewrite(debug=False)
    gateslices = collapse_states(gateslices, if_verify, unitary=unitary1)
    gateslices = gateslices_simplify(gateslices, debug=False)
    gateslices.prune_null_slices()
    gateslices.optimizer_pass(debug=False)
    all_passes = [
      rewritePassManage([passPair("get_all_bookend_broadcasts", "modify_middle_bookends", lambda lst: lst[0] if lst else None, if_spec=True, times=10)]),
      rewritePassManage([passPair(None, "gateslices_simplify", None, if_spec=False, times=2)]),
      rewritePassManage([passPair("get_all_broadcastable_gates", "apply_rule_broadcasts", lambda lst: lst[0] if lst else None, if_spec=False, times=1), 
                    passPair(None, "optimizer_pass", lambda lst: lst[0] if lst else None, if_spec=False, times=1), 
                  ], if_spec=True, times=3),
      rewritePassManage([passPair("get_all_bookend_broadcasts", "modify_middle_bookends", lambda lst: lst[0] if lst else None, if_spec=True, times=10)]),
      rewritePassManage([passPair(None, "gateslices_simplify", None, if_spec=False, times=2)])
    ]
    gateslices = rewritePassesManager(all_passes).applyAllPasses(gateslices)
    gateslices.any_nodeset_preset(debug="length_simplified_orig1" in name)
    gateslices.serialize(name + ".json")
    return gateslices, unitary1

def compile_random(info=None, name="temp", manual_lifting=False, if_return_basic=False, seed=0,
                   json_data=None, itr=0, if_verify=True, compile_mode="", lower_to_qubit_en=False, dimension=3):
    def apply_rule(info, rule, rule_num):
        node_info_list = rule
        if len(node_info_list) > random_value:
          node_info = node_info_list[random_value]
        else:
          node_info = node_info_list[0]
        getattr(info, f'apply_rule_{rule_num}', lambda **kwargs: print("No function mapped for the given variable value"))(target=node_info)
        choices.append("applying " + str(rule_num) + str(node_info))
        assert info.slice_edge_collection.gate_slice_collection == info
        rewrite_applied = str(rule_num)
        return rewrite_applied 
    
    decision_random_values = [seed]
    random.seed(seed)
    unitary1 = None
    if info == "loaded":
      info = RewriteEngine.read_in(name + ".json", dimension)
    else:
      info, unitary1 = compile_basic(name=name, manual_lifting=manual_lifting, json_data=json_data, if_verify=if_verify, dimension=dimension)
    compilation_true = True
    apply_two = False
    cnt = 0
    choices = []
    disabled_counter = 0
    rewrite_applied = ""
    info.unroll_gate_common_phases()
    if if_return_basic:
      circ1, _ = info.build_circuit()
      return circ1
    while compilation_true:
      random_value_orig = random.random()
      random_value = int(random_value_orig * 2)
      decision_random_values.append(random_value)
      rule1 = info.get_all_rule_1_sites()
      rule2 = info.get_all_rule_2_sites()
      rule3 = info.get_all_rule_3_sites()
      compilation_true = rule1 is not None or\
        rule2 is not None or\
          rule3 is not None
      rand2 = random.random()
      rand2 = int(rand2 * 2)
      if apply_two:
        disabled_counter += 1
        if disabled_counter > 2:
          apply_two = False
      if rewrite_applied is None:
         cnt += 1
         if cnt > 5:
            break
      rewrite_applied = None
      if not apply_two and rule2 is not None:
        rewrite_applied = apply_rule(info, rule2, 2)
        apply_two = True
        continue
      if apply_two and rule3 is not None:
        rewrite_applied = apply_rule(info, rule3, 3)
        continue
      if rule1 is not None:
        rewrite_applied = apply_rule(info, rule1, 1)
    rule2 = info.get_all_rule_2_sites(check_control=False)
    if rule2 is not None:
      for r in rule2:
        info = apply_rule_2(info, r)
        info.slice_edge_collection.set_gate_slice_collection(info)
    resulting_merges = info.get_all_mergable_gates()
    for r in resulting_merges:
      info.apply_rule_gate_merge(r)
      info = delete_projected_gates(info, debug=False)
      info.slice_edge_collection.set_gate_slice_collection(info)
    info = gateslices_simplify(info)      
    circ2, qubits = info.build_circuit(if_final_circ=not "qubit" in compile_mode)
    if lower_to_qubit_en:
      info.unroll_multicontrolled_gates()
      info.convert_gateslices()
      verify_gateslices(info, if_verify, unitary1)
      result = info.build_circuit(dimension=2)
      print("check the dimension ", info.dimension)
      print("circuit width ", info.get_circuit_width())
      return result[0], result[1], []
    return circ2, qubits, choices

def count_circuits(circ_names):
  circ_directory, circ_name = circ_names.rsplit("/", 1)
  matching_files = []
  original_circ_name = circ_name
  circ_name = circ_name.replace("_graph", "")
  indexes = []
  for root, dirnames, filenames in os.walk(circ_directory):
    for circ_name in filenames:
      if "_graph.dot" in circ_name and circ_name.replace(".dot", "") != original_circ_name:
        indexes.append(re.findall(r'\d+', circ_name)[-1])
        matching_files.append(os.path.join(root, circ_name.replace(".dot", "")))
  indexes = [int(i) for i in indexes]
  matching_files = [x for _, x in sorted(zip(indexes, matching_files))]
  return matching_files

def compile_compare(name, manual_lifting, lower_to_qubit_en=False, dimension=3):
  info = compile_basic(name=name, manual_lifting=manual_lifting, if_base=True, if_verify=False, dimension=dimension)
  if lower_to_qubit_en:
    info.unroll_multicontrolled_gates()
    info.convert_gateslices()
    result = info.build_circuit(dimension=2)
    return result[0], result[1]
  circ2, qubits = info.build_circuit(if_final_circ=not "qubit" in compile_mode)
  return circ2, qubits

def compile(args, data):
  if "optimized" in args.compile_mode:
    # Preprocessing for nam circuits
    start_total = time.time()
    depth_arr = []
    seed_arr = []
    choices_arr = []
    random.seed(0)
    start_time = time.time()
    for i in range(int(args.runs)):
      seed = int(1000000 * random.random())
      seed_arr.append(seed)
      if not args.loaded:
        circ_toffoli_greedy, qubits, choices = compile_random("loaded" if i > 0 else None, name=args.test_name, manual_lifting=args.use_lifted,
                                                                seed=seed, json_data=data, itr=i, if_verify=args.verify_en,
                                                                compile_mode=args.compile_mode, dimension=args.dimension)
      else:
        circ_toffoli_greedy, qubits, choices = compile_random("loaded", name=args.test_name, manual_lifting=args.use_lifted,
                                                                seed=seed, json_data=data, itr=i, if_verify=args.verify_en,
                                                                compile_mode=args.compile_mode, dimension=args.dimension) 
      estimate_value = depth_estimater(circ_toffoli_greedy, qubits, mode=args.compile_mode)
      print("ESTIMATED DEPTH VALUE ", estimate_value)
      depth_arr.append(estimate_value)
      choices_arr.append(choices)
      data["random_compilation_itr_" + str(i)] = {
        "estimated_depth": depth_arr,
        "seeds": seed,
        "rewrites_applied": choices
      }
      print("After ", i, " runs")
      print("depth ", depth_arr)
      print("seeds ", seed_arr)
      print("choices ", choices)
    if len(depth_arr) > 0:
      min_seed = seed_arr[depth_arr.index(min(depth_arr))]
    else:
        min_seed = int(100000 * random.random())
    circ_toffoli_optimized, qubits, choices = compile_random("loaded", args.test_name, args.use_lifted, seed=min_seed, json_data=data,
                                                              itr="final", if_verify=args.verify_en, compile_mode=args.compile_mode,
                                                              lower_to_qubit_en="qubit" in args.compile_mode, dimension=args.dimension)
    data["random_compilation_itr_final"] = {
        "estimated_depth": min(depth_arr),
        "seeds": min_seed,
        "rewrites_applied": choices
      }
    print("POST compile circ ", args.run_sim)
    print(choices)
    print(circ_toffoli_optimized)
    end_time = time.time()
    assert args.run_sim == False
    post_compile_depth, two_q_gc, one_q_gc, simulation_results_dir = post_rewrite_compile(circ_toffoli_optimized, qubits,
                                                                                          run_sim=args.run_sim and "nam" not in args.test_name,
                                                                                          if_verify=args.verify_en, compile_mode=args.compile_mode)
    lowering_time = time.time()
    print("Final Depth Check ", post_compile_depth)
    print("OPTIMIZING ", end_time - start_time, "lower_time ", "full time ", lowering_time - end_time, end_time-start_total)
    end_total = time.time()
    data["random_compilation_itr_final"] = {
        "estimated_depth": min(depth_arr),
        "seeds": min_seed,
        "rewrites_applied": choices,
        "optimizing_time": end_time - start_time,
        "lowering_time": lowering_time - end_time,
        "depth": post_compile_depth,
        "2_q_gates": two_q_gc,
        "1_q_gates": one_q_gc,
        "fidelity": simulation_results_dir,
        "total_time": end_total - start_total
      }
    if args.qubit_sim and "lifted" in args.test_name and args.use_lifted:
      _, _, res_arr = run_qubit_decomposition(args.test_name.replace("_lifted", "_"), True)
    elif args.qubit_sim and "graph" in args.test_name and not args.use_lifted:
      _, _, res_arr = run_qubit_decomposition(args.test_name.replace("_graph", "_"), True)
    elif args.qubit_sim:
      raise Exception("The Qubit Simulation Has Not Run")
  elif "bench" in args.compile_mode:
      start_time = time.time()
      lowering_time = time.time()
      circ, qubits = compile_compare(args.test_name, args.use_lifted,
                                     lower_to_qubit_en="qubit" in args.compile_mode, dimension=2)
      print(qubits)
      post_compile_depth, two_q_gc, one_q_gc, simulation_results_dir = post_rewrite_compile(circ, qubits,
                                                                                            run_sim=args.run_sim,
                                                                                            compile_mode=args.compile_mode)
      data["manual_rewrites"] = {
          "estimated_depth": "",
          "seeds": "",
          "rewrites_applied": None,
          "optimizing_time": 0,
          "lowering_time": lowering_time - start_time,
          "depth": post_compile_depth,
          "2_q_gates": two_q_gc,
          "1_q_gates": one_q_gc,
          "fidelity": simulation_results_dir
      }
      print("Final Depth Check ", post_compile_depth)
  else:
    if args.test_name != "toffoli_5_lifted":
      raise Exception("This test is not supported yet for manual optimization")
    else:
      # TODO: Patch the required test here
      raise Exception("This test is supported, needs to be patched")  
  with open(args.json_name, 'w') as output_file:
    json.dump(data, output_file, indent=4)
  
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--use-lifted", default=False, action='store_true')
  parser.add_argument("--test-name", default="toffoli_7_lifted")
  parser.add_argument("--json-name", default="toffoli_7_lifted")
  parser.add_argument("--compile-mode", default="optimized")
  parser.add_argument("--runs", default=10)
  parser.add_argument("--seed", default=0)
  parser.add_argument("--dimension", default=2)
  parser.add_argument("--dont-final-replace", default=True, action='store_false')
  parser.add_argument("--qubit-sim", default=False, action='store_true')
  parser.add_argument("--run-sim", default=False, action='store_true')
  parser.add_argument("--verify-en", default=False, action='store_true')
  parser.add_argument("--loaded", default=False, action="store_true")
  parser.add_argument("--no-rccx-lowering", default=False, action="store_true")
  args = parser.parse_args()
  if os.path.exists(args.json_name):
    with open(args.json_name, 'r') as json_file:
      data = {} # json.load(json_file)
  else:
    data = {}  
  compile(args, data)