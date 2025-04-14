# import re
# from cirq.contrib.qasm_import import circuit_from_qasm
# from qiskit import QuantumCircuit
# import cirq
# import collections
# from benchmarks.decomposition_tests import get_circuit_unitaries, attain_symmetry
# import numpy as np
# from src.create_qutrit_circ import get_and_print_dag, save_mapping
# from rewriter.simple_circuit_rewrites import loop_through_rewrites
# from src.cirq_depth_counter import depth_counter
# from src.toffoli import *
# from collections import defaultdict
# from src.create_qutrit_circ import get_and_print_dag, create_qutrit_from_qubit
# from replacement_lib.replace_qubit_circ import replace_qubit_circuit

# class circuit_info:
#     def __init__(self, circ_name):
#         self.circuit = cirq.Circuit()
#         self.circuit_active_qubits = []
#         self.circuit_all_qubits = []
#         self.qubit_map = {}
#         if "check" in circ_name:
#             self.threshold = 5 # 6
#         else:
#             self.threshold = 4
#         self.control_for_target = {}

#     def add_node(self, op):
#         qubits = [x.name for x in op.qubits]
#         if not set(qubits) <= set(self.circuit_all_qubits):
#             self.circuit_all_qubits = list(set([x.name for x in op.qubits]).union(set(self.circuit_all_qubits)))
#         if not set(qubits) <= set(self.circuit_active_qubits):
#             self.circuit_active_qubits = list(set([x.name for x in op.qubits]).union(set(self.circuit_active_qubits)))
#         self.circuit.append(op)
#         if len(op.qubits) > 1 and op.qubits[-1].name in self.control_for_target.keys():
#             controls_arr = [x.name for x in op.qubits[:-1]]
#             running_control = set(controls_arr).union(set(self.control_for_target[op.qubits[-1].name]))
#             for c in controls_arr:
#                 if c in self.control_for_target.keys():
#                     running_control = set(running_control).union(set(self.control_for_target[c]))
#             # Update for target
#             self.control_for_target[op.qubits[-1].name] = list(running_control)
#             # update for controls
#             for c in self.control_for_target:
#                 if c != op.qubits[-1].name:
#                     if op.qubits[-1].name in self.control_for_target[c]:
#                         self.control_for_target[c] = list(set(self.control_for_target[c]).union(controls_arr))
#         elif len(op.qubits) > 1:
#             controls_arr = [x.name for x in op.qubits[:-1]]
#             running_control = controls_arr
#             for c in controls_arr:
#                 if c in self.control_for_target.keys():
#                     running_control = set(running_control).union(set(self.control_for_target[c]))
#             # Update for target
#             self.control_for_target[op.qubits[-1].name] = list(running_control)
#             for c in self.control_for_target:
#                 if c != op.qubits[-1].name:
#                     if op.qubits[-1].name in self.control_for_target[c]:
#                         self.control_for_target[c] = list(set(self.control_for_target[c]).union(controls_arr))


#     def check_if_parent_not_in(self, op):
#         for c in op.qubits[:-1]:
#             # print("CHECK FOR ", c)
#             if c.name in self.control_for_target.keys() and op.qubits[-1].name in self.control_for_target[c.name]:
#                 return False
#         return True

#     def check_active_node(self, op):
#         if not self.check_if_parent_not_in(op):
#             return False
#         qubits = [x.name for x in op.qubits]
#         if set(qubits) <= set(self.circuit_active_qubits):
#             return True
#         else:
#             return False
    
#     def check_non_active_node(self, op):
#         qubits = [x.name for x in op.qubits]
#         if len(self.circuit_all_qubits) < self.threshold and len(list(set(self.circuit_all_qubits).union(set(qubits)))) < self.threshold and self.check_if_parent_not_in(op):
#             for x in qubits:
#                 if x not in self.circuit_active_qubits and x in self.circuit_all_qubits:
#                     print("x and ", x, self.circuit_active_qubits)
#                     return False
#             return True
#         return False
        
#     def get_overlap_score(self, op):
#         qubits = [x.name for x in op.qubits]
#         return len(list(set(self.circuit_active_qubits).intersection(set(qubits))))
    
#     def remove_from_active_list(self, op, target_circuit=cirq.Circuit()):
#         return_val = False
#         for q in op.qubits:
#             if q.name in self.circuit_active_qubits:
#                 assert q.name in self.circuit_all_qubits
#                 self.circuit_active_qubits.remove(q.name)
#                 return_val = True
#         return return_val
    
#     def __str__(self):
#         print(self.circuit)
#         print(self.circuit_active_qubits)
#         return ""
    
#     def all_not_greater(self, cost):
#         for t in self.circuit_all_qubits:
#             if t in self.control_for_target:
#                 for elems in self.control_for_target[t]:
#                     if elems in cost and cost[t] <= cost[elems]:
#                         return False
#         return True        

#     def get_ordering(self):
#         cost = {}
#         for t in self.circuit_all_qubits:
#             cost[t] = 1
#         i = 0
#         while i < 10 and not self.all_not_greater(cost):
#             for t in self.control_for_target:
#                 for elems in self.control_for_target[t]:
#                     cost[t] = cost[t] + cost[elems] + 1
#             i += 1
#         return cost
    

# class circuit_cut:
#     def __init__(self, circuit, circuit_name):
#         self.qubits = circuit.all_qubits()
#         self.qubit_map = {}
#         for q in circuit.all_qubits():
#             self.qubit_map[q.name] = q.name.split("_")[-1]
#         self.circuit_info_list = []
#         self.circuit = circuit
#         self.order = {}
#         self.name = circuit_name

#     def cut_circuit(self):
#         for m_idx, moment in enumerate(self.circuit):
#             for op in moment:
#                 # print("op is here")
#                 # print(op)
#                 self.add_node_to_circuit(op)
#                 # for c in self.circuit_info_list:
#                 #     print("\t", c.control_for_target)
#                 # # print("   ===   ")

#     def add_node_to_circuit(self, op):
#         if self.add_to_active_circ(op):
#             # print("add to active circ")
#             return
#         if self.find_circuit_to_add(op):
#             # print("finding circuit")
#             return
#         self.create_new_circuit(op)
    
#     def add_to_active_circ(self, op):
#         for i, c in enumerate(self.circuit_info_list):
#             if c.check_active_node(op):
#                 # print("add to idx i ", i)
#                 c.add_node(op)
#                 return True
            
#     def has_max_order(self, index, op):
#         result = []
#         result_order = []
#         for i, circ in enumerate(self.circuit_info_list):
#             if i != index:
#                 for q in op.qubits:
#                     q = q.name
#                     if q in circ.circuit_active_qubits:
#                         result.append(circ)
#                         result_order.append(self.order[i])
#                         break
#         if len(result_order) == 0:
#             return True
#         if self.order[index] >= min(result_order):
#             return True
#         else:
#             return False

#     def find_circuit_to_add(self, op):
#         overlap = [c.get_overlap_score(op) for c in self.circuit_info_list]
#         if len(overlap) == 0:
#             return False
#         max_val = max(overlap)
#         overlap_copy = sorted(overlap, reverse=True)
#         indices = [index for index, _ in sorted(enumerate(overlap), key=lambda x: x[1], reverse=True)]
#         for i, max_val in enumerate(overlap_copy):
#             if max_val != 0:
#                 max_idx = indices[i]
#                 if self.circuit_info_list[max_idx].check_non_active_node(op) and self.has_max_order(max_idx, op):
#                     self.circuit_info_list[max_idx].add_node(op)
#                     max_order = -1
#                     for i, c in enumerate(self.circuit_info_list):
#                         if i != max_idx:
#                             if c.remove_from_active_list(op):
#                                 max_order = max(self.order[i], max_order)
#                     self.order[max_idx] = max(max_order + 1, self.order[max_idx])
#                     return True
#             else:
#                 return False
#                 print("not found")
#         return False
    
#     def find_order(self, max_order):
#         for i, c in enumerate(self.circuit_info_list):
#             if max_order == self.order[i]:
#                 max_order += 1
#         return max_order

#     def create_new_circuit(self, op):
#         self.circuit_info_list.append(circuit_info(self.name))
#         self.circuit_info_list[-1].add_node(op)
#         max_idx = len(self.circuit_info_list) - 1         
#         max_order = -1
#         self.order[max_idx] = max_order
#         for i, c in enumerate(self.circuit_info_list):
#             if i < len(self.circuit_info_list) - 1:
#                 if c.remove_from_active_list(op):
#                     max_order = max(self.order[i], max_order)
#         self.order[max_idx] = max_order + 1

#     def greedy_pruning(self):
#         order_list = list(self.order.values())
#         combined_lists = list(zip(order_list, self.circuit_info_list))
#         # Sort the combined list based on the first element of each pair
#         sorted_combined_lists = sorted(combined_lists, key=lambda x: x[0])
#         # Extract the sorted lists back into separate lists
#         order_list, circuit_info_list = zip(*sorted_combined_lists)
#         output_circuit_info = circuit_info(self.name)
#         full_circuit = cirq.Circuit()
#         final_circuits = []
#         major_mappings = []

#         for i, c in enumerate(circuit_info_list):
#             ifmerge = False
#             for _, moment in enumerate(c.circuit):
#                 for op in moment:
#                     if not output_circuit_info.check_non_active_node(op):
#                         ifmerge = False
#                     # Need to merge mappings

#             full_circuit.append(c.circuit)
#             if ifmerge:
#                 output_circuit_info.circuit.append(c.circuit)
#                 output_circuit_info.circuit_all_qubits = list(set(output_circuit_info.circuit_all_qubits).union(set(c.circuit_all_qubits)))
#                 if set(output_circuit_info.control_for_target.keys()) & set(c.control_for_target.keys()):
#                     assert False
#                 print(output_circuit_info.control_for_target, " - ", c.control_for_target)
#                 output_circuit_info.control_for_target = {**output_circuit_info.control_for_target, **c.control_for_target}
#                 print("res ios ", output_circuit_info.control_for_target, " - ", c.control_for_target)

#             else:
#                 final_circuits.append(c.circuit)
#                 major_mappings.append(c.get_ordering())
#                 # print(final_circuits[-1])
#                 # print(major_mappings[-1])
#                 # print("-------------------------------------")
#                 output_circuit_info = circuit_info(self.name)
#                 output_circuit_info.circuit.append(c.circuit)
#                 output_circuit_info.control_for_target = c.control_for_target.copy()
#                 output_circuit_info.circuit_all_qubits = list(set(output_circuit_info.circuit_all_qubits).union(set(c.circuit_all_qubits)))
#                 # print("res is ", output_circuit_info.control_for_target, " - ", c.control_for_target)
#                 # print("   ", c.get_ordering())    
        
#         # print("PRUNING RES")
#         # print(len(self.circuit_info_list))
#         # print(len(final_circuits))
#         # print(len(major_mappings))
#         return final_circuits, major_mappings

#     def verify(self, final_circuits):
#         return True
    
#     def print_circuits(self):
#         for i, c in enumerate(self.circuit_info_list):
#             print("position ", self.order[i])
#             print(c)
    
#     def get_all_orderings(self):
#         result = []
#         for i in self.circuit_info_list:
#             result.append(i.get_ordering())
#         return result
    
# def cut_circuit(circuit, name):
#     cutter = circuit_cut(circuit, name)
#     cutter.cut_circuit()
#     for c in cutter.circuit_info_list:
#        print(c)
#        print(c.control_for_target)
#        print("__________")
#     # assert False
#     orderings = cutter.get_all_orderings()
#     circuit_list, major_mappings = cutter.greedy_pruning()
#     print(major_mappings)
#     return circuit_list, major_mappings

# def circuit_mapper(circuits, ordering):
#     def index_of(mapping, q):
#         return mapping.index(q.name)
    
#     def cost_create(circuit, mapping):
#         cost = [0] * len(mapping)
#         ops  = []
#         qbts = []
#         for m in circuit:
#             for op in m:
#                 if len(op.qubits) > 1:
#                     cost[index_of(mapping, op.qubits[-1])] += 9*(len(op.qubits) - 1) - 8
#                 print("check ", op.qubits, cost)
#                 ops.append(op)
#                 qbts.append(len(op.qubits))
#         paired_elements = list(zip(qbts, ops))
#         sorted_pairs = sorted(paired_elements, key=lambda x: x[0])
#         ops = [x[1] for x in sorted_pairs]
#         print("curr cost is ", cost)
#         decay = len(circuit)
#         for m in circuit:
#             for op in m:
#                 if len(op.qubits) > 1:
#                     max_val = -1
#                     for x in op.qubits[:-1]:
#                         max_val = max(max_val, cost[index_of(mapping, x)])
#                     cost[index_of(mapping, op.qubits[-1])] = max(max_val, cost[index_of(mapping, op.qubits[-1])]) + decay
#             decay -= 1
#         return cost
    
#     res = []
#     mappers = []
#     # print("order ", len(circuits), len(ordering))
#     for i, c in enumerate(circuits):
#         mapping_orig = [q.name for q in c.all_qubits()]
#         mapping2 = [int(re.findall(r'\d+', x)[-1]) for x in mapping_orig]
#         paired_elements = list(zip(mapping2, mapping_orig))
#         sorted_pairs = sorted(paired_elements, key=lambda x: x[0])
#         mapping = [x[1] for x in sorted_pairs]
#         cost_arr = []
#         # USe gate ordering
#         for q in mapping:
#             if q in ordering[i]:
#                 cost_arr.append(ordering[i][q])
#             else:
#                 cost_arr.append(0)
#         # print(c)
#         # print(cost_arr)
#         # print(mapping)
#         # print("@@@@@@@@@@")
#         # DOnt use gate orderimg
#         paired_elements = list(zip(cost_arr, mapping))
#         sorted_pairs = sorted(paired_elements, key=lambda x: x[0])
#         print("pairs - ", sorted_pairs)
#         # Extract the first array elements from the sorted pairs
#         # print(sorted_pairs)
#         mapping = [x[1] for x in sorted_pairs]
#         # print("   ", mapping)
#         new_c = cirq.Circuit()
#         q = cirq.LineQid.range(len(mapping), dimension=2)
#         for m in c:
#             for op in m:
#                 if op.qubits == 1:
#                     new_c.append(cirq.MatrixGate(cirq.unitary(op)).on(q[index_of(mapping, op.qubits[0])]))
#                 else:
#                     matrix = cirq.unitary(op)
#                     matrix = matrix[-2:, -2:]
#                     control = []
#                     for q_ in op.qubits[:-1]:
#                         control.append(q[index_of(mapping, q_)])
#                     new_c.append(cirq.MatrixGate(matrix).on(q[index_of(mapping, op.qubits[-1])]).controlled_by(* control))
#         # print(new_c)
#         # print(mapping)
#         # print("@@@@@")
#         res.append(new_c)
#         mappers.append(mapping)
#     return res, mappers

# def merge_sequential_operations(c):
#     dict_arr = []
#     for m in c:
#         temp_dict = {}
#         for op in m:
#             temp_dict[op.qubits] = op
#         dict_arr.append(temp_dict)

#     for m in range(len(dict_arr)):
#         mark = []
#         for q in dict_arr[m].keys():
#             if m < len(dict_arr) - 1 and q in dict_arr[m+1].keys():
#                 if np.allclose(cirq.unitary(dict_arr[m][q]) @ cirq.unitary(dict_arr[m+1][q]),
#                                np.eye(cirq.unitary(dict_arr[m][q]).shape[0])):          
#                     del dict_arr[m+1][q]
#                     mark.append(q)
#         for q in mark:
#             del dict_arr[m][q]
#     circuit = cirq.Circuit()
#     for m in range(len(dict_arr)):
#         for q in dict_arr[m].keys():
#             circuit.append(dict_arr[m][q])
#     return circuit
            
# def improve_circuits(circuits):
#     result = []
#     for c in circuits:
#         result.append(merge_sequential_operations(c))
#     return result

# def get_qutrit_unitary_circ(circuits, mapping, attain_symmetry_true=False, dimension=3):
#     result = []
#     qs = []
#     for idx, c in enumerate(circuits):
#         size = len(mapping[idx])
#         if size != len(c.all_qubits()):
#             print(c)
#             print(mapping[idx])
#             print(c.all_qubits())
#         assert size == len(c.all_qubits())

#         # Check attain symmetry
#         print("Attain symmertry is enabled: ", attain_symmetry_true)
#         if attain_symmetry_true and False:
#             matrix = attain_symmetry(matrix.copy())
#         circ2 =  cirq.Circuit()
#         q = cirq.LineQid.range(size, dimension=dimension)
#         matrix = get_circuit_unitaries(c, dim=2)
#         circ2.append(cirq.MatrixGate(matrix, qid_shape=[3]*(size)).on(* q))
#         result.append(circ2)
#         # result.append(create_qutrit_from_qubit(c, len(q))[0])
#         qs.append(q)
#     return result, qs

# def assembled_circuit(circuits, circuit_name, dimension=3, not_saved=False, mappings=None):
#     path = circuit_name
#     final_circuit = cirq.Circuit()
#     max_mapping = 0
#     print("assembling circuit")
#     if mappings is None:
#         mappings = []
#         for i, c in enumerate(circuits):
#             file = open(path.replace("_graph", "") + "_" + str(i) + "_mapping.txt", "r")
#             mapping = [int(x.split("_")[-1]) for x in file.read().split(',')]
#             mappings.append(mapping)
#             print(c)
#             print(mapping)
#             max_mapping = max(max(mapping), max_mapping)
#     else:
#         for i, mapping in enumerate(mappings):
#             mappings[i] = [int(re.findall(r'\d+', x)[-1]) for x in mapping]
#             print(circuits[i])
#             print(mapping)
#             max_mapping = max(max(mappings[i]), max_mapping)    
#     # assert False
#     qubits = cirq.LineQid.range(max_mapping+1, dimension=dimension)
#     qmap = {}
#     for q in qubits:
#         qmap[q.x] = q
#     for i, c in enumerate(circuits):
#         for _, m in enumerate(c):
#             for op in m:
#                 q = op.qubits
#                 if len(q) == 1:
#                     final_circuit.append(cirq.MatrixGate(cirq.unitary(op), qid_shape=[dimension]).on(qmap[mappings[i][op.qubits[-1].x]]))
#                 else:
#                     matrix = cirq.unitary(op)
#                     matrix = matrix[(matrix.shape[0] - dimension) // (dimension-1) : min(matrix.shape[0], (matrix.shape[0] + dimension) // (dimension-1)),
#                                     (matrix.shape[0] - dimension) // (dimension-1) : min(matrix.shape[0], (matrix.shape[0] + dimension) // (dimension-1))]
#                     controls = [qmap[mappings[i][x.x]] for x in op.qubits[:-1]]
#                     if np.allclose(matrix, cirq.unitary(QutritX01Gate())[:matrix.shape[0], :matrix.shape[0]]):  
#                         if dimension == 2:
#                             final_circuit.append(cirq.X(qmap[mappings[i][op.qubits[-1].x]]).controlled_by(* controls))
#                         else:
#                             final_circuit.append(QutritX01Gate().on(qmap[mappings[i][op.qubits[-1].x]]).controlled_by(* controls))    
#                     else:
#                         final_circuit.append(cirq.MatrixGate(matrix, qid_shape=[dimension]).on(qmap[mappings[i][op.qubits[-1].x]]).controlled_by(* controls))
#     if not not_saved:
#         folder_name, circuit = circuit_name.rsplit("/", 1)
#         get_and_print_dag(final_circuit, qubits, circuit_name.replace("_graph", ""))
#     return final_circuit

# def verify_circuits(circ1, circ2):
#     qudits = []
#     len1 = 0
#     len2 = 0
#     for i, m in enumerate(circ1):
#         for op in m:
#             if len(op.qubits) > 1:
#                 len1 += 1     
#     for i, m in enumerate(circ2):
#         for op in m:
#             if len(op.qubits) > 1:
#                 len2 += 1
#     assert len1 == len2

# def get_nam_circuits(circuit_name, return_base=False, return_qubit=False):
#     # nam-benchmarks
#     file = open("/nobackupkiwi/rsharma3/QC_work/sara_work/qutrits/cirq_direct_clone/dare_check/quartz/circuit/nam-benchmarks/" + circuit_name.replace("nam_", "").replace("_check", "") + ".qasm", "r")
#     print("file received")
#     lines = file.read()
#     lines = lines.splitlines()
#     print("lines fetched")
#     # Process each line (for example, adding a prefix)
#     processed_lines = []
#     for line in lines:
#         line = line.replace("qubits", "q")
#         if "ccz" in line:
#             target = last_integer = re.findall(r'\b\d+\b|\[\d+\]', line)[-1]
#             processed_lines.append("h q" + target + ";")
#         processed_lines.append(line.replace("ccz", "ccx"))
#         if "ccz" in line:
#             # tagret = last_integer = re.findall(r'\b\d+\b|\[\d+\]', line)[-1]
#             processed_lines.append("h q" + target + ";")
#     # Merge the processed lines back into the same string using join() method
#     merged_string = '\n'.join(processed_lines)
#     circuit = circuit_from_qasm(merged_string)
#     original_circ = circuit.copy()
#     print("original")
#     print(original_circ)
#     print("merging")
#     circuit = merge_sequential_operations(circuit)
#     orig2 = circuit.copy()
    
#     if return_base:
#         return circuit, circuit.all_qubits()
#     if return_qubit:
#         dummy_circ2 = cirq.Circuit()
#         dummy_circ2.append(replace_qubit_circuit(original_circ, original_circ.all_qubits()))
#         return dummy_circ2, dummy_circ2.all_qubits()

#     # Convert ccx to toffoli impl to get depth
#     depth_all, _, ttl_2, ttl_1 = depth_counter(circuit, cirq.LineQid.range(len(circuit.all_qubits()), dimension=2))
#     print("BENCHMARKING ", depth_all, _, ttl_2, ttl_1)
#     circuit = loop_through_rewrites(circuit, qutrits=None, checking_blocked=False, d=2)
#     print("cut")
#     circuit_list, ordering = cut_circuit(circuit, circuit_name)
#     print("map")
#     circuit_list, mapping = circuit_mapper(circuit_list, ordering)
#     print("ass")
#     circ_assembled = assembled_circuit(circuit_list, None, dimension=2, not_saved=True, mappings=mapping)
#     print("assembled circ ")
#     print(original_circ)
#     print(orig2)
#     print(circ_assembled)
#     verify_circuits(circ_assembled, original_circ)    
#     print("print orderings are ", len(ordering))
#     print("print circuis are ", len(circuit_list))

#     for i, c in enumerate(circuit_list):
#         print("Partial circuit is - ", mapping[i])
#         print("\t\t", ordering[i])
#         print(c)

#     if "check" in circuit_name:
#         circuit_list, qs = get_qutrit_unitary_circ(circuit_list, mapping, attain_symmetry_true=True)
#     else:
#         circuit_list, qs = get_qutrit_unitary_circ(circuit_list, mapping, attain_symmetry_true=False)    
#     # assert False
#     return circuit_list, (qs, mapping)
