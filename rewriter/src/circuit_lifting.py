import cirq
import scipy
import pydot
import argparse
import ast
import networkx as nx
import numpy as np
import scipy.linalg
import numpy as np
import collections
import re
from src.qutrit_decompositions import *
from src.toffoli import *
from rewriter.src.rustxgraph_circuit import if_mat_preserves_binary
from rewriter.src.node import Node
from rewriter.src.node_directory import NodeDirectory
# from rustxgraph_circuit import GraphInfo  

def get_qubits(node_entry, debug=False):
    node_entry = node_entry.get_label()
    pattern = r'q\(.+?\)'
    matches = re.findall(pattern, node_entry)
    result_qubits = []
    for match in matches:
        result_qubits.append(int(match.replace("q(", "").replace(")", "")))
    new_arr = result_qubits[:-1].copy()
    new_arr.sort()
    new_arr.append(result_qubits[-1])
    return result_qubits

def get_matrix(node_entry):
    node_entry = node_entry.get_label()  # ["label"]
    node_entry = node_entry.replace("\n", "")
    matrix_entry = r'matrix.+\[\[.+\]\]'
    matches = re.findall(matrix_entry, node_entry)
    if len(matches) == 0:
        if "QutritX01Gate" in node_entry:
            return cirq.unitary(QutritX01Gate())
        if "QutritX12Gate" in node_entry:
            return cirq.unitary(QutritX12Gate())
        if "QutritX02Gate" in node_entry:
            return cirq.unitary(QutritX02Gate())
        if "QutritPlusGate" in node_entry:
            return cirq.unitary(QutritPlusGate())
        if "QutritCC1PlusGate" in node_entry:
            return cirq.unitary(QutritPlusGate())
        if "QutritCC1PlusGate" in node_entry:
            return cirq.unitary(QutritPlusGate())
        if "QutritMinusGate" in node_entry:
            return cirq.unitary(QutritMinusGate())
        if "QutritC1PlusGate" in node_entry:
          return cirq.unitary(QutritPlusGate())
        if "QutritC1MinusGate" in node_entry:
            return cirq.unitary(QutritMinusGate())
        if "QutritCC1MinusGate" in node_entry:
            return cirq.unitary(QutritMinusGate())
        if "QutritCCX101Gate" in node_entry:
            return cirq.unitary(QutritX01Gate())
        if "QutritCX101Gate" in node_entry:
            return cirq.unitary(QutritX01Gate())
        if "QutritCX112Gate" in node_entry:
            return cirq.unitary(QutritX12Gate())
        if "QutritCX102Gate" in node_entry:
            return cirq.unitary(QutritX02Gate())
        if "QutritCCX102Gate" in node_entry:
            return cirq.unitary(QutritX02Gate())
        if "QutritCCX112Gate" in node_entry:
            return cirq.unitary(QutritX12Gate())
        if "QutritCCZ112Gate" in node_entry:
            return cirq.unitary(QutritZ12Gate())
        if "QutritCCZ101Gate" in node_entry:
            return cirq.unitary(QutritZ01Gate())
        else:
          pass          #  print(matches)
          #  if len(matches) != 0:
          # assert False
    # print(matches)
    if len(matches) > 0:
        match = matches[0]
        matrix_entry = r'\[\[.+\]\]'
        match = re.findall(matrix_entry, match)[0]
        if match == "None":
            return None
        matrix = np.array(ast.literal_eval(match), dtype=np.complex128)
        V, __, Wh = scipy.linalg.svd(matrix)
        U = np.matrix(V.dot(Wh))
        return U
    return None

def get_controls(node_entry):
    node_entry = node_entry.get_label()  # ["label"]
    matrix_entry = r'C+'
    entry = re.findall(matrix_entry, node_entry)
    if len(entry) > 0:
        entry = entry[0]
        return [1] * len(entry)
    return []

def lift_dag_circuit(graph_object, if_base=False):
    for nodes in graph_object.node_info:
        graph_object.node_info[nodes].set_qubits(get_qubits(graph_object.node_info[nodes]))
        graph_object.node_info[nodes].set_matrix(get_matrix(graph_object.node_info[nodes]))
        graph_object.node_info[nodes].set_controls(get_controls(graph_object.node_info[nodes]))
        # It has no parents and is a qubit
        if len(list(graph_object.dag.predecessors(nodes))) == 0 and graph_object.node_info[nodes].get_matrix() is None:
            graph_object.node_info[nodes].set_root()
        elif len(list(graph_object.dag.predecessors(nodes))) == 0:
            # print(nodes, list(graph_object.dag.predecessors(nodes)), graph_object.node_info[nodes])
            raise Exception("This gate has no qubits feeding it") 
        if len(list(graph_object.dag.successors(nodes))) == 0 and graph_object.node_info[nodes].get_matrix() is None:
            graph_object.node_info[nodes].set_sink()
        elif len(list(graph_object.dag.successors(nodes))) == 0:
           raise Exception("This gate has no qubits outputting to")
    if not if_base:
      graph_object.create_mapping()
    return graph_object

def apply_control_rewrites(graph_object):
    qubits = graph_object.get_all_qubits()
    nxt_node = None
    for qubit_node in qubits:
        last_gate_marked = None
        qubit_num = graph_object.node_info[qubit_node].get_qubits()[0]
        nxt_node = None
        curr_node = None
        modify_gates_table = []
        modify_matrix_vals = {}
        # print("QUBIT NUM ", qubit_num)
        while True:
            if nxt_node is None:
                curr_node = qubit_node
                curr_node = graph_object.get_nxt_qubit_node(curr_node, qubit_num)
            else:
                curr_node = nxt_node
            if curr_node is None and nxt_node is None:
                break
            nxt_node = graph_object.get_nxt_qubit_node(curr_node, qubit_num)
            if graph_object.get_node_type(curr_node) == "qubit":
                # print("EXIT QUBIT")
                if graph_object.node_info[curr_node].is_sink():
                    # if qubit_num == 3:
                    #    print("insert for last gate ", last_gate_marked)
                    if last_gate_marked is not None and not np.allclose(np.eye(3), last_gate_marked):
                        graph_object.insert_new_gate(last_gate_marked, curr_node)
                    break
            elif graph_object.get_target(curr_node) == qubit_num and len(graph_object.get_control(curr_node)) == 0:
                graph_object.node_info[curr_node].set_status("delete")
                # if qubit_num ==3:
                #    print("st")
                #    print(last_gate_marked, curr_node)
                #    print(graph_object.node_info[curr_node])
                if last_gate_marked is None or np.allclose(last_gate_marked, np.eye(3)):
                    last_gate_marked = graph_object.node_info[curr_node].get_matrix()
                    modify_gates_table.append(curr_node)
                else:
                  # Decide if we want to use only this or others too
                  if True:
                    # if_mat_preserves_binary(graph_object.node_info[curr_node]["matrix"] @ last_gate_marked):
                    last_gate_marked = graph_object.node_info[curr_node].get_matrix() @ last_gate_marked
                    modify_gates_table.append(curr_node)
                    modify_matrix_vals[curr_node] = None

            elif graph_object.get_target(curr_node) != qubit_num:
              # Make last gate marked a pauli gate
              if last_gate_marked is None or if_mat_preserves_binary(last_gate_marked):
                 pass
              elif if_mat_preserves_binary(last_gate_marked @ cirq.unitary(QutritX02Gate())):
                pass
              elif if_mat_preserves_binary(last_gate_marked @ cirq.unitary(QutritX12Gate())):
                pass
              elif if_mat_preserves_binary(last_gate_marked @ cirq.unitary(QutritX01Gate())):
                 pass
              else:
                candidate = np.zeros((3, 1))
                candidate[1] = 1
                cadidate = np.reshape(last_gate_marked @ candidate, (3, 1))
                skip = False
                if not skip:    
                  gate_node = modify_gates_table[-1]
                  modify_matrix_vals[gate_node] = last_gate_marked
                  graph_object.node_info[gate_node].set_matrix(last_gate_marked)
                  last_gate_marked = np.eye(3)
                  graph_object.node_info[gate_node].set_status(None)

              graph_object.set_control_value(curr_node, qubit_num, last_gate_marked)

            elif graph_object.get_target(curr_node) == qubit_num and len(graph_object.get_control(curr_node)) > 0:
                if last_gate_marked is not None:
                    graph_object.node_info[curr_node].set_matrix(last_gate_marked.conjugate().T @
                                                                 graph_object.node_info[curr_node].get_matrix() @
                                                                 last_gate_marked)
                    if not np.allclose(graph_object.node_info[curr_node].get_matrix() @
                                       graph_object.node_info[curr_node].get_matrix().conjugate().T, np.eye(3)):
                      print("ERROR")
                      print(graph_object.node_info[curr_node].get_matrix())
                      assert False
    graph_object.delete_marked_gates()
    return graph_object

def mark_and_delete_gates(graph_obj):
    graph_obj.reset_status()
    status_dict = collections.defaultdict(dict)
    not_processed = False
    wire_status_condition = collections.defaultdict(dict)
    final_result_arr = []
    print("check start")
    # Working from Start
    for nodes in graph_obj.node_info.keys():
        status_dict[nodes]["node"] = graph_obj.node_info[nodes]
        status_dict[nodes]["status"] = graph_obj.node_info[nodes].is_root()
        # This node is not the node and entire graph has not yet been set to not_processed
        if not status_dict[nodes]["status"] and not not_processed:
            not_processed = True

    while not_processed:
      for nodes in graph_obj.node_info.keys():
        if not status_dict[nodes]["status"]:
          need_to_process = not status_dict[nodes]["status"]
          for parents in graph_obj.dag.predecessors(nodes):
            if status_dict[parents]["status"] and need_to_process:
              need_to_process = True
            else:
              need_to_process = False
          if need_to_process:
            status_dict[nodes]["status"] = True
            for parents in graph_obj.dag.predecessors(nodes):
              for k in graph_obj.dag.get_edge_data(parents, nodes).keys():
                for succs in graph_obj.dag.successors(nodes):
                  for k2 in graph_obj.dag.get_edge_data(nodes, succs).keys():
                    if int(graph_obj.edge_info[(parents, nodes, k)]["label"]) == \
                      int(graph_obj.edge_info[(nodes, succs, k2)]["label"]):
                      # Update the qubits and the corresponding result matrix
                      # If not the taget dont need to change the mapping
                      # I am assuming this node is always active what if not?
                      q = graph_obj.edge_info[(parents, nodes, k)]["qubits"]
                      if q == graph_obj.get_qubit(nodes)[-1] and \
                        len(graph_obj.get_qubit(nodes)) == 1:
                        pass
                      elif q == graph_obj.get_qubit(nodes)[-1] and \
                        len(graph_obj.get_qubit(nodes)) > 1:
                        # We have two options in this case
                        # Option 1 the result is always false
                        # result_matrix1 = mapping
                        temp_tuple = tuple([tuple(graph_obj.get_qubit(nodes)), tuple(graph_obj.get_all_control_value(nodes, if_ordered=True))])
                        del_arr = []
                        
                        # if temp_tuple in wire_status_condition[q]:
                        #   for qubits_and_ctrls in wire_status_condition[q]:
                        #     if wire_status_condition[q][qubits_and_ctrls][-1] < nodes and\
                        #       wire_status_condition[q][temp_tuple][-1] < wire_status_condition[q][qubits_and_ctrls][-1]:
                        #         if not if_mat_preserves_binary(graph_obj.get_matrix(nodes)):
                        #           qubits = list(qubits_and_ctrls[0])
                        #           ctrls = list(qubits_and_ctrls[1])
                        #           if len(set(qubits[:-1]).intersection(set(list(temp_tuple[0])[:-1]))) > 0 and set(temp_tuple[0]) <= set(qubits):
                        #             q_idx = [qubits.index(x) for x in temp_tuple[0][:-1]]
                        #             c_vals = [ctrls[x] for x in q_idx]
                        #             if c_vals == list(temp_tuple[1]):
                        #               mat_check = graph_obj.get_matrix(nodes) @ \
                        #                           graph_obj.get_matrix(wire_status_condition[q][qubits_and_ctrls][-1]) @ \
                        #                           graph_obj.get_matrix(wire_status_condition[q][temp_tuple][-1])
                        #               if if_mat_preserves_binary(mat_check):
                        #                 del_arr.append(qubits_and_ctrls)
                        # for arr in del_arr:
                        #   if len(wire_status_condition[q][arr]) > 1:
                        #     wire_status_condition[q][arr].pop(-1)
                        #   else:
                        #     del wire_status_condition[q][arr]

                        if temp_tuple not in wire_status_condition[q] and not \
                          (graph_obj.check_gate_type(nodes) == "01"): # or graph_obj.check_gate_type(nodes) == "10"):
                          wire_status_condition[q][temp_tuple] = [nodes]
                        else:
                          if temp_tuple in wire_status_condition[q]:
                            mat_check = graph_obj.get_matrix(nodes) @ graph_obj.get_matrix(wire_status_condition[q][temp_tuple][-1])
                            if if_mat_preserves_binary(mat_check):  #  or if_mat_preserves_binary(mat_check @ cirq.unitary(QutritX01Gate())):                  
                              if len(wire_status_condition[q][temp_tuple]) > 1:
                                wire_status_condition[q][temp_tuple].pop(-1)
                              else:
                                del wire_status_condition[q][temp_tuple]

                            elif not (graph_obj.check_gate_type(nodes) == "01"): # or graph_obj.check_gate_type(nodes) == "10"):
                              wire_status_condition[q][temp_tuple].append(nodes)
                      else:
                        q = graph_obj.edge_info[(parents, nodes, k)]["qubits"]
                        ctrls = graph_obj.get_all_control_value(nodes, if_ordered=True)
                        ctrl_q = graph_obj.get_control(nodes)
                        ctrl_stat = {ctrl_q[i]: ctrls[i] for i in range(len(ctrl_q))}
                        dont_change = False
                        
                        for elems in wire_status_condition[q]:
                          elem_qs = list(elems[0])
                          elem_cs = list(elems[1])
                        
                          if set(elem_qs) <= set(ctrl_q):
                            if elem_cs == [ctrl_stat[x] for x in elem_qs[:-1]]:
                              # print("got here ", graph_obj.node_info[wire_status_condition[q][elems]])
                              for i in wire_status_condition[q][elems]:
                                if graph_obj.check_gate_type(i) == "12" and ctrl_stat[elem_qs[-1]] > 0:
                                  ctrl_stat[elem_qs[-1]] = 3 - ctrl_stat[elem_qs[-1]]
                                elif graph_obj.check_gate_type(i) == "+":
                                    ctrl_stat[elem_qs[-1]] = (ctrl_stat[elem_qs[-1]] + 1) % 3
                                elif graph_obj.check_gate_type(i) == "-":
                                    ctrl_stat[elem_qs[-1]] = (ctrl_stat[elem_qs[-1]] - 1) % 3                                
                                elif graph_obj.check_gate_type(i) == "10" and ctrl_stat[elem_qs[-1]] < 2:
                                  ctrl_stat[elem_qs[-1]] = 1 - ctrl_stat[elem_qs[-1]]
                                elif graph_obj.check_gate_type(i) == "02" and ctrl_stat[elem_qs[-1]] != 1:
                                  ctrl_stat[elem_qs[-1]] = 2 - ctrl_stat[elem_qs[-1]]
                          elif elem_qs[-1] in ctrl_q:
                            # print("got here with check ", graph_obj.node_info[wire_status_condition[q][elems]])
                            for i in wire_status_condition[q][elems]:
                              if (graph_obj.check_gate_type(i) == "01" or
                                  graph_obj.check_gate_type(i) == "10") and ctrl_stat[elem_qs[-1]] == 2:
                                pass
                              elif graph_obj.check_gate_type(i) == "12" and ctrl_stat[elem_qs[-1]] == 0:
                                pass
                              elif graph_obj.check_gate_type(i) == "02" and ctrl_stat[elem_qs[-1]] == 1:
                                pass
                              else:
                                dont_change = True
                                break
                        # print(ctrl_q, dont_change)
                        if ctrl_stat[q] == 2 and not dont_change and nodes not in final_result_arr:
                          # print("APPENDING NODE ", nodes, graph_obj.node_info[nodes])
                          final_result_arr.append(nodes)
                        #   print("CAN DELETE node ", nodes, q, graph_obj.node_info[nodes], wire_status_condition[q])
      not_processed = False
      for nodes in status_dict.keys():
        if not status_dict[nodes]["status"]:
          not_processed = True

    status_dict = collections.defaultdict(dict)
    wire_status_condition = collections.defaultdict(dict)
    # Working from Emd
    for nodes in graph_obj.node_info.keys():
        status_dict[nodes]["node"] = graph_obj.node_info[nodes]
        status_dict[nodes]["status"] = status_dict[nodes]["node"].is_sink()
        # This node is not the node and entire graph has not yet been set to not_processed
        if not status_dict[nodes]["status"] and not not_processed:
            not_processed = True
    while not_processed:
      for nodes in graph_obj.node_info.keys():
        if not status_dict[nodes]["status"]:
          need_to_process = not status_dict[nodes]["status"]
          for succs in graph_obj.dag.successors(nodes):
            if status_dict[succs]["status"] and need_to_process:
              need_to_process = True
            else:
              need_to_process = False
          if need_to_process:
            status_dict[nodes]["status"] = True
            for succs in graph_obj.dag.successors(nodes):
              for k in graph_obj.dag.get_edge_data(nodes, succs).keys():
                for prevs in graph_obj.dag.predecessors(nodes):
                  for k2 in graph_obj.dag.get_edge_data(prevs, nodes).keys():
                    # print((prevs, nodes, k2), self.edge_info[(prevs, nodes, k2)], self.edge_info[(nodes, succs, k)])
                    if int(graph_obj.edge_info[(prevs, nodes, k2)]["label"]) == \
                      int(graph_obj.edge_info[(nodes, succs, k)]["label"]):
                      # Update the qubits and the corresponding result matrix
                      # If not the taget dont need to change the mapping
                      # I am assuming this node is always active what if not?
                      q = graph_obj.edge_info[(prevs, nodes, k2)]["qubits"]
                      
                      if q == graph_obj.get_qubit(nodes)[-1] and \
                        len(graph_obj.get_qubit(nodes)) == 1:
                        pass
                      elif q == graph_obj.get_qubit(nodes)[-1] and \
                        len(graph_obj.get_qubit(nodes)) > 1:
                        # We have two options in this case
                        # Option 1 the result is always false
                        temp_tuple = tuple([tuple(graph_obj.get_qubit(nodes)), tuple(graph_obj.get_all_control_value(nodes, if_ordered=True))])
                        del_arr = []
                        # if temp_tuple in wire_status_condition[q]:
                        #   for qubits_and_ctrls in wire_status_condition[q]:
                        #     if wire_status_condition[q][qubits_and_ctrls][-1] != wire_status_condition[q][temp_tuple][-1] and \
                        #       wire_status_condition[q][qubits_and_ctrls][-1] < wire_status_condition[q][temp_tuple][-1] and \
                        #         nodes < wire_status_condition[q][qubits_and_ctrls][-1]:
                        #         if not if_mat_preserves_binary(graph_obj.get_matrix(nodes)):
                        #           qubits = list(qubits_and_ctrls[0])
                        #           ctrls = list(qubits_and_ctrls[1])
                        #           if len(set(qubits[:-1]).intersection(set(list(temp_tuple[0])[:-1]))) > 0 and set(temp_tuple[0]) <= set(qubits):
                        #             q_idx = [qubits.index(x) for x in temp_tuple[0][:-1]]
                        #             c_vals = [ctrls[x] for x in q_idx]
                        #             if c_vals == list(temp_tuple[1]):
                        #               mat_check = graph_obj.get_matrix(wire_status_condition[q][temp_tuple][-1]) @ \
                        #                           graph_obj.get_matrix(wire_status_condition[q][qubits_and_ctrls][-1]) @ \
                        #                           graph_obj.get_matrix(nodes)
                        #               if if_mat_preserves_binary(mat_check): # or if_mat_preserves_binary(mat_check @ cirq.unitary(QutritX01Gate())):
                        #                   # print(temp_tuple, qubits_and_ctrls)
                        #                   del_arr.append(qubits_and_ctrls)
                        # for arr in del_arr:
                        #   if len(wire_status_condition[q][arr]) > 1:
                        #     wire_status_condition[q][arr].pop(-1)
                        #   else:
                        #     del wire_status_condition[q][arr]
                        if temp_tuple not in wire_status_condition[q] and not \
                          (graph_obj.check_gate_type(nodes) == "01"): #  or graph_obj.check_gate_type(nodes) == "10"):
                          wire_status_condition[q][temp_tuple] = [nodes]
                        else:
                          if temp_tuple in wire_status_condition[q]:
                            mat_check = graph_obj.get_matrix(wire_status_condition[q][temp_tuple][-1]) @ graph_obj.get_matrix(nodes)
                            if if_mat_preserves_binary(mat_check):
                              if len(wire_status_condition[q][temp_tuple]) > 1:
                                wire_status_condition[q][temp_tuple].pop(-1)
                              else:
                                del wire_status_condition[q][temp_tuple]

                            elif not (graph_obj.check_gate_type(nodes) == "01"): # or graph_obj.check_gate_type(nodes) == "10"):
                              wire_status_condition[q][temp_tuple].append(nodes)
                      else:
                        if nodes not in final_result_arr:
                          q = graph_obj.edge_info[(prevs, nodes, k2)]["qubits"]                          
                          ctrls = graph_obj.get_all_control_value(nodes, if_ordered=True)
                          ctrl_q = graph_obj.get_control(nodes)
                          ctrl_stat = {ctrl_q[i]: ctrls[i] for i in range(len(ctrl_q))}
                          dont_change = False
                          # if nodes == "55":
                          #    print(q, ctrls, ctrl_q)
                          #    print(wire_status_condition[q])

                          for elems in wire_status_condition[q]:
                            elem_qs = list(elems[0])
                            elem_cs = list(elems[1])
                            # if nodes == "55":
                            #    for a in wire_status_condition[q][elems]:
                            #     #  print(graph_obj.node_info[a])
                            if set(elem_qs) <= set(ctrl_q):
                              if elem_cs == [ctrl_stat[x] for x in elem_qs[:-1]]:
                                for i in wire_status_condition[q][elems]:
                                  if graph_obj.check_gate_type(i) == "12" and ctrl_stat[elem_qs[-1]] > 0:
                                    ctrl_stat[elem_qs[-1]] = 3 - ctrl_stat[elem_qs[-1]]
                                  elif graph_obj.check_gate_type(i) == "+":
                                    ctrl_stat[elem_qs[-1]] = (ctrl_stat[elem_qs[-1]] + 1) % 3
                                  elif graph_obj.check_gate_type(i) == "-":
                                    ctrl_stat[elem_qs[-1]] = (ctrl_stat[elem_qs[-1]] - 1) % 3                                
                                  elif graph_obj.check_gate_type(i) == "10" and ctrl_stat[elem_qs[-1]] < 2:
                                    ctrl_stat[elem_qs[-1]] = 1 - ctrl_stat[elem_qs[-1]]
                                  elif graph_obj.check_gate_type(i) == "02" and ctrl_stat[elem_qs[-1]] != 1:
                                    ctrl_stat[elem_qs[-1]] = 2 - ctrl_stat[elem_qs[-1]]
                            elif elem_qs[-1] in ctrl_q:
                              for i in wire_status_condition[q][elems]:
                                if (graph_obj.check_gate_type(i) == "01" or
                                    graph_obj.check_gate_type(i) == "10") and ctrl_stat[elem_qs[-1]] == 2:
                                  pass
                                elif graph_obj.check_gate_type(i) == "12" and ctrl_stat[elem_qs[-1]] == 0:
                                  pass
                                elif graph_obj.check_gate_type(i) == "02" and ctrl_stat[elem_qs[-1]] == 1:
                                  pass
                                else:
                                  dont_change = True
                          if ctrl_stat[q] == 2 and not dont_change:
                            # print("APPENDING AGAIN- ", q,  nodes, ctrl_stat, graph_obj.node_info[nodes])
                            final_result_arr.append(nodes)
                            # pass
      not_processed = False
      for nodes in status_dict.keys():
        if not status_dict[nodes]["status"]:
          not_processed = True
    for nodes in final_result_arr[:1]:
      graph_obj.node_info[nodes].set_status("delete")
    graph_obj.delete_marked_gates()
    return graph_obj
    
def delete_high_control_gates(graph_obj):
    # print("Start deleting")
    skip_to_end = False
    for nodes in graph_obj.node_info:
        if 2 in graph_obj.node_info[nodes].get_control_values():
            for control_idx in graph_obj.get_control(nodes):
                qubit = control_idx
                prev, nxt, _, s1, s2 = graph_obj.parent_and_child_node(nodes, q=qubit)
                if graph_obj.get_control_value(nodes, qubit) == 2 and (graph_obj.edge_info[(prev, nodes, s1)]["wire"] == 'b' or
                                                                       graph_obj.edge_info[(nodes, nxt, s2)]["wire"] == 'b'):
                    graph_obj.node_info[nodes].set_status("delete")
                    skip_to_end = True
                    break
        if skip_to_end:
           break
    graph_obj.delete_marked_gates()
    return graph_obj

def delete_high_phase_gates(graph_obj):
    for nodes in graph_obj.node_info:
        if graph_obj.node_info[nodes].get_matrix() is not None and \
            graph_obj.node_info[nodes].is_high_phase_gate():
            q = graph_obj.node_info[nodes].get_target()
            prev, nxt, _, s1, s2 = graph_obj.parent_and_child_node(nodes, q=q)
            if (graph_obj.edge_info[(prev, nodes, s1)]["wire"] == 'b' and graph_obj.edge_info[(nodes, nxt, s2)]["wire"] == 'b'):
                # print("node to be deleted")
                # print(graph_obj.node_info[nodes])
                graph_obj.node_info[nodes].set_status("delete")
    graph_obj.delete_marked_gates()
    return graph_obj

def apply_gate_flip_rewrites(graph_object):
    qubits = graph_object.get_all_qubits()
    nxt_node = None
    node_status = graph_object.node_info.copy()
    for qubit_node in qubits:
        qubit_num = graph_object.node_info[qubit_node].get_qubits()
        nxt_node = None
        curr_node = None
        last_gate_dir = dict()
        while True:
            if nxt_node is None:
                curr_node = qubit_node
                curr_node = graph_object.get_nxt_qubit_node(curr_node, qubit_num)
            else:
                curr_node = nxt_node
            if curr_node is None and nxt_node is None:
                break
            nxt_node = graph_object.get_nxt_qubit_node(curr_node, qubit_num)
            if graph_object.node_info[curr_node].is_sink():
                break
            if node_status[curr_node].get_status() == "Done":
                continue
            # if len(graph_object.get_qubit(curr_node)) > 1:
            #     # print(curr_node)
            #     # print(node_status[curr_node])
            #     curr_node_qubits = tuple(self.get_qubits(curr_node))
            #     # controls = graph_object.get_control(curr_node)
            #     # controls.sort()
            #     # curr_node_qubits = tuple([* controls, node_status[curr_node].get_target()])
            # else:
            curr_node_qubits = tuple(graph_object.get_qubit(curr_node))
            if curr_node_qubits not in last_gate_dir:
                node_status[curr_node].set_status("Seen")
                last_gate_dir[curr_node_qubits] = curr_node
            elif curr_node_qubits in last_gate_dir and \
                np.allclose(graph_object.get_matrix(curr_node) @ graph_object.get_matrix(last_gate_dir[curr_node_qubits]), np.eye(3)) and\
                    set(graph_object.get_control(last_gate_dir[curr_node_qubits])) == set(graph_object.get_control(curr_node)) and \
                        graph_object.get_target(curr_node) == graph_object.get_target(curr_node) and \
                            graph_object.get_all_control_value(curr_node) == \
                                graph_object.get_all_control_value(last_gate_dir[curr_node_qubits]):
                # Change status
                # print(node_status[last_gate_dir[curr_node_qubits]])
                # print(node_status[last_gate_dir[curr_node_qubits]])
                if node_status[last_gate_dir[curr_node_qubits]].get_status() == "Seen":
                    del last_gate_dir[curr_node_qubits]
                elif node_status[last_gate_dir[curr_node_qubits]].get_status() == "Done":
                    last_gate_dir[curr_node_qubits] = curr_node
                else:
                    assert False    
                continue
            elif curr_node_qubits in last_gate_dir and \
                not np.allclose(node_status[curr_node].get_matrix() @ node_status[last_gate_dir[curr_node_qubits]].get_matrix(), np.eye(3)) and\
                    set(graph_object.get_control(last_gate_dir[curr_node_qubits])) == set(graph_object.get_control(curr_node)) and \
                        graph_object.get_target(curr_node) == graph_object.get_target(curr_node) and \
                            graph_object.get_all_control_value(curr_node) == \
                                graph_object.get_all_control_value(last_gate_dir[curr_node_qubits]):
                if node_status[last_gate_dir[curr_node_qubits]].get_status() == "Done":
                    node_status[curr_node].set_status("Seen")
                    last_gate_dir[curr_node_qubits] = curr_node
                    continue
                assert node_status[last_gate_dir[curr_node_qubits]].get_status() == "Seen"
                # If these nodes are neighbours we shouldnt do this optimization
                dont_optimize = False
                new_received_gate = []
                for nodes in graph_object.dag.predecessors(curr_node):
                    if nodes == last_gate_dir[curr_node_qubits]:
                        dont_optimize = True
                        break
                if not dont_optimize: 
                    new_received_gate = graph_object.insert_identity_gate(last_matrix=node_status[last_gate_dir[curr_node_qubits]].get_matrix(),
                                                                          location=curr_node)
                    # print("INSERTING NODE ", curr_node)
                    node_status[new_received_gate[0]] = graph_object.node_info[new_received_gate[0]]
                    node_status[new_received_gate[1]] = graph_object.node_info[new_received_gate[0]]
                    node_status[new_received_gate[0]].set_status("Done")
                    node_status[new_received_gate[1]].set_status("Done")
                node_status[last_gate_dir[curr_node_qubits]].set_status("Done")
                last_gate_dir[curr_node_qubits] = curr_node
                node_status[curr_node].set_status("Seen")
    return graph_object

if __name__ == "__main__":
    # Get the non lifted dot file
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--name", help="name of the circuit to be lifted in this")
    # args = parser.parse_args()
    # dot = pydot.graph_from_dot_file(args.name + ".dot")[0]
    # # Get the graph
    # networkx_graph = nx.nx_pydot.from_pydot(dot)
    # info = GraphInfo()
    # info.create_node_directory(networkx_graph)
    # # Need to add support to see if a graph has not been lifted
    # # circ1, _ = info.build_circuit()
    # print("Start lifting")
    # lift_dag_circuit(info)
    print("main has been removed")
