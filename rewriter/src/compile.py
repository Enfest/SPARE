import numpy as np
# from rewriter.src.node import Node
from rewriter.src.sliceedges import EdgeCollection, EdgeSliceNodes
from rewriter.src.slice import GateSlice
from rewriter.src.slicecollection import GateSliceCollection
from src.toffoli import *
from rewriter.utils.graph_utils import *
from rewriter.utils.edge_types import EdgeType
import json
from rewriter.utils.node_gate_types import *
import collections

class RewriteEngine(GateSliceCollection):
    def serialize(self, file_path):
        data = {
            "qubits": self.num_qubits,
            "slices_len": len(self.slice_collection),
            "slices": self.slice_collection_serialize(),
            "slice_collection": self.slice_edge_collection.serialize()
        }
        result  = json.dumps(data, indent=4)
        with open(file_path, "w") as json_file:
            json_file.write(result)
        return
    
    @classmethod
    def read_in(cls, json_name, dimension):
        with open(json_name, "r") as json_file:
            data = json.load(json_file)
        num_qubits = data.get("qubits")
        len_slices = data.get("slices_len")
        slices_data = data.get("slices")
        # print(slices_data)
        slice_list = []
        for i in range(len_slices):
            slice_list.append(GateSlice.read_in(slices_data.get(str(i)), dimension=dimension))
        slice_collection_data = data.get("slice_collection")
        slice_collection_edges = EdgeSliceNodes.read_in(slice_collection_data)
        return cls(None, num_qubits, {"slice_collection":slice_list,
                                      "slice_edge_collection":slice_collection_edges,
                                      "num_qubits":num_qubits, "dimension":dimension})
        
    def deepcopy(self):
        '''Save and prepare a deepcopy of a slicecollection object'''
        self.serialize("temp_new1.json")
        return RewriteEngine.read_in("temp_new1.json", self.dimension)

    def slice_collection_serialize(self):
        data = collections.defaultdict()
        for i in range(len(self.slice_collection)):
            data[i] = self.slice_collection[i].serialize()
        return data
    
    def get_rule_candidates(self, idx, if_inverse=False, phase_ok=False, check_status=False, check_control=False):
        gates = self.slice_collection[idx].get_all_nodes()
        res = []
        for g in gates:
            if g.get_node_type().is_gate() and not g.get_node_type().is_gateset():
                if len(g.get_control()) > 0 and \
                    ((self.slice_edge_collection.get_edge_slice(idx-1).get_wire_status(g.get_target()) == EdgeType.Basis01 and \
                    self.slice_edge_collection.get_edge_slice(idx).get_wire_status(g.get_target()) != EdgeType.Basis01) or check_control):
                    second = self.get_rule_end(idx, g, if_inverse=if_inverse, phase_ok=phase_ok,
                                            check_status=check_status, check_basis=not check_control)
                    if not isinstance(second, bool):
                        assert self.find_gate(g)[1] < self.find_gate(second)[1]
                        res.append((g, second))
        return res
    
    def get_all_candidates(self, idx):
        gates = self.slice_collection[idx].get_all_nodes()
        res = []
        for g in gates:
            second = self.get_rule_end(idx, g, if_inverse=False, phase_ok=True, check_status=False, check_basis=False)
            if not isinstance(second, bool):
                res.append((g, second))
        return res

    def get_rule_end(self, idx, startgate, if_inverse=False, phase_ok=True, check_status=False, check_basis=True):
        slice_len = len(self.slice_collection)
        res = []
        running_node = startgate.deepcopy()
        for i in range(idx + 1, slice_len-1):
            gates = self.slice_collection[i].get_all_nodes(startgate.get_qubits())
            for g in gates:
                if g.get_node_type().is_gate() and not g.get_node_type().is_gateset():
                    if not if_inverse and g.is_synchronized(running_node):
                        if not check_basis or \
                            (self.slice_edge_collection.get_edge_slice(i-1).get_wire_status(g.get_target()) == EdgeType.Basis012 and \
                            self.slice_edge_collection.get_edge_slice(i).get_wire_status(g.get_target()) == EdgeType.Basis01):
                            return g
                    if if_inverse and g.is_inverse(running_node, phase_ok=phase_ok, check_status=check_status):
                        if not check_basis or \
                            (self.slice_edge_collection.get_edge_slice(i-1).get_wire_status(g.get_target()) == EdgeType.Basis012 and \
                            self.slice_edge_collection.get_edge_slice(i).get_wire_status(g.get_target()) == EdgeType.Basis01):
                            return g
        return False
    
    def get_all_rule_1_sites(self):
        slice_len = len(self.slice_collection)
        result_arr = []
        for i in range(1, slice_len-1):
            gates_pairs = self.get_rule_candidates(i, if_inverse=False, phase_ok=True)
            for pair in gates_pairs:
                real_gate_list = self.get_all_nodes_between(pair[0], pair[1])
                gate_list = self.get_all_nodes_between(pair[0], pair[1], qubits=pair[0].get_control())
                gate_list = self.remove_mutable_gates(pair, gate_list)
                gate_list = self.remove_inv_pairs(gate_list, qubits=pair[0].get_control(), phase_ok=True)
                res_arr = []
                dont_add = False
                for g in gate_list:
                    if g.get_target() in pair[0].get_control():
                        # check for immutability or inverse
                        dont_add = True
                        pass
                    if g.get_target() == pair[0].get_target():
                        pass
                for g in real_gate_list:
                    if g.get_target() not in pair[0].get_qubits() and \
                        len(list(set(g.get_control()).intersection(set(pair[0].get_control())))) > 0 \
                        and pair[0].get_target() in g.get_control() and g.get_control_value_for_qubit(pair[0].get_target()) == 2:
                        res_arr.append(g)
                if not dont_add and len(res_arr) > 0:
                    result_arr.append((pair[0], pair[1], res_arr))
        if len(result_arr) == 0:
            return None
        return result_arr

    def get_all_rule_2_sites(self, check_control=True):
        def is_basic_rewrite(gate):
            assert gate.get_node_type().is_gate()
            if gate.get_gate_type().state_permutation_gates():
                return True
            return False

        slice_len = len(self.slice_collection)
        result_arr = []
        for i in range(1, slice_len-1):
            gates_pairs = self.get_rule_candidates(i, if_inverse=False, phase_ok=True, check_control=check_control)
            for pair in gates_pairs:
                return_true = True
                if not is_basic_rewrite(pair[0]):
                    return_true = False
                orig_gatelist = self.get_all_nodes_between(pair[0], pair[1])
                gate_list = self.remove_mutable_gates(pair, orig_gatelist)
                gate_list = self.remove_inv_pairs(gate_list, qubits=pair[0].get_qubits())                
                for g in orig_gatelist:
                    if pair[0].get_target() in g.get_control() and (not g.is_control_superset(pair[0])) and\
                          (not pair[0].is_immutable(g.get_control_value_for_qubit(pair[0].get_target()))):
                        return_true = False
                if return_true:     
                    return_true2 = False   
                    for g in orig_gatelist:
                        if check_control and pair[0].get_control_value_for_qubit(pair[0].get_control()[0]) == 2:
                            orig_gatelist = self.get_all_nodes_between(pair[0], pair[1])
                            gate_list = self.get_all_nodes_between(pair[0], pair[1], qubits=pair[0].get_control()) # gate_list
                            gate_list = self.remove_mutable_gates(pair, gate_list)
                            gate_list = self.remove_inv_pairs(gate_list, qubits=pair[0].get_qubits())
                            running_node = pair[0].deepcopy()
                            for g in orig_gatelist:
                                if not running_node.is_control_subset(g) and running_node.get_target() in g.get_control() and \
                                    not running_node.is_immutable(g.get_control_value_for_qubit(running_node.get_target())):
                                    return_true2 = False
                                    break
                                elif running_node.is_control_subset(g) and running_node.get_target() == g.get_target():
                                    pass
                                elif running_node.is_control_subset(g) and running_node.get_target() in g.get_control() and not \
                                    running_node.is_immutable(g.get_control_value_for_qubit(running_node.get_target())):
                                    return_true2 = True
                        elif not check_control:
                            return_true2 = True

                        if self.get_if_nodes_mutable(pair[0], pair[1], gate_list, just_binary=True) and return_true2:
                            result_arr.append((pair[0], pair[1], orig_gatelist))
                            break
        if len(result_arr) == 0:
            return None
        return result_arr

    def get_all_rule_3_sites(self):
        slice_len = len(self.slice_collection)
        result_arr = []
        for i in range(1, slice_len-1):
            # print("rule 3 check ", i)            
            gates_pairs = self.get_rule_candidates(i, if_inverse=False, phase_ok=True, check_status=True)
            for pair in gates_pairs:
                if "e" in pair[0].get_rule_3_status() and "e" in pair[1].get_rule_3_status() and\
                    pair[0].get_control_value_for_qubit(pair[0].get_control()[0]) == 2: 
                    gate_list = self.get_all_nodes_between(pair[0], pair[1], qubits=pair[0].get_control())
                    orig_gatelist = self.get_all_nodes_between(pair[0], pair[1])
                    rem_list = []
                    for g in gate_list:
                        if g.get_target() == pair[0].get_target():
                            rem_list.append(g)
                    gate_list = list(set(gate_list).intersection(rem_list))
                    gate_list = self.remove_mutable_gates(pair, gate_list)
                    gate_list = self.remove_inv_pairs(gate_list, qubits=pair[0].get_qubits())
                    if self.get_if_nodes_mutable(pair[0], pair[1], gate_list, just_binary=True):
                        result_arr.append((pair[0], pair[1], orig_gatelist))
        if len(result_arr) == 0:
            return None
        return result_arr

    def get_all_swapdel_sites(self, debug=False):
        resulting_matches = []
        position_tracker = []
        all_swaps = []
        for i, sl in enumerate(self.slice_collection):
            for g1 in sl.get_all_nodes():
                if g1.get_node_type().is_gate():
                    result, gates, positions  = self.parse_for_matching_gates(g1.deepcopy(), i, forward_direction=True,
                                                                             swap_through_gates=True, debug=debug, desired_result=None)
                    result2, gates2, positions2 = self.parse_for_matching_gates(g1.deepcopy(), i, forward_direction=False,
                                                                                swap_through_gates=True, debug=False, desired_result=None) # "GateMatch")
                    if (result2[0] and not result[0]):
                        result = result2
                        gates = gates2
                        positions = positions2
                    elif (result[0] and result[1][0] != "GateMatch" and result2[0]):
                        result = result2
                        gates = gates2
                        positions = positions2
                    
                    if not (result[0]):
                        if debug:
                            print(f"\t{g1.get_qubits()}, break here, {result}, {gates}, {positions}")
                        break
                    elif result[1][0] == "GateMatch":
                        gate_added = False
                        for p in positions:
                            if p in position_tracker:
                                gate_added = True
                        if debug:
                            print(f"\t{g1.get_qubits()}, Gate Match, {result}, {gates[1].get_qubits()}, {positions}, {gate_added}")
                        if not gate_added:
                            resulting_matches.append({"result":result, "targets":gates, "orig":g1, "positions":result[2]})
                            for p in positions:
                                if p not in position_tracker:
                                    position_tracker.append(p)
                            all_swaps.append(positions)
        return resulting_matches

    def get_all_matrix_mergable_gates(self, debug=False):
        ''' 
        Find all gates that can be merged with respect to one control in one sweep
        '''
        gates_to_del = []
        for i, sl in enumerate(self.slice_collection):
            for g1 in sl.get_all_nodes():
                if g1.get_node_type().is_gate() and g1 not in gates_to_del:    
                    result, gates, positions = self.parse_for_matching_gates(g1.deepcopy(), i, forward_direction=True, debug=debug)
                    if not (result[0]):
                        # break
                        continue                    
                    elif result[1][0] == "GateMatch":
                        if debug:
                            print("-------------------------")
                        assert len(result[1][1]) == 0
                        gates_to_del.append(g1)
                        gates[1].set_matrix(gates[1].get_matrix() @ gates[0].get_matrix())
                        break
        for gates in gates_to_del:   
            self.delete_gate(gates)    
        return
    
    def any_nodeset_preset(self, debug=False):
        if debug:
          for i, s in enumerate(self.slice_collection):
            for g in s.get_all_nodes():
                if g.get_node_type().is_gateset():
                    print("gateset found at index ", i)
        #   print(self.build_circuit()[0])
    
    def get_all_bookend_broadcasts(self):
        '''Find all gates that do not have any possible mergeable gates'''
        all_targets = {}
        all_gates = []
        for i, s in enumerate(self.slice_collection):
          for g in s.get_all_nodes():
            gates_in_between = []
            middle_broadcast_gates = []
            if not g.get_node_type().is_gateset() and g.get_node_type().is_gate() and g not in all_gates:
              result, gates, positions, CtrandGates = self.parse_for_matching_gates(g.deepcopy(), i, forward_direction=True,
                                                                         swap_through_gates=True, invalid_swap=True, return_middle=True)
              if result[0]:
                key = len(CtrandGates)
                pair_and_middle_gate = (result, gates, g, result[2], CtrandGates)
                if key in all_targets:
                  all_targets[key].insert(0, pair_and_middle_gate)
                else:
                  all_targets[key] = [pair_and_middle_gate]
                all_gates.append(gates[1])
        result = all_targets[min(list(all_targets.keys()))] if len(all_targets.keys()) > 0 else []
        return all_targets[min(list(all_targets.keys()))] if len(all_targets.keys()) > 0 else []

    def get_all_broadcastable_gates(self):
        '''Find all gates that do not have any possible mergeable gates'''
        all_gates = []
        gate_pairs = []
        gates_of_interest = []
        bookend_pos = []
        for i, s in enumerate(self.slice_collection):
            for g in s.get_all_nodes():
                if not g.get_node_type().is_gateset() and g.get_node_type().is_gate() \
                    and g not in all_gates and g not in gate_pairs:
                  result, gates, positions = self.parse_for_matching_gates(g.deepcopy(), i, forward_direction=True,swap_through_gates=True, invalid_swap=True)
                  if result[0]:
                    gate_pairs.append(g)
                    gate_pairs.append(gates[1])
                    bookend_pos.append((self.find_gate(g)[1], self.find_gate(gates[1])[1]))
                    gates_of_interest.append(g)
                if not g.get_node_type().is_gateset() and g.get_node_type().is_gate():
                  all_gates.append(g)
        gates_not_in_pairs = [x for x in all_gates if x not in gate_pairs]
        gates_not_in_pairs = [x for x in gates_not_in_pairs if len(x.get_qubits()) != self.num_qubits]
        prune_list = []
        bookend_start_gates = []
        for g in gates_not_in_pairs:
            curr_bookend_of_itr = []
            self_pos = self.find_gate(g)[1]
            for i, p in enumerate(bookend_pos):
                if p[0] < self_pos and self_pos < p[1]:
                    curr_bookend_of_itr.append(gates_of_interest[i])
            if len(curr_bookend_of_itr) == 0:
                prune_list.append(g)
            else:
                bookend_start_gates.append(curr_bookend_of_itr)
        for g in prune_list:
            idx = gates_not_in_pairs.index(g)
            gates_not_in_pairs.pop(idx)
        for g in gates_not_in_pairs:
            assert g.get_node_type().is_gate()
        assert len(gates_not_in_pairs) == len(bookend_start_gates)
        return [(gates_not_in_pairs[i], bookend_start_gates[i]) for i in range(len(gates_not_in_pairs))]
    
    def get_all_mergable_gates(self, debug=False):
        '''Find all gates that have same control when swapped with respect to all gates in between:
        This pass goes in both directions and has repeated implementation'''
        resulting_matches = []
        position_tracker = []
        all_swaps = []
        for i, sl in enumerate(self.slice_collection):
            for g1 in sl.get_all_nodes():
                if g1.get_node_type().is_gate():
                    result, gates, positions  = self.parse_for_matching_gates(g1.deepcopy(), i, forward_direction=True,
                                                                              swap_through_gates=True, debug=debug, desired_result=None) # "GateMerge")
                    result2, gates2, positions2 = self.parse_for_matching_gates(g1.deepcopy(), i, forward_direction=False,
                                                                                swap_through_gates=True, debug=debug, desired_result=None) # "GateMerge")
                    if (result2[0] and not result[0]):
                        result, gates, positions = (result2, gates2, positions2)
                    elif (result[0] and result[1][0] != "GateMerge" and result2[0]):
                        result, gates, positions = (result2, gates2, positions2)
                    if not (result[0]):
                        break
                    elif result[1][0] == "GateMerge":
                        gate_added = False
                        for p in positions:
                            if p in position_tracker:
                                gate_added = True
                        if not gate_added:
                            resulting_matches.append({"result":result, "targets":gates, "orig":g1, "positions":result[2]})
                            for p in positions:
                                if p not in position_tracker:
                                    position_tracker.append(p)
                            all_swaps.append(positions)
        return resulting_matches
    
    def get_all_mergeswap_gates(self, debug=False):
        resulting_swaps = []
        position_tracker = []
        all_pos = []
        for i, sl in enumerate(self.slice_collection):
            for g1 in sl.get_all_nodes():
                if g1.get_node_type().is_gate():
                    result, gates, positions  = self.parse_for_matching_gates(g1.deepcopy(), i, forward_direction=True, desired_result=None, debug=debug)
                    result2, gates2, positions2 = self.parse_for_matching_gates(g1.deepcopy(), i, forward_direction=False, desired_result=None)
                    if (result2[0] and not result[0]):
                        result, gates, positions = (result2, gates2, positions2)
                    elif (result[0] and result[1][0] != "ControlSwapper" and result2[0]):
                        result, gates, positions = (result2, gates2, positions2)
                    if not (result[0]):
                        break
                    elif result[1][0] == "ControlSwapper":
                        positions = list(positions)
                        positions.sort()
                        already_added = False
                        for p in positions:
                            if p in all_pos:
                                already_added = True
                        if not already_added:
                            resulting_swaps.append({"result":result, "targets":gates, "orig":g1, "positions":result[2]})
                            position_tracker.append(positions)
                            all_pos.append(positions[0])
                            all_pos.append(positions[1])
        if debug:
            print("position track ", debug, position_tracker)
        all_pos2 = []
        for r in position_tracker:
            if r[0] in all_pos2 or r[1] in all_pos2:
                print("\t", r, " data is already present")
                assert False
            all_pos2.append(r[0])
            all_pos2.append(r[1])           
        return resulting_swaps
    
    def get_all_partial_merges(self, debug=False):
        resulting_swaps = []
        position_tracker = []
        all_pos = []
        for i, sl in enumerate(self.slice_collection):
            for g1 in sl.get_all_nodes():
                if g1.get_node_type().is_gate():
                    result, gates, positions  = self.parse_for_matching_gates(g1.deepcopy(), i, forward_direction=True, desired_result="PartialMerge", debug=debug)
                    result2, gates2, positions2 = self.parse_for_matching_gates(g1.deepcopy(), i, forward_direction=False, desired_result="PartialMerge")
                    if (result2[0] and not result[0]):
                        result, gates, positions = (result2, gates2, positions2)
                    elif (result[0] and result[1][0] != "PartialMerge" and result2[0]):
                        result, gates, positions = (result2, gates2, positions2)
                    if not (result[0]):
                        break
                    elif result[1][0] == "PartialMerge":
                        positions = list(positions)
                        positions.sort()
                        already_added = False
                        for p in positions:
                            if p in all_pos:
                                already_added = True
                        if not already_added:
                            resulting_swaps.append({"result":result, "targets":gates, "orig":g1, "positions":result[2]})
                            position_tracker.append(positions)
                            all_pos.append(positions[0])
                            all_pos.append(positions[1])
        all_pos2 = []
        return resulting_swaps