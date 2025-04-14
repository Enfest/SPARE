from rewriter.src.nodeset import NodeSet
from rewriter.src.compile import RewriteEngine
from rewriter.src.compile_basic import gateslices_simplify, delete_projected_gates, verify_gateslices
import numpy as np
from src.toffoli import *

gate_name_history = []

def is_common_phase_present(self):
    slice_max, edge_max = self.return_iterator()
    return_val = False
    for i in range(1, slice_max-1):
        return_val = return_val or self.get_slice(i).is_common_phase_present(gateslice_handle=self, index=i)
    return return_val
    
def unroll_gate_common_phases(self):
    index = 0 
    while self.is_common_phase_present():
        slice_max, edge_max = self.return_iterator()
        self.slice_edge_collection.update_all_mappings()
        gate_pairs_to_add = []
        for i in range(1, slice_max-1):
            gate_pairs_to_add.append(self.get_slice(i).unroll_gate_common_phases(gateslice_handle=self, index=i))
        for gate_pairs in gate_pairs_to_add:
            for gate_pair in gate_pairs:
                self.add_a_gate_after_gate(gate_pair[1], gate_pair[0])
        self.slice_edge_collection.update_all_mappings()
        self.merge_greedily()
        slice_max_new, edge_max_new = self.return_iterator()
        index += 1

def apply_mergeswap_ops(self, target, debug):
    assert target["result"][1][0] == "ControlSwapper"
    gates = target["targets"]
    if not gates[0].get_node_type().is_gateset():
        assert not gates[1].get_node_type().is_gateset()
        frst_gate = gates[0]
        snd_gate = gates[1]
        frst_node = target["orig"]
        snd_node = snd_gate
        interface_gates, interface_gatesT = (None, None)
    else:
        assert gates[1].get_node_type().is_gateset()
        interface_gates = gates[0].getNodeSet()[0]
        interface_gatesT = gates[0].getNodeSet()[2]
        frst_gate = gates[0].getActiveSet()[0]
        snd_gate = gates[1].getActiveSet()[0]
        frst_node = target["orig"]
        snd_node = gates[1]
    
    if debug:
        print("target is ", target)
    all_mismatches = target["result"][1][1]
    global gate_name_history
    start = True
    sample_bookend = snd_gate.deepcopy(use_old_name=True)
    bookend_lists = [[], []]
    all_controls = frst_gate.get_control()
    while len(all_mismatches) > 1:
        mismatches = all_mismatches[-2:]
        assert len(mismatches) == 2
        all_mismatch_controls = np.array([[frst_gate.get_control_value_for_qubit(x) for x in mismatches],
                                          [snd_gate.get_control_value_for_qubit(x) for x in mismatches]]).T
        mismatch_controls = np.sum(all_mismatch_controls, axis=1)
        if mismatch_controls[0] == 2 or mismatch_controls[0] == 1:
            control = mismatch_controls[0]
        elif mismatch_controls[0] == 3:
            control = 1
        else:
            raise Exception
        for a in range(all_mismatch_controls.shape[1]):
            if all_mismatch_controls[0, a] == control:
                other_control = all_mismatch_controls[1, a]
                break
        other_control = mismatch_controls[1] - other_control
        # TODO: Change this workaround
        if mismatch_controls[1] == 2:
            matrix = cirq.unitary(QutritX02Gate())[:self.dimension, :self.dimension]
        elif mismatch_controls[1] == 1:
            matrix = cirq.unitary(QutritX01Gate())[:self.dimension, :self.dimension]
        elif mismatch_controls[1] == 3:
            matrix = cirq.unitary(QutritX12Gate())[:self.dimension, :self.dimension]
        else:
            raise Exception
        assert np.allclose(matrix @ matrix.T, np.eye(self.dimension))
        new_gate1 = sample_bookend.deepcopy()
        new_gate1.set_control_value_for_qubit(mismatches[0], control)
        new_gate1.delete_last_qubit(sample_bookend.get_target(), mismatches[1],
                                    new_matrix=matrix)
        for q in all_controls:
            if q not in mismatches:
                new_gate1.delete_qubit(q)
        if (new_gate1.get_qubits(), new_gate1.get_control_values()) not in gate_name_history:
            gate_name_history.append((new_gate1.get_qubits(), new_gate1.get_control_values()))
        new_gate1.setName("A"+str(len(gate_name_history)-1))
        new_gate4 = new_gate1.return_inverse()
        sample_bookend.delete_last_qubit(sample_bookend.get_target(), mismatches[1],
                                         new_matrix=matrix)
        bookend_lists[0].append(new_gate1)
        bookend_lists[1] = [new_gate4, * bookend_lists[1]]
        if start:
            new_gate2 = snd_gate.deepcopy()
            new_gate3 = snd_gate.deepcopy()
            new_gate2.set_control_value_for_qubit(mismatches[1], other_control)
            new_gate3.set_control_value_for_qubit(mismatches[1], other_control)
            new_gate2.setName("Ga"+str(len(gate_name_history)-1))
            new_gate3.setName("Gb"+str(len(gate_name_history)-1))
        else:
            new_gate2.set_control_value_for_qubit(mismatches[1], other_control)
            new_gate3.set_control_value_for_qubit(mismatches[1], other_control)
        start = False
        all_mismatches.pop(-1)
        all_controls.remove(mismatches[1])
    new_gate2.delete_control(mismatches[0])
    new_gate3.set_control_value_for_qubit(mismatches[0], 3-mismatch_controls[0])
    if interface_gates is None:
        gates = [bookend_lists[0], [new_gate2, new_gate3], bookend_lists[1]]
    else:
        gates = [[* interface_gates, *bookend_lists[0]], [new_gate2, new_gate3], [*bookend_lists[1], * interface_gatesT]]
    self.add_a_gate_after_gate(NodeSet(nodes=gates), snd_node)
    self.delete_gate(frst_node)
    self.delete_gate(snd_node)
    return True

def apply_rule_gate_merge(self, target):
    '''Merges gates found using gate merge function'''
    def get_gate_modifications(elements, control_vals, qtrt, dim=3):
        result_array = []
        element = ("Delete_control", [qtrt])
        result_array.append(element)
        if len(elements) == dim:
            result_array.append(("Modify_phase_to", elements[-1].get_matrix() @ elements[0].get_matrix().conjugate().T))
            return result_array
        elif len(elements) == dim - 1:
            result_array.append(("Modify_control_to", list(set(range(dim)).difference(control_vals))[0]))
            result_array.append(("Modify_phase_to", elements[-1].get_matrix() @ elements[0].get_matrix().conjugate().T))
            return result_array
        assert False
    
    positions = (self.find_gate(target["orig"])[1], self.find_gate(target["targets"][1])[1])
    result, nodes, orig_gate, gateids = (target["result"], target["targets"], target["orig"], target["positions"])
    if nodes[0].get_node_type().is_gateset():
        assert nodes[1].get_node_type().is_gateset()
        gates = (nodes[0].getActiveSet()[gateids[0]], nodes[1].getActiveSet()[gateids[1]])
    else:
        gates = nodes
    assert result[0] and result[1][0] == "GateMerge"
    mismatches = target["result"][1][1]
    assert len(mismatches) == 1
    all_mismatch_controls = [gates[0].get_control_value_for_qubit(mismatches[0]), gates[1].get_control_value_for_qubit(mismatches[0])]
    result_arr = get_gate_modifications(gates, all_mismatch_controls, qtrt=mismatches, dim=self.dimension)
    assert len(all_mismatch_controls) == len(result_arr)
    new_gates = []        
    for i in range(len(result_arr)):
        if result_arr[i][0] == "Modify_phase_to":
            # Need a wrapper with gateslice for this
            gates[1].set_matrix(result_arr[i][1])    
            new_gates.append(gates[1].deepcopy(use_old_name=True))
        elif result_arr[i][0] == "Modify_control_to":
            invctrlgate = gates[0].return_inverse(use_old_name=True)
            invctrlgate.set_control_for_qubit(mismatches[0], result_arr[i][1])
            new_gates.append(invctrlgate)
        elif result_arr[i][0] == "Delete_control":
            delctrlgate = gates[0].deepcopy(use_old_name=True)
            delctrlgate.delete_control(mismatches[0])
            new_gates.append(delctrlgate)
    self.add_gate_list_after_node(new_gates, nodes[1], gateids[1])
    self.delete_gate(orig_gate)
    self.prune_null_slices()
    return True

def apply_swapdel_ops(self, target, debug=False):
    '''This operation should return the swapped and modified gate maybe we remove the return gateslice case'''
    result, targets, orig_gate, gateids = (target["result"], target["targets"], target["orig"], target["positions"])
    positions = (self.find_gate(target["orig"])[1], self.find_gate(target["targets"][1])[1])
    assert result[0] and result[1][0] == "GateMatch"
    start_gate = targets[0]
    end_gate = targets[1]
    running_node = orig_gate.deepcopy(use_old_name=True)
    if positions[1] > positions[0]:
        start, end, diff = self.find_gate(orig_gate)[1], self.find_gate(end_gate)[1], 1
    else:
        start, end, diff = self.find_gate(orig_gate)[1], self.find_gate(end_gate)[1], -1
    for i in range(start+diff, end, diff):
        for gate in self.slice_collection[i].get_all_nodes():
            if gate.get_node_type().is_gate():
                running_node = gate.swap_and_modify_gates(running_node, forward_direction=positions[1] > positions[0])
    if positions[0] > positions[1]:
        targets[1].set_matrix(targets[0].get_matrix() @ targets[1].get_matrix())
    else:
        targets[1].set_matrix(targets[1].get_matrix() @ targets[0].get_matrix())
    self.delete_gate(orig_gate)
    self.prune_null_slices()
    return True
    
def apply_rule_broadcasts(self, target):
    '''Broadcasts a gate in the circuit'''
    target, target_bookends = target
    assert target.get_node_type().is_gate()
    target_copy = target.deepcopy()
    position = self.find_gate(target)[1]
    qubits_not_controlled = []
    for book_pair0 in target_bookends:
        qubits_not_controlled = list(set(qubits_not_controlled).union(set(book_pair0.get_qubits())))
    qubits_not_controlled = sorted(list(set(qubits_not_controlled) - set(target.get_qubits())))
    # assert False
    assert isinstance(qubits_not_controlled, list)
    # print("qubits not controlled - ", qubits_not_controlled)
    # qubits_not_controlled = list(range(self.slice_collection[0].qubits))
    gateList = target.broadcastQubits(qubits_not_controlled)
    for gate in gateList:
        self.add_a_gate_after_gate(gate, target)
    self.delete_gate(target)
    possible_sites = self.get_all_swapdel_sites(debug=False)  # True)
    bookends_deleted = False
    while len(possible_sites) > 0:
        for sites in possible_sites:
            if sites["orig"] not in gateList:
                bookends_deleted = True
                self.apply_swapdel_ops(sites)
        possible_sites = self.get_all_swapdel_sites(debug=False)  # True)
    #print("broadcasts done ", bookends_deleted)
    #print(self.build_circuit()[0])
    self.slice_edge_collection.update_all_mappings()
    return bookends_deleted
    
def apply_rule_1(self, target):
    start = target[0]
    end = target[1]
    target_list = target[2]
    assert end.is_synchronized(start)
    assert len(start.get_control()) > 0
    assert len(target_list) > 0
    for g in target_list:
        target_qubits = list(set(g.get_control()).intersection(set(start.get_control())))
        # Need to rempve protected controls from this target list
        for q in target_qubits:
            if g.get_rule_3_status_qubit(q) == "p":
                target_qubits.remove(q)
        self.delete_controls_on_gate(g, target_qubits)
    return True

def apply_rule_2(self, target):
    start = target[0]
    end = target[1]
    target_list = target[2]
    assert end.is_synchronized(start)
    assert len(start.get_control()) > 0
    assert len(target_list) > 0
    # Need to modify gates that have a target on this qubit
    running_matrix = start.get_matrix()    
    for g in target_list:
        if g.get_target() in start.get_control():
            pass
        elif g.get_target() == start.get_target():
            g.set_matrix(running_matrix.conjugate().T @ g.get_matrix() @ running_matrix)
        elif start.get_target() in g.get_control():
            if g.get_control_value_for_qubit(start.get_target()) == 2:
                for i in list(set(start.get_control()).intersection(g.get_control())):
                    # Set to enforced TODO: Make this redundant
                    g.set_rule_3_status("e", i)
            self.control_modify(g, running_matrix, start.get_target())
            # is this still 2
            if g.get_control_value_for_qubit(start.get_target()) == 2:
                self.delete_gate(g)
    if np.allclose(end.get_matrix() @ start.get_matrix(), np.eye(3)):
        self.delete_gate(start)
        self.delete_gate(end)
    else:
        end.set_matrix(end.get_matrix() @ start.get_matrix())
        self.delete_gate(start)
    self.slice_edge_collection.update_all_mappings()    
    return True

def apply_rule_3(self, target):
    start = target[0]
    end = target[1]
    end_copy = end.deepcopy()
    matrix2 = end.get_matrix() @ start.get_matrix()
    end.set_matrix(start.get_matrix().conjugate().T)
    end_copy.set_matrix(matrix2)
    self.add_a_gate_after_gate(end_copy, end)
    target_list = target[2]
    assert end.is_synchronized(start, check_status=True)
    assert len(start.get_control()) > 0
    assert len(target_list) > 0
    assert "e" in end.get_rule_3_status()
    assert np.allclose(np.eye(3), end.get_matrix() @ start.get_matrix())
    # Need to modify gates that have a target on this qubit
    control = end.get_control()
    rule_3_status = end.get_rule_3_status()
    i = rule_3_status.index("e")
    c = control[i]
    for g in target_list:
        if i in g.get_control() and g.get_target() != start.get_target():
            g.set_rule_3_status("p", c)
        else:
            pass
    self.delete_controls_on_gate(start, [c])
    self.delete_controls_on_gate(end, [c])    
    return True

def parse_and_merge_different_control_gates(self, if_verify, unitary1, by_one=False):
    results, controls, controls_vals = self.get_all_mergable_gates()
    while (len(results)) > 0:
        for i, r in enumerate(results):
            self = apply_rule_gate_merge(self, r, controls[i], controls_vals[i])
            self.merge_greedily()
            self.slice_edge_collection.set_gate_slice_collection(self)
            if by_one:
                break
        self.merge_greedily()
        self.get_all_matrix_mergable_gates()
        verify_gateslices(self, if_verify, unitary1)
        self = delete_projected_gates(self, debug=False)
        self.control_rewrite()
        results, controls, controls_vals = self.get_all_mergable_gates()
        verify_gateslices(self, if_verify, unitary1)
    return True

def modify_middle_bookends(self, target):
    result, targets, orig_gate, gateids, CtrandGates = target
    for ctr, gates in CtrandGates:
      new_gate = gates.deepcopy()
      add_gate = False
      for c in ctr:
        if c["ctr"] not in new_gate.get_qubits():
          new_gate.add_control(c["ctr"], c["val"])
          add_gate = True
      if orig_gate.get_target() in new_gate.get_control():
        new_gate.delete_qubit(orig_gate.get_target())
      if add_gate:
        self.add_a_gate_after_gate(new_gate, gates)
    self.delete_gate(orig_gate)
    self.delete_gate(targets[1])
    self.slice_edge_collection.update_all_mappings()
    return True

def partial_merge_nodesets(self, target):
    result, targets, orig_gate, gateids = (target["result"], target["targets"], target["orig"], target["positions"])
    positions = (self.find_gate(target["orig"])[1], self.find_gate(target["targets"][1])[1])
    if positions[0] < positions[1]:
        mismatched_bookends = targets[0].getMismatchedBookends(targets[1], forward_dir=True)
        targets[1].mergeMismatchedBookends(targets[0], mismatched_bookends, gateids)
    elif positions[0] > positions[1]:
        mismatched_bookends = targets[1].getMismatchedBookends(targets[0], forward_dir=False)
        targets[1].mergeMismatchedBookends(targets[0], mismatched_bookends, gateids)
    self.delete_gate(target["orig"])
    return True    

def optimizer_pass(self, debug=False):
    optimization_possible = True
    counter = 0    
    c, _ = self.build_circuit(breakdown="vvv")
    while optimization_possible:
        resulting_merges = self.get_all_swapdel_sites(debug=False)
        for r in resulting_merges:
            self.apply_swapdel_ops(r, debug=False)
        while True:
            resulting_merges = self.get_all_mergable_gates(debug=debug)
            for r in resulting_merges:
                self.apply_rule_gate_merge(r)
            self = delete_projected_gates(self, debug=False)
            self.slice_edge_collection.set_gate_slice_collection(self)
            if len(resulting_merges) == 0:
                break
        while True:
            resulting_merge_swaps = self.get_all_mergeswap_gates(debug=False)
            for r in resulting_merge_swaps:
                self.apply_mergeswap_ops(r, debug=False)
            self = delete_projected_gates(self, debug=False)
            break
        self.prune_null_slices()
        counter += 1  
        status1 = len(self.get_all_swapdel_sites()) > 0
        status2 = len(self.get_all_mergable_gates()) > 0
        status3 = len(self.get_all_mergeswap_gates()) > 0
        optimization_possible = status1 or status2 or status3 
    
    self.slice_edge_collection.update_all_mappings()
    self = gateslices_simplify(self)
    self.get_all_matrix_mergable_gates()    
    self.merge_greedily()
    while True:
        partial_merges = self.get_all_partial_merges()
        for m in partial_merges:
            self.partial_merge_nodesets(m)
        # print("partial merges")
        if len(partial_merges) == 0:
            break
    self.slice_edge_collection.update_all_mappings()
    return True


RewriteEngine.is_common_phase_present = is_common_phase_present
RewriteEngine.unroll_gate_common_phases = unroll_gate_common_phases
RewriteEngine.apply_mergeswap_ops = apply_mergeswap_ops
RewriteEngine.apply_rule_gate_merge = apply_rule_gate_merge
RewriteEngine.apply_swapdel_ops = apply_swapdel_ops
RewriteEngine.apply_rule_broadcasts = apply_rule_broadcasts
RewriteEngine.partial_merge_nodesets = partial_merge_nodesets
RewriteEngine.apply_rule_1 = apply_rule_1
RewriteEngine.apply_rule_2 = apply_rule_2
RewriteEngine.apply_rule_3 = apply_rule_3
RewriteEngine.optimizer_pass = optimizer_pass
RewriteEngine.modify_middle_bookends = modify_middle_bookends
