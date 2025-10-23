from rewriter.utils.check_subspace_equality import *
from rewriter.simple_circuit_rewrites import *
from rewriter.utils.node_gate_types import *
from rewriter.utils.edge_types import EdgeType

def merge_slices(gateslice):
    gateslice.merge_greedily()
    return gateslice

def delete_projected_gates(gateslice, debug=False):
    '''Modify this code to use the new syntactic sugar for discovering gates?'''
    def swap_with_gate(running_gate, g2, forward_dir=False):
        if g2.is_control_subset(running_gate) and g2.get_target() in running_gate.get_control():
            if g2.get_gate_type().state_permutation_gates():
                copy_temp = g2.return_inverse() if forward_dir else g2.deepcopy()
                running_gate.control_update_qubit(copy_temp, copy_temp.get_target())
            else:
                return False, running_gate
        elif g2.is_control_subset(running_gate) and g2.get_target() == running_gate.get_target():
            return True, running_gate
        elif g2.is_control_mismatch(running_gate) or running_gate.is_control_mismatch(g2):
            return True, running_gate
        elif not g2.is_control_subset(running_gate) and g2.get_target() in running_gate.get_control():
            if not forward_dir and not g2.is_immutable(running_gate.get_control_value_for_qubit(g2.get_target())):
                return False, running_gate
            return False, running_gate
        return True, running_gate
    
    def swap_with_node(gateslice, running_gate, j, forward_dir=False):
        offset = 0 if not forward_dir else -1
        if 2 in running_gate.get_control_values():
            for q in running_gate.get_control():
                if running_gate.get_control_value_for_qubit(q) == 2 and \
                    gateslice.slice_edge_collection.get_edge_slice(j + offset).get_wire_status(q) == EdgeType.Basis01:
                    return False, running_gate
        for g2 in gateslice.get_slice(j).get_all_nodes():
            if g2.get_node_type().is_qubit():
                return True, running_gate
            elif g2.get_node_type().is_gateset():
                can_swap = True
                for g2_ in g2.getActiveSet():
                    res, running_gate = swap_with_gate(running_gate, g2_, forward_dir=forward_dir)
                    can_swap = can_swap and res
            elif g2.get_node_type().is_gate():
                can_swap, running_gate = swap_with_gate(running_gate, g2, forward_dir=forward_dir)
            if not can_swap:
                return False, running_gate
        return True, running_gate

    gateslice.slice_edge_collection.update_all_mappings()
    slice_max, edge_max = gateslice.return_iterator()
    for i in range(1, slice_max-1):
        gates = gateslice.get_slice(i).get_all_nodes()
        for g in gates:
            cant_swap_to_boundary = False
            delete_gate = False
            if g.get_node_type().is_gateset():
                g.simplify_active_gates(gateslice.slice_edge_collection.get_edge_slice(i-1),
                                        gateslice.slice_edge_collection.get_edge_slice(i))
                continue
            running_gate = g.deepcopy()
            for j in reversed(range(i)):
                result, running_gate = swap_with_node(gateslice, running_gate, j)
                if not result and 2 in running_gate.get_control_values():
                    gateslice.delete_gate(g)
            if 2 in running_gate.get_control_values():
                gateslice.delete_gate(g)
            running_gate = g.deepcopy()
            cant_swap_to_boundary = False
            for j in range(i+1, slice_max-1):
                result, running_gate = swap_with_node(gateslice, running_gate, j, forward_dir=True)
                if not result and 2 in running_gate.get_control_values():
                    gateslice.delete_gate(g)
            if 2 in running_gate.get_control_values():
                gateslice.delete_gate(g)
    gateslice.slice_edge_collection.update_all_mappings()
    return gateslice

def verify_gateslices(gateslices, if_verify, unitary1):
  assert gateslices == gateslices.slice_edge_collection.gate_slice_collection
  if if_verify:
    print("ACTUALLY VERIFYING")
    c, _ = gateslices.build_circuit(check_circ=True, breakdown="vvv")
    unitary2 = get_circuit_unitaries(c, dimension=gateslices.dimension)
    check_circuit_equality(unitary1, unitary2, dimension=gateslices.dimension,
                           qubits_considered=[x for x in gateslices.slice_collection[0].node_map if not x._is_ancilla])
  else:
    return
  gateslices.slice_edge_collection.set_gate_slice_collection(gateslices)

def collapse_states(gateslice, if_verify=True, unitary=None):
    def basic_pass(gateslice, pass_name, if_verify, unitary1):
        res = False
        newdepth = None
        depth, _ = gateslice.return_iterator()
        while newdepth != depth or res:
            if pass_name == "collapse_states":
                gateslice.collapse_states()
            elif pass_name == "rem_high_ctrl":
                gateslice.remove_high_control_gates()
            verify_gateslices(gateslice, if_verify, unitary1)
            gateslice.merge_greedily()
            if newdepth is not None:
                depth = newdepth
            newdepth, _ = gateslice.return_iterator()
            print("\tverify current parse")
            verify_gateslices(gateslice, if_verify, unitary1)
    
    c, _ = gateslice.build_circuit()
    if if_verify:
        unitary1 = get_circuit_unitaries(c, dimension=gateslice.dimension)
        if unitary is not None:
            a, b = check_circuit_equality(unitary1, unitary, dimension=gateslice.dimension,
                                          qubits_considered=[x for x in gateslice.slice_collection[0].node_map if not x._is_ancilla])
    else:
        unitary1 = None
    basic_pass(gateslice, "collapse_states", if_verify, unitary1)
    gateslice.slice_edge_collection.update_all_mappings()
    basic_pass(gateslice, "rem_high_ctrl", if_verify, unitary1)    
    gateslice.prune_null_slices()
    gateslice.control_rewrite()
    gateslice.slice_edge_collection.update_all_mappings()
    basic_pass(gateslice, "collapse_states", if_verify, unitary1)
    gateslice.slice_edge_collection.update_all_mappings()
    for _ in range(2):
        gateslice = delete_projected_gates(gateslice, debug=False)
        verify_gateslices(gateslice, if_verify, unitary1)
        gateslice.get_all_matrix_mergable_gates(debug=False)
        verify_gateslices(gateslice, if_verify, unitary1)
        c, _ = gateslice.build_circuit()
    '''TODO: Figure out how to best prune phases of the computation '''
    gateslice.merge_greedily()
    for _ in range(2):
        gateslice.get_all_matrix_mergable_gates(debug=False)
    verify_gateslices(gateslice, if_verify, unitary1)
    return gateslice

def gateslices_simplify(info, debug=False):
    if debug:
        info.print()
    assert info == info.slice_edge_collection.gate_slice_collection
    info.remove_high_control_gates()
    info.remove_high_phase_gates()
    info.prune_null_slices()
    info.merge_greedily()
    info.slice_edge_collection.update_all_mappings(debug=debug)
    return info
