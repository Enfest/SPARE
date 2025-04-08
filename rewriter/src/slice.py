import numpy as np
from rewriter.src.nodeset import LabelledNode, NodeSet, control_config
from rewriter.src.sliceedges import EdgeCollection, EdgeSliceNodes
from src.toffoli import *
import cirq
from rewriter.utils.graph_utils import *
import collections


def renormalize_matrix(gate, dimension, force_convert=True, debug=False):
    if debug:
        print(force_convert, (not gate.get_node_type().is_gateset() and gate.get_node_type().is_gate()))
    if force_convert and (not gate.get_node_type().is_gateset() and gate.get_node_type().is_gate()):
        if debug:
            print(gate.get_matrix().shape, dimension)
        if gate.get_matrix().shape[0] < dimension:
            new_matrix = np.eye(dimension, dtype=np.complex128)
            new_matrix[:gate.get_matrix().shape[0], :gate.get_matrix().shape[1]] = gate.get_matrix()
            gate.set_matrix(new_matrix, assert_disabled=force_convert)
                    
class GateSlice:
    def __init__(self, qubits:int, dimension:int, gate=None, gatelist=None, force_convert=False, print_det=False):
        self.node_map = [None] * qubits
        self.slice_status = False
        self.dimension = dimension
        self.qubits = qubits
        if gate is not None and gatelist is None:
            for q in gate.get_qubits():
                if isinstance(q, int) and isinstance(self.node_map, list) and q >= len(self.node_map):
                    print("gate is ", gate, " qubits are ", len(self.node_map))
                    if print_det:
                        print(gate)
                        print(len(self.node_map))
                self.node_map[q] = gate
            renormalize_matrix(gate, dimension, force_convert)
                
        if gate is None and gatelist is not None:
            assert qubits == len(gatelist)
            for g in gatelist:
                if g is not None:
                    for q in g.get_qubits():
                        self.node_map[q] = g
                    renormalize_matrix(g, dimension, force_convert)
    
    def renormalize(self, dimension):
        for g in self.node_map:
            if g is not None and not g.get_node_type().is_gateset():
                renormalize_matrix(g, dimension, debug=False)
            elif g is not None and g.get_node_type().is_gateset():
                raise NotImplementedError

        self.dimension = dimension

    def __str__(self):
        assert self.node_map is not None
        c = cirq.Circuit()
        assert len(self.node_map) > 0
        c.append(self.realize(cirq.LineQid.range(len(self.node_map), dimension=self.dimension), dimension=self.dimension))
        print(c)
        return ""

    def realize(self, qubits, dimension=3):
        nodes_realized = []
        circ =  cirq.Circuit()
        print_null_slice = True
        string_print = []
        for g in self.node_map:
            if g is None:
                string_print.append("None")
            elif g.get_node_type().is_qubit():
                string_print.append(str(self.node_map.index(g)) + " anc: " + str(g.is_ancilla()))
            else:
                print_null_slice = False
        if print_null_slice:
            print(string_print)
        for g in self.node_map:
            if g not in nodes_realized and g is not None and g.get_node_type().is_gate():
                '''Realize nodes here'''
                circ.append(g.realizeNode(qubits, breakdown="vvv"))
                nodes_realized.append(g)
        return circ
                
    def serialize(self):
        data = collections.defaultdict()
        data["qubits"] = len(self.node_map)
        for i in range(len(self.node_map)):
            if self.node_map[i] is None:
                data["node_" + str(i)] = ""
            else:
                data["node_" + str(i)] = self.node_map[i].serialize()
        return data

    @classmethod
    def read_in(cls, data, dimension):
        gatelist = []
        for i in range(data["qubits"]):
            if data.get("node_" + str(i)) != "":
                if "gate data" in data.get("node_"+ str(i)):
                    gatelist.append(NodeSet.read_in(data.get("node_" + str(i))))
                else:
                    gatelist.append(LabelledNode.read_in(data.get("node_" + str(i))))
            else:
                gatelist.append(None)
        return cls(data["qubits"], dimension=dimension, gate=None, gatelist=gatelist) 

    def if_qubit_annotation(self):
        for q in self.node_map:
            if q is not None and q.get_node_type().is_qubit():
                pass
            else:
                return False
        return True

    def get_qubits(self):
        return range(len(self.node_map))
    
    def get_num_qubits(self):
        return len(self.node_map)
    
    def set_node(self, node):
        assert type(node) == LabelledNode
        for q in node.get_qubits():
            self.node_map[q] = node

    def remove_node(self, node):
        for q in node.get_qubits():
            self.node_map[q] = None

    def modify_control_on_gate(self, node):
        for q in self.node_map:
            if self.node_map[q] == node and q not in node.get_qubits():
                self.node_map[q] = None
    
    def get_node(self, qubit):
        assert isinstance(qubit, int)
        return self.node_map[qubit]
    
    def is_root_slice(self):
        if self.is_null_slice():
            return False
        for q in self.node_map:
            if q is not None:
                if not q.is_root():
                    return False
        return True

    def is_null_slice(self):
        gates_to_del = []
        for q in self.node_map:
            if q is not None and q.get_node_type().is_gateset() and q not in gates_to_del:
                if q.prune_if_null():
                    gates_to_del.append(q)
            elif q is not None and q.get_node_type().is_gate() and q not in gates_to_del:
                if np.allclose(q.get_matrix(), np.eye(self.dimension)):
                    gates_to_del.append(q)
        for g in gates_to_del:
            for q in g.get_qubits():
                self.node_map[q] = None            
        for q in self.node_map:
            if q is not None:
                if q.get_node_type().is_qubit():
                    return False
                elif q.get_node_type().is_gateset():
                    if len(q.getActiveSet()) > 0:
                        return False
                elif q.get_node_type().is_gate():
                    if not np.allclose(q.get_matrix(), np.eye(self.dimension)):
                        return False
                else:
                    assert False
        return True
    
    def is_sink_slice(self):
        if self.is_null_slice():
            return False
        for q in self.node_map:
            if q is not None:
                if not q.is_sink():
                    return False
        return True

    def get_all_nodes(self, qubits=None):
        result = []
        if qubits is not None:
            for q in qubits:
                if self.node_map[q] is not None and self.node_map[q] not in result:
                    result.append(self.node_map[q])
            return result
        for nodes in self.node_map:
            if nodes is not None and nodes not in result:
                result.append(nodes)
        return result
    
    def reroll_cnots(self, previous_maps, control_maps, previous_states):
        # other_gates = []
        # for i, n in enumerate(self.node_map):
        #     if n is not None and i == n.get_target():
        #         if n.get_node_type().is_gate() and not n.get_node_type().is_gateset() and len(n.get_control()) == 1:
        #             if previous_states[n.get_qubits()[1]].get_wire_mapping().is_0_basis() and previous_states[n.get_qubits()[0]].get_wire_mapping().is_binary_state():
        #                 if previous_maps[n.get_qubits()[1]] == n.get_qubits()[0]:
        #                     pass
        #                 elif previous_maps[n.get_qubits()[1]] == n.get_qubits()[1] and control_maps[n.get_qubits()[1]]:
        #                     previous_maps[n.get_qubits()[1]] = n.get_qubits()[0]
        #                     if control_maps[n.get_qubits()[1]] is None or control_maps[n.get_qubits()[1]] == n.get_control_values()[0]:
        #                         control_maps[n.get_qubits()[1]] = n.get_control_values()[0]
        #                     else:
        #                         self.maps_per_slice = None
        #                         return False, None
        #                 else:
        #                     self.maps_per_slice = None
        #                     return False, None
        #         elif n.is
        # self.control_maps = control_maps
        # self.maps_per_slice = previous_maps
        pass
        # return True, previous_maps, 
    
    def reset_ancilla_maps(self):
        self.maps_per_slice = None
        self.control_maps = None

    def parse_and_roll_cnots(self):
        assert self.control_maps is not None and self.maps_per_slice is not None
        gates_to_del = []
        for i, n in enumerate(self.node_map):
            if i == n.get_target():
                if n is not None and n.get_node_type().is_gate() and not n.get_node_type().is_gateset() and len(n.get_control()) == 1:
                    # Mark it identity for deletion
                    self.node_map[q].set_matrix(np.eye(self.get_matrix().shape[0]))
                elif n is not None:
                    if n.get_node_type().is_gate():
                        n.reset_qubits(self.maps_per_slice, self.control_maps)
                    else:
                        assert False
        self.reset_ancilla_maps()

    def if_slices_mergeable(self, prev_slice):
        assert len(self.node_map) == len(prev_slice.node_map)
        if prev_slice.get_slice_status():
            return False
        for i, n in enumerate(self.node_map):
            if n is None and prev_slice.get_node(i) is None:
                pass
            elif n is None and prev_slice.get_node(i).get_node_type().is_qubit():
                return False
            elif n is None and prev_slice.get_node(i).get_node_type().is_gate():
                pass
            elif n.get_node_type().is_qubit():
                if prev_slice.get_node(i) is not None and prev_slice.get_node(i).get_node_type().is_qubit():
                    pass
                elif prev_slice.get_node(i) is not None and prev_slice.get_node(i).get_node_type().is_gate():
                    return False
            elif n.get_node_type().is_gateset():
                if prev_slice.get_node(i) is not None:
                    if prev_slice.get_node(i).get_node_type().is_gateset():
                        return n.is_mergable(prev_slice.get_node(i))
                    else:
                        return False
            elif n.get_node_type().is_gate():
                if prev_slice.get_node(i) is not None:
                    if prev_slice.get_node(i).get_node_type().is_gate():
                        if not n.is_mergable(prev_slice.get_node(i)):
                            return False
                    else:
                         return False
        return True
    
    def get_slice_status(self):
        return self.slice_status
    
    def merge(self, prev_slice, binary=False):
        '''Merge happens for the slice with the gates coming in'''
        assert self.if_slices_mergeable(prev_slice=prev_slice)
        nodes_modified = []
        if self.if_qubit_annotation() or prev_slice.if_qubit_annotation():
            assert False
        for i, n in enumerate(self.node_map):
            if n is None and prev_slice.get_node(i) is not None:
                self.node_map[i] = prev_slice.get_node(i)
            elif n is None and prev_slice.get_node(i) is None:
                self.node_map[i] = None
            elif n is not None and prev_slice.get_node(i) is None:
                pass
            elif n.get_node_type().is_gateset():
                if n not in nodes_modified:
                    if prev_slice.get_node(i).get_node_type().is_gateset():
                        n.merge(prev_slice.get_node(i))
                    nodes_modified.append(n)
            elif n.get_node_type().is_gate():
                if n not in nodes_modified:
                    if prev_slice.get_node(i).get_node_type().is_gate():
                        n.set_matrix(n.get_matrix() @ prev_slice.get_node(i).get_matrix(), assert_disabled=binary)
                    nodes_modified.append(n)
        return
    
    def get_single_qubit_gates(self):
        result = []
        for q, g in enumerate(self.node_map):
            if q is not None and q.get_node_type().is_gate():
                if len(q.get_qubits()) == 1:
                    result.append((q, g))
        return result
    
    def update_gate_slice(self, single_qubit_gate):
        assert len(single_qubit_gate) == len(self.node_map)
        new_slice_to_add = [None] * len(self.node_map)
        valid = True
        q_collect = []
        # print("_______________________")
        for q in range(len(single_qubit_gate)):
            if self.node_map[q] is not None:
                single_qubit_gate[q], if_valid = self.node_map[q].swap_or_merge_gate(single_qubit_gate[q], self.dimension)
                # print("Validity check ", if_valid)
                valid = valid and if_valid
                if not if_valid:
                    q_collect.append(q)
        # print("q_collect ", q_collect, if_valid)
        if not valid:
            # Need to add an element here before the existing node
            for q in q_collect:
                new_slice_to_add[q] = single_qubit_gate[q]
                single_qubit_gate[q] = None
            for q in range(len(single_qubit_gate)):
                if self.node_map[q] is not None and self.node_map[q].to_delete():
                    self.node_map[q] = None
            # print("We get the following lists; ", single_qubit_gate, new_slice_to_add)
            return single_qubit_gate, new_slice_to_add
        for q in range(len(single_qubit_gate)):
            if self.node_map[q] is not None and self.node_map[q].to_delete():
                self.node_map[q] = None
        return single_qubit_gate, new_slice_to_add

    def reset_node_map(self, node, old_qubits):
        new_qubits = node.get_qubits()
        for n in list(set(old_qubits) - set(new_qubits)):
            self.node_map[n] = None
    
    def slice_control_simplify(self, slice_before, slice_after):
        def resolve_gate(gate, deactivate_qubits=[]):
            for q in list(set(gate.get_control())-set(deactivate_qubits)):
                states_prev, states_next = slice_before.get_wire_status(q), slice_after.get_wire_status(q)
                if (not states_prev.contains(gate.get_control_value_for_qubit(q))) or \
                    (not states_next.contains(gate.get_control_value_for_qubit(q))):
                    gate.set_matrix(np.eye(gate.get_matrix().shape[0]))
                elif states_prev.is_equal(gate.get_control_value_for_qubit(q)) or states_next.is_equal(gate.get_control_value_for_qubit(q)):
                    gate.delete_control(q)
            deactivate_qubits.append(gate.get_target())
            return deactivate_qubits
        
        for nodes in self.get_all_nodes():
            old_qubits = nodes.get_qubits()
            if nodes.get_node_type().is_gateset():
                nodes_deact = []
                for gates in nodes.getAllNodes():
                    nodes_deact = resolve_gate(gates, nodes_deact)            
            elif nodes.get_node_type().is_gate():
                resolve_gate(nodes)
            self.reset_node_map(nodes, old_qubits)

    def slice_phase_simplify(self, slice_before, slice_after):
        def resolve_gate(gate, deactivate_qubits=[]):
            q = gate.get_target()
            states_prev, states_next = slice_before.get_wire_status(q), slice_after.get_wire_status(q)
            if gate.is_phase_gate() and (states_prev.decimate_gate(gate.get_matrix()) or states_next.decimate_gate(gate.get_matrix())):
                gate.set_matrix(np.eye(self.dimension))
            if not gate.is_phase_gate():
                deactivate_qubits.append(gate.get_target())
            return deactivate_qubits

        for nodes in self.get_all_nodes():
            old_qubits = nodes.get_qubits()
            if nodes.get_node_type().is_gateset():
                deactivate_qubits = []
                for gates in nodes.getAllNodes():
                    deactivate_qubits = resolve_gate(gates, deactivate_qubits)
            elif nodes.get_node_type().is_gate():
                resolve_gate(nodes)
            # self.reset_node_map(node, old_qubits)

    def remove_gate_from_qubit(self, gate, idx):
        assert gate == self.node_map[idx]
        self.node_map[idx] = None

    def collapse_nodes(self, edges_prev, edges_nxt):
        processed = []
        for i in range(len(self.node_map)):
            if self.node_map[i] not in processed and self.node_map[i] is not None \
                and self.node_map[i].get_node_type().is_gate() and self.node_map[i].get_gate_type().is_phase_gate()\
                and len(self.node_map[i].get_control()) > 0:
                    if self.node_map[i].is_uniform_phase_gate(edges_prev={"status": edges_prev.get_wire_status(self.node_map[i].get_target())},
                                                              edges_nxt={"status":edges_nxt.get_wire_status(self.node_map[i].get_target())})[0]:
                        # print("We delete the matrix ", i, self.node_map[i].get_matrix())
                        node = self.node_map[i]
                        self.node_map[node.get_target()] = None
                        if len(node.get_qubits()) > 1:
                            node.delete_last_qubit(node.get_target())
                            for q in node.get_qubits():
                                self.node_map[q] = self.node_map[i]
                            
                    processed.append(self.node_map[i])
                    for x in range(len(self.node_map)):
                        if self.node_map[x] is not None:
                            if x not in self.node_map[x].get_qubits():
                                assert False
            if self.node_map[i] is not None and self.node_map[i].to_delete():
                for q in self.node_map[i].get_qubits():
                    self.node_map[i] = None
            if self.node_map[i] is not None and np.allclose(self.node_map[i].get_matrix(), np.eye(self.dimension)):
                assert False
    
    def is_common_phase_present(self, gateslice_handle, index):
        gate_keys = []
        for gates in self.node_map:
            if gates is not None and not gates.get_node_type().is_gateset() and gates.get_node_type().is_gate():
                modified_matrix, phase_diff = check_matrix_phase(gates.get_matrix())
                if not np.allclose(np.exp(1j*phase_diff), 1):
                    return True
        return False      

    def unroll_gate_common_phases(self, gateslice_handle, index):
        gate_keys = []
        for i, gates in enumerate(self.node_map):
            if gates is not None and not gates.get_node_type().is_gateset() and gates.get_node_type().is_gate() and gates.get_target() == i:
                modified_matrix, phase_diff = check_matrix_phase(gates.get_matrix())
                if np.allclose(phase_diff, 0):
                    gates.set_matrix(gateslice_handle.slice_edge_collection.node_works_on_qutrit(index, gates)[1].normalize_matrix(gates.get_matrix()))
                else: 
                    gates.set_matrix(modified_matrix)
                    gate_copy = gates.deepcopy()
                    gate_copy.set_matrix(np.eye(self.dimension, dtype=np.complex128)*np.exp(1j*phase_diff))
                    gate_copy.delete_last_qubit(gates.get_target())
                    gate_keys.append((gates, gate_copy))
        return gate_keys
    
    def remove_phases(self):
        for i, gates in enumerate(self.node_map):
            if gates is not None and gates.get_node_type().is_gate():
                gates.remove_phase()
