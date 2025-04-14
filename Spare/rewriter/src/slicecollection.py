import numpy as np
from rewriter.src.nodeset import LabelledNode
from rewriter.src.slice import GateSlice, control_config
from rewriter.src.sliceedges import EdgeCollection, EdgeSliceNodes
from rewriter.utils.graph_utils import if_mat_preserves_binary
from rewriter.utils.edge_types import EdgeType
from src.toffoli import *
import networkx as nx
import cirq
from rewriter.utils.graph_utils import *
import collections

def add_gate_with_interface(status_binary, slice):
    assert len(status_binary) == len(slice.get_qubits())
    dic = {}
    for g in slice.get_all_nodes():
        if g.get_node_type().is_qubit():
            pass
        else:
            for c in g.get_control():
                if g.get_control_value_for_qubit(c) == 2:
                    if status_binary[c] is None:
                        dic[g] = False
                elif g.get_control_value_for_qubit(c) == 1:
                    if status_binary[c] == "Anc":
                        dic[g] = False
            if g not in dic or dic[g]:
                if g.get_gate_type().turns_state_to_t():
                    status_binary[g.get_target()] = g
                elif g.get_gate_type().is_binary() and status_binary[g.get_target()] == "Anc":
                    status_binary[g.get_target()] = None
    for g in dic:
        if not dic[g]:
            slice.remove_node(g)
    return slice, status_binary


class GateSliceCollection:
    def __init__(self, clir=None, qubits=0, dict={}):
        self.slice_collection = []
        self.slice_edge_collection = EdgeSliceNodes(gate_slice_collection=self)
        if clir is None and qubits==0 and len(dict.keys()) == 0:
            assert False
        elif clir is None and len(dict.keys()) == 0:
            self.num_qubits = qubits
            return
        elif len(dict.keys()) != 0 and clir is None:
            self.slice_collection = dict["slice_collection"]
            self.slice_edge_collection = dict["slice_edge_collection"]
            self.num_qubits = dict["num_qubits"]
            self.slice_edge_collection.set_gate_slice_collection(self)
            self.dimension = dict["dimension"]
            return
        elif clir is not None:
            pass
        else:
            raise ModuleNotFoundError
        self.dimension=clir.dimension
        all_qubits = clir.get_all_qubits()
        for q in all_qubits:
            self.num_qubits = max(1, clir.node_info[q].get_target() + 1)
        self.slice_edge_collection = EdgeSliceNodes(gate_slice_collection=self)
        startgateslice = GateSlice(qubits=self.num_qubits, gate=None, dimension=clir.dimension)
        endgateslice = GateSlice(qubits=self.num_qubits, gate=None, dimension=clir.dimension)
        for nodes in nx.topological_sort(clir.dag):
            if nodes in clir.node_info:
                if clir.node_info[nodes].get_node_type().is_qubit() and clir.node_info[nodes].is_root():
                    if startgateslice.get_node(clir.node_info[nodes].get_qubits()[-1]) is None:
                        startgateslice.set_node(clir.node_info[nodes])
        self.add_slice(startgateslice)
        for nodes in nx.topological_sort(clir.dag):
            if nodes in clir.node_info:
                if clir.node_info[nodes].get_node_type().is_qubit() and clir.node_info[nodes].is_sink():
                    if startgateslice.get_node(clir.node_info[nodes].get_qubits()[-1]) != clir.node_info[nodes]:
                        endgateslice.set_node(clir.node_info[nodes])
                elif clir.node_info[nodes].get_node_type().is_gate():
                    gateslice = GateSlice(qubits=self.num_qubits, gate=clir.node_info[nodes], dimension=clir.dimension)
                    self.add_slice(gateslice)
        self.add_slice(endgateslice, no_edges=True)
        assert startgateslice.if_qubit_annotation() and endgateslice.if_qubit_annotation()
        assert self.slice_collection[0].if_qubit_annotation() and self.slice_collection[-1].if_qubit_annotation()
        assert self == self.slice_edge_collection.gate_slice_collection
        self.slice_edge_collection.update_all_mappings()
        return
    
    def reset(self, qubits, dimension=None):
        # Resets the datastructure
        curr_copy = self.slice_edge_collection
        self.slice_collection = []
        self.slice_edge_collection = EdgeSliceNodes(gate_slice_collection=self)
        self.num_qubits = qubits
        if dimension is not None:
            self.dimension = dimension
        return curr_copy
    
    def update_qubits(self):
        # Used with resets should be subsumed together
        self.num_qubits = len(self.slice_collection[0].node_map)
    
    def remove_high_control_gates(self):
        for idx, sl in enumerate(self.slice_collection):
            if idx != 0 and idx != len(self.slice_collection)-1:
                # print("index number ", idx, end="")
                sl.slice_control_simplify(self.slice_edge_collection.get_slice_edge(idx-1), self.slice_edge_collection.get_slice_edge(idx))

    def remove_high_phase_gates(self):
        for idx, sl in enumerate(self.slice_collection):
            if idx != 0 and idx != len(self.slice_collection)-1:
                sl.slice_phase_simplify(self.slice_edge_collection.get_slice_edge(idx-1), self.slice_edge_collection.get_slice_edge(idx))            

    def return_iterator(self):
        # Handler to get size of object
        return len(self.slice_collection), len(self.slice_edge_collection.list_of_edges)
    
    def print(self):
        # Visualizer
        for i, s in enumerate(self.slice_collection):
            print(s, end="")
            print(i, len(self.slice_collection), len(self.slice_edge_collection.list_of_edges))
            if i != len(self.slice_collection) - 1:
                print(self.slice_edge_collection.list_of_edges[i])
        print("__________________________________________")

    def merge_greedily(self, dont_create_mappings=False, binary=False):
        '''Merge and fuses all slices of the circuit for the circuit
            TODO: Modify the merging to be able to handle a sequence of nodes etc'''
        old_len = len(self.slice_collection)
        new_len = old_len + 1
        while new_len != old_len:
            slices_to_delete = []
            s_idx = []
            old_len = new_len
            for sl in range(1, len(self.slice_collection)):
                if self.slice_collection[sl].if_slices_mergeable(self.slice_collection[sl - 1]):
                    assert not self.slice_collection[sl].if_qubit_annotation() and not self.slice_collection[sl - 1].if_qubit_annotation() 
                    self.slice_collection[sl].merge(self.slice_collection[sl-1], binary=binary)
                    self.slice_edge_collection.merge_index(sl, sl-1, dont_create_mappings=dont_create_mappings)
                    slices_to_delete.append(self.slice_collection[sl-1])
                    s_idx.append(sl - 1)
            for s in slices_to_delete:
                self.slice_collection.remove(s)
            self.slice_edge_collection.prune_edgeslices()
            new_len = len(self.slice_collection)
        return

    def collapse_states(self):
        for i in range(1, len(self.slice_collection) - 1):
            slice = self.slice_collection[i]
            slice.collapse_nodes(self.slice_edge_collection.get_edge_slice(i - 1), self.slice_edge_collection.get_edge_slice(i))
    
    def prune_null_slices(self):
        slices = []
        for i in range(1, len(self.slice_collection)-1):
            if self.slice_collection[i].is_null_slice():
                # print("Null ", end="")
                slices.append(self.slice_collection[i])
        for s in slices:
            assert s in self.slice_collection
            idx = self.slice_collection.index(s)
            self.slice_collection.remove(s)
            self.slice_edge_collection.merge_index(idx, idx-1)
            self.slice_edge_collection.prune_edgeslices()
        return
    
    def get_circuit_width(self):
        '''
            Returns the number of qubits in the circuit including the ancilas
        '''
        circuit_with_tracker = [0]*self.num_qubits  
        for i in range(1, len(self.slice_collection)-1):
            for q in range(self.num_qubits):
                if self.slice_collection[i].get_node(q) is not None:
                    circuit_with_tracker[q] = 1
        return sum(circuit_with_tracker)
        
    def control_modify(self, gate, matrix, qubit):
        '''Push a gate through the gate control'''
        node = LabelledNode()
        node.set_qubits([qubit])
        node.set_matrix(matrix)
        assert node.get_node_type().is_gate()
        if node.get_gate_type().state_permutation_gates():
            gate.control_update_qubit(node, qbit=qubit)
        else:
            # Cant just modify control value need to add gates
            raise NotImplementedError
    
    def remove_anc_cnots(self):
        previous_maps = list(range(self.slice_collection[0].qubits))
        controlvalues = [None]*len(self.slice_collection[0].qubits)
        for i in range(1, len(self.slice_collection)-1):
            result, control_values, previous_maps = self.slice_collection[i].reroll_cnots(previous_maps, control_values=control_values,
                                                                                          previous_states=self.slice_edge_collection.get_edge_slice(i-1))
            if not result:
                for j in range(1, i):
                    self.slice_collection[i].reset_ancilla_maps()
                return
        for i in range(1, len(self.slice_collection)-1):
            self.slice_collection[i].parse_and_roll_cnots()
        self.self.prune_null_slices()

    def control_rewrite(self, debug=False):
        '''Push all the 1 qubit gates through ther circuit'''
        def smarter_status_bin(first_slice):
            default_status_bin = [None]*self.num_qubits
            for i, q in enumerate(first_slice.node_map):
                if not q.get_node_type().is_qubit():
                    return [None]*self.num_qubits
                elif q._is_ancilla:
                    default_status_bin[i] = "Anc"
            return default_status_bin
            
        single_gates = [None] * self.num_qubits
        status_binary = smarter_status_bin(self.slice_collection[0])
        len_slice_coll = len(self.slice_collection) - 1
        new_slice_collection = GateSliceCollection(None, self.num_qubits)
        i = 0
        if debug:
            print(self.build_circuit()[0])
        while i < len_slice_coll:
            single_gates, new_slice_checked = self.slice_collection[i].update_gate_slice(single_gates)
            if new_slice_checked != [None] * self.num_qubits:
                # We need to modify this here by adding the new_slice we get here
                slice, status_binary = add_gate_with_interface(status_binary, GateSlice(qubits=self.num_qubits, gatelist=new_slice_checked, dimension=self.dimension))
                if not slice.is_null_slice():
                    new_slice_collection.add_slice(slice)
            if not self.slice_collection[i].is_null_slice():
                slice, status_binary = add_gate_with_interface(status_binary, self.slice_collection[i])
                if not slice.is_null_slice():
                    new_slice_collection.add_slice(slice) 
            i += 1
        new_slice = GateSlice(qubits=self.num_qubits, gatelist=single_gates, dimension=self.dimension)
        if not new_slice.is_null_slice():
            new_slice_collection.add_slice(new_slice)
        new_slice_collection.add_slice(self.slice_collection[-1], no_edges=True)
        self.slice_collection = new_slice_collection.slice_collection
        self.slice_edge_collection = new_slice_collection.slice_edge_collection
        self.slice_edge_collection.set_gate_slice_collection(self)
        assert len(self.slice_collection) - 1 == len(self.slice_edge_collection.list_of_edges)
        assert self == self.slice_edge_collection.gate_slice_collection
        return

    def add_slice(self, sl, no_edges=False):
        assert self == self.slice_edge_collection.gate_slice_collection
        assert type(sl) == GateSlice
        if self.num_qubits != -1 and sl.get_num_qubits() != self.num_qubits:
            raise Exception
        self.slice_collection.append(sl)
        if not no_edges:
            self.slice_edge_collection.add_slice(sl, len(self.slice_collection) - 1)

    def add_a_gate_after_gate(self, end_copy, end):
        _, idx = self.find_gate(end)
        assert self == self.slice_edge_collection.gate_slice_collection
        gateslice = GateSlice(qubits=self.num_qubits, gate=end_copy, dimension=self.dimension)
        self.slice_collection = [*self.slice_collection[:idx+1], gateslice, *self.slice_collection[idx+1:]]
        self.slice_edge_collection.add_slice_after_idx(gateslice, idx)
        self.slice_edge_collection.update_all_mappings()
        
    def add_gate_list_after_node(self, gates_list, node, position):
        _, idx = self.find_gate(node)
        assert self == self.slice_edge_collection.gate_slice_collection
        if node.get_node_type().is_gateset():
            node.addActiveGate(gates_list, position)
            node.delActiveGate(position)
        else:
            for g in reversed(gates_list):
                self.add_a_gate_after_gate(g, node)
            self.delete_gate(node)
        
    def get_slice(self, index):
        assert isinstance(index, int)
        return self.slice_collection[index]
    
    def delete_slice(self, slice):
        self.slice_collection.remove(slice)
        
    def pop_slice(self, slice_idx):
        self.slice_collection.pop(slice_idx)
        self.slice_edge_collection.delete_slice(self.slice_collection[slice_idx]) 

    def find_gate(self, node):
        if node is None:
            raise Exception()
        qubit = node.get_qubits()[0]
        for s_idx in range(len(self.slice_collection)):
            slices = self.slice_collection[s_idx]
            if node == slices.get_node(qubit):
                return slices, s_idx
        raise Exception("Node not found")
    
    def parse_for_matching_gates(self, gate, i, forward_direction=True, swap_through_gates=False,
                                 invalid_swap=False, desired_result=None, return_middle=False, debug=False):
        if not forward_direction:
            start, end, diff = (i-1, 0, -1)
        else:
            start, end, diff = (i+1, len(self.slice_collection), 1)
        middle_gates = []
        for j in range(start, end, diff):
          for g2 in self.slice_collection[j].get_all_nodes():
            if g2.get_node_type().is_gate():
              if debug:
                  print("current slices ", i, j)
                  print(self.slice_collection[i])
                  print(self.slice_collection[j])
                  
              result = gate.speculatively_swap_with_gate(g2, forward_dir=forward_direction, swap_through_gates=swap_through_gates,
                                                         invalid_swap=invalid_swap, desiredOp=desired_result, debug=debug)
              if debug:
                  print(g2, gate)
                  print(result)
                  print("------------------")
              if not result[0] or len(result[1][1]) > 0 or len(result[1][0]) > 0:
                if (len(result[1][1]) > 0 or len(result[1][0]) > 0) and return_middle and (desired_result is None or result[1][0] == desired_result):
                  return result, [gate, g2], (i, j), middle_gates
                if (len(result[1][1]) > 0 or len(result[1][0]) > 0) and (desired_result is None or result[1][0] == desired_result):
                  return result, [gate, g2], (i, j)
                elif not result[0] and return_middle:
                  return result, [gate, g2], (i, j), middle_gates
                elif not result[0]:
                  return result, [gate, g2], (i, j)
                elif (desired_result is not None and result[1][0] != desired_result):
                  pass
                else:
                  print(result, desired_result, len(result[1][1]) > 0 or len(result[1][0]) > 0, (desired_result is None or result[1][0] == desired_result))
                  raise Exception
              elif result[0] and return_middle and not g2.get_node_type().is_gateset() and (gate.get_target() in g2.get_control() or g2.get_target() in gate.get_control()):
                middle_gates.append(([{"ctr": x, "val": gate.get_control_value_for_qubit(x)} for x in gate.get_control()], g2))
        if return_middle:
          return [False], None, None, None
        return [False], None, None
        
    def delete_gate(self, node):
        assert node.get_node_type().is_gate()
        slice, s_idx = self.find_gate(node)
        slice.remove_node(node)
        del node
    
    def delete_controls_on_gate(self, gate, controls):
        assert len(controls) > 0
        slice, s_idx = self.find_gate(gate)
        for c in controls:
            gate.delete_qubit(c, control_only=True)
            slice.remove_gate_from_qubit(gate, c)
    
    def set_control_values(self):
        starting_index = 1
        while starting_index < len(self.slice_collection) - 1:
            idx_and_gates = self.slice_collection[starting_index].get_single_qubit_gates()
            new_index = starting_index + 1
            if len(idx_and_gates) > 0:
                while new_index < len(self.slice_collection):
                    idx_and_gates = self.slice_collection[new_index].update(idx_and_gates)
        for i, s in enumerate(self.slice_collection):
            s.prune_gates()
        self.slice_edge_collection.update_all_mappings()
    
    def remove_phases(self):
        for i, slice in enumerate(self.slice_collection):
            slice.remove_phases()
        
    @classmethod
    def create_copy(cls, slicecollection, sliceEdges, numQubits, dimension):
        '''Save and prepare a deepcopy of a slicecollection object'''
        slicecollection = [x.create_copy() for x in slicecollection]
        sliceEdges = EdgeCollection.create_copy()
        data = {"slice_collection":slicecollection, "slice_edge_collection":sliceEdges, "num_qubits":numQubits, "dimension":dimension}
        return cls(data=data)
            
    def build_circuit(self, if_final_circ=False, check_circ=False, breakdown="vvv", dimension=3):
        # Builds the circuit for visualization
        # self.dimension=2
        dimension = self.dimension
        qubits = cirq.LineQid.range(self.num_qubits, dimension=dimension)
        circuit = cirq.Circuit()
        if check_circ:
            for q in qubits:
                circuit.append(cirq.IdentityGate(num_qubits=1).on(q))
        for idx, s in enumerate(self.slice_collection):
            gate_arr = []
            for g in s.node_map:
                if g is not None and g not in gate_arr and g.get_node_type().is_gate():
                    control_qubits = []
                    if if_final_circ:
                        for c in self.get_control():
                            control_qubits.append(qubits[c])
                            circuit.append(control_config(qubits, self, c, dimension=dimension))
                        circuit.append(self._build_controlled_gate(self, idx, qubits, dimension=dimension))
                        for c in g.get_control():
                            circuit.append(control_config(qubits, g, c, dimension=dimension))                    
                    else:
                        circuit.append(g.realizeNode(qubits, breakdown=breakdown))
                    gate_arr.append(g)
        return circuit, qubits
        
    def _build_controlled_gate(self, g, s, qubits, dimension=3):
        # Unrolls gates for purpose of visulaization
        # Redundant with use of gate unrolling pass in ancilla convertor
        # We need information about the gate status to build this gate in reality 
        # Cannot simply delete 
        g_copy_1 = g.deepcopy()
        g_copy_2 = g.deepcopy()
        g_copy_1.set_matrix(cirq.unitary(QutritX12Gate()))
        g_copy_2.set_matrix(cirq.unitary(QutritX12Gate()))
        circuit = cirq.Circuit()
        targets = g.get_control()
        targets = [x for i, x in enumerate(targets) if i%2 == 1]
        if len(g.get_control()) > 1: 
            for c in targets:  # g.get_control()
                if self.slice_edge_collection.get_edge_slice(s-1).get_wire_status(c) == EdgeType.Basis01 or\
                  self.slice_edge_collection.get_edge_slice(s).get_wire_status(c) == EdgeType.Basis01:
                    assert g.get_control_value_for_qubit(c) != 2
                    assert self.slice_edge_collection.get_edge_slice(s-1).get_wire_status(c) ==\
                          self.slice_edge_collection.get_edge_slice(s).get_wire_status(c)
                    g_copy_1.delete_qubit(c)
                    g_copy_2.delete_qubit(c)
                    g_copy_1.set_target(c)
                    g_copy_2.set_target(c)
                    circuit.append(self.build_controlled_gate(g_copy_1, s, qubits))
                    circuit.append(control_config(qubits, None, c))
                    if if_mat_preserves_binary(g.get_matrix() @ cirq.unitary(QutritX01Gate()) @ cirq.unitary(QutritX12Gate())) or\
                        if_mat_preserves_binary(g.get_matrix() @ cirq.unitary(QutritX12Gate()) @ cirq.unitary(QutritX01Gate())):
                        circuit.append(cirq.MatrixGate(cirq.unitary(QutritX12Gate()) @ g.get_matrix(),
                                                       qid_shape=[dimension]).on(qubits[g.get_target()]).controlled_by(qubits[c]))
                        circuit.append(cirq.MatrixGate(cirq.unitary(QutritX12Gate()),
                                                       qid_shape=[dimension]).on(qubits[g.get_target()]).controlled_by(qubits[c]))
                    else:
                        circuit.append(cirq.MatrixGate(g.get_matrix(), qid_shape=[dimension]).on(qubits[g.get_target()]).controlled_by(qubits[c]))
                    circuit.append(control_config(qubits, None, c))
                    circuit.append(self.build_controlled_gate(g_copy_2, s, qubits))
                    return circuit
        control_qubits = []
        for c in g.get_control():
            control_qubits.append(qubits[c])
        circuit.append(cirq.MatrixGate(g.get_matrix(), qid_shape=[dimension]).on(qubits[g.get_target()]).controlled_by(* control_qubits))
        return circuit
    
    def find_similar_gates(self, node):
        s, sidx = self.find_gate(node)
        qubit = node.get_qubits()[0]
        result = {}
        for s_idx in range(len(self.slice_collection[sidx+1:])):
            slices = self.slice_collection[s_idx]
            if node == slices.get_node(qubit):
                result = {}
            elif slices.get_node(qubit) is not None and \
                tuple(slices.get_node(qubit).get_qubits()) == tuple(node.get_qubits()) and \
                tuple(slices.get_node(qubit).get_control_values()) == tuple(node.get_control_values()):
                if if_mat_preserves_binary(slices[qubit].get_matrix() @ node.get_matrix()):      
                    result[s_idx] = slices.get_node(qubit)
                if if_mat_preserves_binary(slices[qubit].get_matrix() @ node.get_matrix() @ cirq.unitary(QutritX01Gate())):      
                    result[s_idx] = slices.get_node(qubit)
        return result
    
    def find_inv_gates(self, node, set=None, end_limit=None, exact=False, just_binary=False, check_idx=True):
        def node_cases(orig_node, curr_node, dummy_node, index, result, reverse=False):
            if curr_node is None or curr_node.get_node_type().is_qubit():
                pass
            elif curr_node == orig_node:
                return False, result, dummy_node.get_matrix()
            elif curr_node.get_node_type().is_gateset():
                True, result, dummy_node.get_matrix()
            elif curr_node.is_control_subset(dummy_node) and \
                curr_node.get_target() in dummy_node.get_control():
                if curr_node.get_gate_type().state_permutation_gates():
                    if reverse:
                        dummy_node.control_update_qubit(curr_node, curr_node.get_target())
                    else:    
                        copy_ = curr_node.deepcopy()
                        copy_.set_matrix(copy_.get_matrix().conjugate().T)
                        dummy_node.control_update_qubit(copy_, curr_node.get_target())
                else:
                    return False, result, dummy_node.get_matrix()
            elif curr_node.is_control_subset(dummy_node) and curr_node.get_target() == dummy_node.get_target():
                if reverse:
                    dummy_node.set_matrix(dummy_node.get_matrix() @ curr_node.get_matrix())    
                else:
                    dummy_node.set_matrix(curr_node.get_matrix() @ dummy_node.get_matrix())
                # Wont be an inverse unless.. 
                if curr_node.get_qubits() == dummy_node.get_qubits():
                    if exact and dummy_node.get_gate_type().is_phase_gate() and np.allclose(dummy_node.get_matrix(), np.eye(self.dimension)):
                        result[index] = curr_node
                        return True, result, dummy_node.get_matrix()
                    elif not exact and dummy_node.get_gate_type().is_binary():
                        result[index] = curr_node
                        return True, result, dummy_node.get_matrix()
            elif curr_node.is_control_mismatch(dummy_node):
                pass
            elif not curr_node.is_control_subset(dummy_node) and curr_node.get_target() == dummy_node.get_target():
                if dummy_node.is_control_subset(curr_node):
                    dummy_copy = curr_node.deepcopy()
                    if reverse:
                        dummy_copy.set_matrix(dummy_node.get_matrix() @ curr_node.get_matrix() @ dummy_node.get_matrix().conjugate().T)
                    else:
                        dummy_copy.set_matrix(dummy_node.get_matrix().conjugate().T @ curr_node.get_matrix() @ dummy_node.get_matrix())
                    if not dummy_copy.get_gate_type().is_binary():
                        return False, result, dummy_node.get_matrix()
                else:
                    return False, result, dummy_node.get_matrix()
            elif not curr_node.is_control_subset(dummy_node) and curr_node.get_target() in dummy_node.get_control()\
                  and not curr_node.is_immutable(dummy_node.get_control_value_for_qubit(curr_node.get_target())):
                return False, result, dummy_node.get_matrix()
            return True, result, dummy_node.get_matrix()
        
        assert node is not None
        assert node.get_node_type().is_gate()
        s, sidx = self.find_gate(node)
        result = {}
        dummy_node = node.deepcopy()
        if set is not None:
            n, idx_ = self.find_gate(node)
            if len(set) == 0:
                return result
            sidx_ = sidx
            if self.find_gate(set[0])[1] > sidx_:
                reverse = False
            else:
                reverse = True
                set = reversed(set)
            for s in set:
                s_, s_idx_ = self.find_gate(s)
                res, result, mat = node_cases(node, s, dummy_node, s_idx_, result, reverse=reverse)
                assert np.allclose(mat, dummy_node.get_matrix())
                if not res:
                    return result
        else:
            if end_limit is None:
                end_limit = len(self.slice_collection)
            assert end_limit > sidx
            for s_idx in range(sidx+1, end_limit):
                slices = self.slice_collection[s_idx]
                for n in slices.get_all_nodes():
                    res, result, mat = node_cases(node, n, dummy_node, s_idx, result)
                    assert np.allclose(mat, dummy_node.get_matrix())
                    if not res:
                        return result
        return result

    def get_nodes_on_qubit(self, q, last_slice=-1):
        assert isinstance(q, int)
        nodes = []
        for slice in range(len(self.slice_collection)):
            if last_slice > 0 and slice == last_slice:
                break
            slice = self.slice_collection[slice]
            node = slice.get_node(q)
            if node is not None:
                nodes.append(node)
        return nodes

    def get_nodes_before_after_gate(self, gate):
        slice, sidx = self.find_gate(gate)
        result_before = []
        for s_ in range(len(self.slice_collection[:sidx])):
            for g in self.slice_collection[s_].get_all_nodes():
                if g is not None and g not in result_before and g.get_node_type().is_gate():
                    result_before.append(g)        
        result_after = []
        for s_ in range(sidx + 1, len(self.slice_collection)):
            for g in self.slice_collection[s_].get_all_nodes():
                if g is not None and g not in result_after and g.get_node_type().is_gate():
                    result_after.append(g)
        return result_before, result_after

    def get_all_nodes_between(self, node1, node2, qubits=[]):
        slice1, s_idx1 = self.find_gate(node1)
        slice2, s_idx2 = self.find_gate(node2)
        if s_idx1 > s_idx2:
            node_ = node2
            node2 = node1
            node1 = node_
            s_idx_ = s_idx2
            s_idx2 = s_idx1
            s_idx1 = s_idx_
        assert s_idx2 > s_idx1
        result = []
        nodes_processed = []
        for slices in self.slice_collection[s_idx1+1: s_idx2]:
            for q in list(set(node1.get_qubits()).union(set(node2.get_qubits()))):
                if slices.get_node(q) is not None and slices.get_node(q) not in nodes_processed:
                    node = slices.get_node(q)
                    if node.get_node_type().is_gateset():
                        for g in node.getActiveSet():
                            if g.get_target() in qubits:
                                result.append(g)
                    elif node.get_node_type().is_gate():
                        result.append(node)
                    nodes_processed.append(slices.get_node(q))
        return result
    
    def get_if_nodes_mutable(self, node1, node2, intersect, just_binary=False, check_idx=True, set_={}):
        if not check_idx and self.find_gate(node1)[1] > self.find_gate(node2)[1]:
            # If node 2 comes later need to swap them to check the order is correct
            # COmpare indexes
            tmp = node1
            node1 = node2
            node2 = tmp
        assert self.find_gate(node2)[1] > self.find_gate(node1)[1]
        running_node = node1.deepcopy()
        # Need to remove pairs of gates
        remove_arr = []
        for i in intersect:
            if i not in remove_arr:
                end_limit = self.find_gate(node2)
                result = self.find_inv_gates(i, end_limit=end_limit[1], exact=True)
                for r in result:
                    assert self.get_if_nodes_mutable(i, result[r], self.get_all_nodes_between(i, result[r]))
                    remove_arr.append(i)
                    remove_arr.append(result[r])
                    break
        for x in remove_arr:
            if x in intersect:
                intersect.remove(x)
        for i in intersect:
            assert self.find_gate(i)[1] <= self.find_gate(node2)[1]
            if i is None or i.get_node_type().is_qubit():
                pass
            elif i == node2:
                return node2.is_synchronized(running_node)
            elif i.is_control_mismatch(running_node):
                pass
            elif i.is_control_subset(running_node) and \
                i.get_target() in running_node.get_control():
                if i.get_gate_type().state_permutation_gates():
                    running_node.control_update_qubit(i, i.get_target())
                else:
                    return False
            elif not i.is_control_subset(running_node) and i.get_target() in running_node.get_control() and \
                not self.if_state_mutable(i, running_node.get_control_value_for_qubit(i.get_target())):
                return False
        return True

    def remove_mutable_gates(self, bookends, list_arr):
        rem_arr = []
        for b in list_arr:
            if b not in rem_arr and b.get_target() in bookends[0].get_control():
                if self.if_state_mutable(b, bookends[0].get_control_value_for_qubit(b.get_target())):
                    rem_arr.append(b)
        for x in rem_arr:
            list_arr.remove(x)
        return list_arr
    
    def if_state_mutable(self, gate, control_value):
        if gate.is_phase_gate():
            return True
        elif gate.get_gate_type().on_01_basis() and control_value == 2:
            return True
        elif gate.get_gate_type().on_12_basis() and control_value == 0:
            return True
        elif gate.get_gate_type().on_02_basis() and control_value == 1:
            return True
        else:
            return False

    def remove_inv_pairs(self, gateset, qubits=None, phase_ok=False, just_binary=False):
        deleting_nodes = []
        st = -1
        for x in gateset:
            y = self.find_gate(x)[1] 
            assert y >= st
            st = y
        del_arr = []
        for node in gateset:
            del_nodes = []
            if qubits is not None and node is not None and node.get_target() not in qubits:
                pass
            elif node is not None and node not in deleting_nodes and node.get_node_type().is_gate():
                gateset_ = []
                for x in gateset:
                    if self.find_gate(x)[1] > self.find_gate(node)[1]:
                        gateset_.append(x)
                res = self.find_inv_gates(node, set=gateset_,
                                          exact=not phase_ok, just_binary=just_binary)
                early_exit = False
                for r in res:
                    assert self.find_gate(node)[1] < self.find_gate(res[r])[1]
                    if res[r] not in gateset:
                        assert False
                        continue
                    dummy_node = node.deepcopy()
                    for gate in self.get_all_nodes_between(node, res[r]):
                        if gate in gateset:
                            if gate.is_control_mismatch(dummy_node):
                                pass
                            elif gate.is_control_subset(dummy_node) and gate.get_target() == dummy_node.get_target() and\
                                  gate.get_qubits() == dummy_node.get_qubits():
                                dummy_node.set_matrix(gate.get_matrix() @ dummy_node.get_matrix())
                                del_arr.append(gate)
                            elif gate.is_control_subset(dummy_node) and gate.get_target() == dummy_node.get_target() and\
                                  gate.get_qubits() != dummy_node.get_qubits():
                                dummy_node.set_matrix(gate.get_matrix().conjugate().T @ dummy_node.get_matrix() @ gate.get_matrix())
                                del_arr.append(gate)
                            elif gate.is_control_subset(dummy_node) and gate.get_target() in dummy_node.get_control():
                                if gate.get_gate_type().state_permutation_gates():
                                    dummy_node.control_update_qubit(gate, gate.get_target())
                                else:
                                    early_exit = True
                        else:
                            pass
                        if early_exit:
                            del_arr = []
                            break
                    if not early_exit:
                        del_arr.append(node)
                        del_arr.append(res[r])
                        for x in del_nodes:
                            if x in gateset:
                                gateset.remove(x)
                        break
        return gateset

    def get_inv_pairs(self, gateset, gate, qubits=None, phase_ok=False, just_binary=False, mode="normal", ogs=None):
        deleting_nodes = []        
        if len(gateset) == 1:
            return deleting_nodes
        if qubits is not None and gate is not None and gate.get_target() not in qubits:
            pass
        elif gate is not None and gate not in deleting_nodes and gate.get_node_type().is_gate():
            gateset_ = gateset.copy()
            gateset_.remove(gate)
            res = self.find_inv_gates(gate, set=gateset_,
                                      exact=not phase_ok, just_binary=just_binary, check_idx=False)
            for r in res:
                if self.get_if_nodes_mutable(gate, res[r],
                                             self.get_all_nodes_between(gate, res[r],
                                                                        qubits=gate.get_control()),
                                             just_binary=just_binary, check_idx=False, set_=gateset_):
                    if res[r] not in deleting_nodes:
                      # Delete here from the set and also add to deleting nodes
                      if self.find_gate(gate)[1] < self.find_gate(res[r])[1]:
                          dummy_node = gate.deepcopy()
                          start = gate
                          end = res[r]
                          assert mode == "reverse"
                      else:
                          dummy_node = res[r].deepcopy()
                          start = res[r]
                          end = gate
                          assert mode == "forward"                        
                      for g in self.get_all_nodes_between(start, end):
                          if g.is_control_mismatch(dummy_node):
                            pass
                          elif dummy_node.is_control_subset(g) and g.get_target() == dummy_node.get_target():
                            if g in gateset:
                                gateset.remove(g)
                                deleting_nodes.append(g)
                            else:
                                pass
                          elif g.is_control_subset(dummy_node) and g.get_target() in dummy_node.get_control():
                              if g.get_gate_type().state_permutation_gates():
                                  dummy_node.control_update_qubit(g, g.get_target())
                              else:
                                  print(g)
                                  assert False
                      gateset.remove(res[r])
                      gateset.remove(gate)
                      if mode == "reverse":
                        deleting_nodes.insert(0, gate)
                        deleting_nodes.append(res[r])
                      else:
                        deleting_nodes.insert(0, res[r])
                        deleting_nodes.append(gate)
                      mat = np.eye(self.dimension, dtype=np.complex128)
                      for g in deleting_nodes:
                        mat = g.get_matrix() @ mat
                      if just_binary and self.dimension == 3:
                        assert if_mat_preserves_binary(mat) or if_mat_preserves_binary(mat @ cirq.unitary(QutritX01Gate())) or is_binary_mapping(mat)
                      elif self.dimension == 3:
                        assert if_mat_preserves_binary(mat) or if_mat_preserves_binary(mat @ cirq.unitary(QutritX01Gate()))
                      return deleting_nodes
                else:
                    pass
        else:
            assert False
        return deleting_nodes