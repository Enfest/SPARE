import cirq
from ancilla_convertor.gate_convertor import gate_convertor
from rewriter.src.slice import GateSlice
from rewriter.src.slicecollection import GateSliceCollection
import copy
import numpy as np
from rewriter.src.nodeset import LabelledNode
from rewriter.utils.edge_types import EdgeType
from src.toffoli import *
from rewriter.utils.graph_utils import check_matrix_phase

def generate_ancilla_maps(gateslices):
    '''
        We use a less efficient circuit mapper implementation maybe we should use a stack for the next version
    '''
    slice_max, edge_max = gateslices.return_iterator()
    gateslices.slice_edge_collection.update_all_mappings(debug=False)
    # Cant be first two slices
    ancilla_maps = [0] * len(gateslices.get_slice(0).get_qubits())
    currently_set = [False] * len(gateslices.get_slice(0).get_qubits())
    next_ancilla = len(ancilla_maps)
    for i in range(1, slice_max-1):
        pos = []
        ctrl2s = []
        
        for q in gateslices.get_slice(i).get_qubits():
            if gateslices.slice_edge_collection.get_edge_slice(i).get_wire_status(q) == EdgeType.Basis012:
                pos.append(q)
        for q in gateslices.get_slice(i).get_qubits():
            if gateslices.get_slice(i).node_map[q] is not None and len(gateslices.get_slice(i).node_map[q].get_qubits()) > 2:
                if 2 in gateslices.get_slice(i).node_map[q].get_control_values() and gateslices.get_slice(i).node_map[q] not in ctrl2s:
                    ctrl2s.append(gateslices.get_slice(i).node_map[q])
        
        for q in gateslices.get_slice(i).get_qubits():
            if ancilla_maps[q] != 0 and currently_set[q] and gateslices.slice_edge_collection.get_edge_slice(i).get_wire_status(q) == EdgeType.Basis01:
                next_ancilla -= 1
                currently_set[q] = False
        
        for q in gateslices.get_slice(i).get_qubits():
            if ancilla_maps[q] == 0 and gateslices.slice_edge_collection.get_edge_slice(i).get_wire_status(q) == EdgeType.Basis012:
                ancilla_maps[q] = next_ancilla
                next_ancilla += 1
                currently_set[q] = True
            elif ancilla_maps[q] != 0 and not currently_set[q] and gateslices.slice_edge_collection.get_edge_slice(i).get_wire_status(q) == EdgeType.Basis012:
                if next_ancilla < ancilla_maps[q]:
                    currently_set[q] = True                
                else:    
                    ancilla_maps[q] = next_ancilla
                    next_ancilla += 1
                    currently_set[q] = True                
                    
        ancilla_pos = [(i, x) for i, x in enumerate(ancilla_maps) if i in pos]
        if len(pos) > 1:
            for i, q_a_tup in enumerate(ancilla_pos):
                if i > 0:
                    if q_a_tup[1] in [x[1] for x in ancilla_pos[:i]]:
                        ancilla_maps[q_a_tup[0]] = next_ancilla
                        ancilla_pos[i] = (q_a_tup[0], next_ancilla)
                        next_ancilla += 1
    return ancilla_maps

def slice_convertor(slice, ancilla_maps, edge_data, slice_idx, ancilla_data):
    def check_gate(x):
        for i, q in enumerate(x.get_qubits()[:-1]):
            if q in x.get_qubits()[i+1:]:
                print("gate has multiple same qubit ", x)
                print(gates)
                assert False    
    new_slices = []
    new_qubits = max(0, max(ancilla_maps) - len(slice.node_map) + 1)
    new_node_map = copy.deepcopy(slice.node_map)
    if slice.is_root_slice() or slice.is_sink_slice():
        for x in range(new_qubits):
            new_node = LabelledNode(data={"qubits":[x+len(new_node_map)], "dimension":2, "ancilla":True, "matrix":None})
            new_node.set_root() if slice.is_root_slice() else new_node.set_sink()
            new_node_map.append(new_node)
    else:
        new_node_map.extend([None]*max(0, new_qubits))
    assert isinstance(new_node_map, list)
    for i, gates in enumerate(slice.node_map):
        if gates is not None and gates.get_node_type().is_gate():
            gatetype_info = edge_data.node_works_on_qutrit(slice_idx, gates)[1]
            check_matrix_phase(gatetype_info.normalize_matrix(gates.get_matrix()))
        if gates is not None and gates.get_node_type().is_gate() and not gates.get_gate_type().is_binary() and gates.get_target() == i:
            for q in gates.get_qubits():
                new_node_map[q] = None
            new_gates = gate_convertor(gates, ancilla_maps, gatetype_info, ancilla_data)
            for x in new_gates:
                check_gate(x)
            if len(new_gates) > 1:
                for j, new_gate in enumerate(new_gates[1:]):
                    new_slices.append(GateSlice(qubits=len(new_node_map), gate=new_gate, dimension=2, print_det=True))
            for q in new_gates[0].get_qubits():
                new_node_map[q] = new_gates[0]
        elif gates is not None and gates.get_node_type().is_gate() and gates.get_gate_type().is_binary() and gates.get_target() == i:
            if new_node_map[i].get_matrix().shape[0] == 2:
                temp = np.eye(3, dtype=np.complex128)
                temp[:2, :2] = new_node_map[i].matrix
                new_node_map[i].matrix = temp
            assert np.allclose(new_node_map[i].get_matrix()[2, 2] * new_node_map[i].get_matrix()[2, 2].conjugate(), 1)
            for q in gates.get_qubits():
                new_node_map[q] = None
            new_gate = gate_convertor(gates, ancilla_maps, gatetype_info, ancilla_data)
            check_gate(new_gate)
            for q in new_gate.get_qubits():
                new_node_map[q] = new_gate
    slice.node_map = new_node_map
    slice.dimension = 2
    return [slice, *new_slices]

def unroll_multicontrolled_slice(self, slices, sliceid=0, edge_slice_copy=None, debug=False, works_on_two=False, use_controls=None):
    '''
        TODO: Need to fix this function for the bug in the unrolling in the ancilla
        The idea should be to generate new gates that are @----X gates i.e. 1 control and act on two qubits that are oth not controlled on state 2 directly already.
        Maybe confirm if these gates can use 2 states in the first lapce? Maybe they dont exist with the new commute and delete operation (it uses rule 2 all the time)
    '''
    def unroll_a_gate(node, new_slices, flag=False, mk_new_slice=False, use_num_controls=None):        
        # print(node, "flag is ", flag)
        if len(node.get_qubits()) < 3:
            if mk_new_slice:
                new_slices.append(GateSlice(qubits=len(slices.node_map), gate=node, dimension=3, force_convert=True))
                return new_slices
        elif len(node.get_qubits()) == 3 and (not flag):
            if mk_new_slice:
                new_slices.append(GateSlice(qubits=len(slices.node_map), gate=node, dimension=3, force_convert=True))
                return new_slices
        else:
            # print(node)
            num_controls = len(node.get_control())
            controls_used = {}
            gates_list = []
            if flag:
                active_controls = len(node.get_control()) - 1
            else:
                active_controls = max(len(node.get_qubits()) // 2 - 1, 1)    
            while (num_controls > 2 and not flag) or (num_controls > 1 and flag):
                for q in node.get_qubits():
                    slices.node_map[q] = None
                    
                any_non_2_found = False
                for q in node.get_control():
                    if node.get_control_value_for_qubit(q) != 2:
                        any_non_2_found = True
                assert any_non_2_found
                for q in node.get_control():
                    if node.get_control_value_for_qubit(q) != 2:
                        assert active_controls != 0
                        controls_for_new_gate = [_ for _ in node.get_control()]
                        controls_for_new_gate.remove(q)
                        for q_c in controls_used:
                            controls_for_new_gate.remove(q_c)
                        controls_for_new_gate = controls_for_new_gate[:active_controls]
                        if node.get_control_value_for_qubit(q) == 0:
                            new_gate = LabelledNode(data={"qubits":[*controls_for_new_gate, q],
                                                "controls":[node.get_control_value_for_qubit(q_c) for q_c in controls_for_new_gate],
                                                "matrix":{"real": cirq.unitary(QutritX02Gate()), "imaginary": np.zeros((3, 3))},
                                                "status_rule_3":[]})
                        elif node.get_control_value_for_qubit(q) == 1:
                            new_gate = LabelledNode(data={"qubits":[*controls_for_new_gate, q],
                                                "controls":[node.get_control_value_for_qubit(q_c) for q_c in controls_for_new_gate],
                                                "matrix":{"real": cirq.unitary(QutritX12Gate()), "imaginary": np.zeros((3, 3))},
                                                "status_rule_3":[]})                        
                        gates_list.append(new_gate.deepcopy())
                        controls_used[q] = controls_for_new_gate
                        new_slices = new_slices + \
                            self.unroll_multicontrolled_slice(GateSlice(qubits=len(slices.node_map),
                                                                        gate=new_gate, dimension=3, force_convert=True),
                                                              sliceid=sliceid, edge_slice_copy=edge_slice_copy, works_on_two=True, use_controls=active_controls)
                        num_controls -= len(controls_for_new_gate)
                        node.set_control_value_for_qubit(q, 2)
                        for q_c in controls_for_new_gate:
                            node.delete_qubit(q_c)
                        break
            new_slices.append(GateSlice(qubits=len(slices.node_map), gate=node, dimension=3, force_convert=True))
            for gates in gates_list:
                new_slices = new_slices + \
                    self.unroll_multicontrolled_slice(GateSlice(qubits=len(slices.node_map), gate=gates, dimension=3, force_convert=True),
                                                        sliceid=sliceid, edge_slice_copy=edge_slice_copy, works_on_two=True)
        return new_slices

    nodes_realized = []
    new_slices = []
    for i, node in enumerate(slices.node_map):
        if node is not None and node.get_node_type().is_gateset():
            for q in node.get_qubits():
                slices.node_map[q] = None
            flag = works_on_two or edge_slice_copy.node_works_on_qutrit(sliceid, node)[0]
            if node not in nodes_realized:                
                for gates in node.getAllNodes():
                    new_slices = unroll_a_gate(gates, new_slices, flag=flag, mk_new_slice=True)
                nodes_realized.append(node)
        elif node is not None and node.get_node_type().is_gate() and i == node.get_target():
            flag = works_on_two or edge_slice_copy.node_works_on_qutrit(sliceid, node)[0]
            new_slices = unroll_a_gate(node, new_slices, flag=flag, use_num_controls=use_controls)
            nodes_realized.append(node)            
    slices.renormalize(dimension=3)
    slice_result = [slices, *new_slices]
    return slice_result

def convert_gateslices(self):
    self.slice_edge_collection.update_all_mappings()
    ancilla_maps = generate_ancilla_maps(self)
    if_ancilla = [x._is_ancilla for x in self.slice_collection[0].node_map]
    list_of_slices = copy.deepcopy(self.slice_collection)
    edge_data = self.reset(max(len(self.slice_collection[0].node_map), max(ancilla_maps)+1))
    self.dimension = 2
    for i, curr_slice in enumerate(list_of_slices):
        j = 0
        for new_slice in slice_convertor(curr_slice, ancilla_maps, edge_data=edge_data, slice_idx=i, ancilla_data=if_ancilla):
            j += 1
            self.add_slice(new_slice, no_edges=i==len(list_of_slices)-1)
    self.update_qubits()
    self.merge_greedily(binary=True)
    return

def unroll_multicontrolled_gates(self):
    self.slice_edge_collection.update_all_mappings()
    old_slices = copy.deepcopy(self.slice_collection) 
    slice_list = []
    curr_copy = self.reset(len(old_slices[0].node_map))
    for s, slices in enumerate(old_slices):
        new_slices = self.unroll_multicontrolled_slice(slices, sliceid=s, edge_slice_copy=curr_copy, debug=True)
        slice_list = [*slice_list, *new_slices]
        if s == len(old_slices) - 1:
            assert len(new_slices) == 1
            for i, node in enumerate(new_slices[0].node_map):
                if node is not None:
                    assert i in node.get_qubits()
            self.add_slice(new_slices[0], no_edges=True)
        else:        
            for new_slice in new_slices:
                for i, node in enumerate(new_slice.node_map):
                    if node is not None:
                        assert i in node.get_qubits()
                self.add_slice(new_slice)
    self.update_qubits()
    self.merge_greedily()
    self.dimension = 3
    curr = -1
    while(curr != len(self.slice_collection)):
        curr = len(self.slice_collection)
        self.get_all_matrix_mergable_gates(debug=False)
        self.prune_null_slices()
    self.get_all_matrix_mergable_gates(debug=False)
    self.prune_null_slices()
    self.get_all_matrix_mergable_gates(debug=False)
    self.prune_null_slices()       
    return

GateSliceCollection.convert_gateslices = convert_gateslices
GateSliceCollection.unroll_multicontrolled_gates = unroll_multicontrolled_gates
GateSliceCollection.unroll_multicontrolled_slice = unroll_multicontrolled_slice
