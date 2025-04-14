import json
import numpy as np
import cirq
from src.toffoli import *
from rewriter.utils.node_gate_types import *
from rewriter.src.node import Node
from rewriter.utils.edge_types import EdgeType
from rewriter.utils.node_gate_types import GateOpType
from collections import Counter


def control_config(qubits, g, c, dimension=3):
    circuit = cirq.Circuit()
    if g is not None:
        if dimension == 3:
            if g.get_control_value_for_qubit(c) == 2:
                circuit.append(QutritX12Gate().on(qubits[c]))
            elif g.get_control_value_for_qubit(c) == 0:
                circuit.append(QutritX01Gate().on(qubits[c]))
            elif g.get_control_value_for_qubit(c) == 1:
                pass
            else:
                assert False
        elif dimension == 2:
            if g.get_control_value_for_qubit(c) == 0:
                circuit.append(cirq.X(qubits[c]))
            elif g.get_control_value_for_qubit(c) == 1:
                pass
            else:
                print("control val is ", g.get_control_value_for_qubit(c))
                assert False
    else:
        circuit.append(QutritX12Gate().on(qubits[c]))
    return circuit

# Define the custom gate in Cirq
class ArbitraryGate(cirq.Gate):
    def __init__(self, unitary, qubits, controls, name, dimension=2):
        super(ArbitraryGate, self)
        self._unitary = unitary
        self.qs = qubits
        self.cs = controls
        self._dimension = dimension
        self.set_name(name)
        # print("dim set is ", dimension)
        if self._unitary is not None:
            assert self._unitary.shape[0] == self._dimension ** (len(self.qs)-len(self.cs))
    
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return tuple([self._dimension]*len(self.qs))
    
    def set_name(self, label):
        self._display_name = label

    def _num_qubits_(self):
        return len(self.qs)  # This is a 3-qubit gate

    def _unitary_(self):
        return self._unitary

    def _circuit_diagram_info_(self, args):
        '''In case a label is defined we need to return the gate with the label'''
        if self._display_name is not None and isinstance(self._display_name, str):
            return ["["+self._display_name+"]"] * len(self.qs)
        elif self._display_name is not None and isinstance(self._display_name, list):
            result_name = []
            for x in self._display_name:
                if isinstance(x, str):
                    result_name.append(f"[{x}]")
                else:
                    result_name.append(f"({x})")
            return result_name
        else:
            '''Trust but verify'''
            super()
        
class LabelledNode(Node):
    def __init__(self, data={}, node_info={}):
        '''Need to switch over to use Labelled Nodes from Nodes'''
        super().__init__(data, node_info)
        self._display_name = None
    
    def setName(self, name):
        self._display_name = name
    
    def getName(self):
        return self._display_name
    
    def realizeNode(self, qubits, if_final_circ=False, breakdown=""):
        '''
            We should be realizing a gate only here
            This realizes a node should be defined for a labelled node
        '''
        assert not if_final_circ
        # self.dimension=2
        circuit = cirq.Circuit()
        if (self._display_name is not None) and len(breakdown)==0:
            # print("node set qubits are ", self.get_qubits(), qubits, self.get_control(),  self.get_qubits())
            qubits = [qubits[q] for q in self.get_qubits()]
            '''If a label is specified add the special gate'''
            # print("mat check is ", self.get_matrix())
            circuit.append(ArbitraryGate(self.get_matrix(), qubits, self.controls, self._display_name, self.dimension).on(* qubits))
            return circuit
        elif len(breakdown)==1:
            controls = [qubits[q] for q in self.get_control()]
            qubits = [qubits[q] for q in self.get_qubits()]
            '''If a label is specified add the special gate'''
            name = [c for c in self.get_control_values()]
            name.append("G")
            circuit.append(ArbitraryGate(self.get_matrix(), qubits, self.controls, name, self.dimension).on(* qubits))
            return circuit
        if self.get_matrix().shape[0] > self.dimension:
            print(self.get_matrix())
            if np.allclose(self.matrix[2, 2], 1):
                self.matrix = self.matrix[:2, :2]
            else:
                raise Exception
        control_qubits = []
        for c in self.get_control():
            control_qubits.append(qubits[c])
            circuit.append(control_config(qubits, self, c, dimension=self.dimension))
        if self.dimension == 2 and np.allclose(self.get_matrix(), cirq.unitary(cirq.X)):
            circuit.append(cirq.X(qubits[self.get_target()]).controlled_by(* control_qubits))
        elif self.dimension == 3 and np.allclose(self.get_matrix(), cirq.unitary(QutritX01Gate())):
            circuit.append(QutritX01Gate().on(qubits[self.get_target()]).controlled_by(* control_qubits))
        else:    
            circuit.append(cirq.MatrixGate(self.get_matrix(), qid_shape=[self.dimension]).on(qubits[self.get_target()]).controlled_by(* control_qubits))
        for c in self.get_control():
            circuit.append(control_config(qubits, self, c, dimension=self.dimension))
        return circuit
    
    @classmethod
    def NodeSetVisualizer(cls, display_name, gate_qubits, dimension):
        circuit = cirq.Circuit()
        circuit.append(ArbitraryGate(None, gate_qubits, [], display_name, dimension).on(* gate_qubits))
        return circuit
    
    def deepcopy(self, setname=None, use_old_name=False):
        new_node = LabelledNode()
        if setname is not None:
            new_node.setName(setname)
        elif use_old_name:
            new_node.setName(self._display_name)    
        new_node.set_qubits(self.qubits.copy())
        new_node.set_control(self.controls.copy())
        new_node.set_matrix(self.get_matrix())
        new_node.status_rule_3 = self.status_rule_3.copy()
        return new_node
    
    def return_inverse(self, use_old_name=False):
        assert self.get_node_type().is_gate()
        new_node = self.deepcopy(use_old_name=use_old_name)
        new_node.set_matrix(new_node.get_matrix().conjugate().T)
        if self.getName() is not None:
            new_node.setName(self.getName()+"*")
        return new_node
    
    def get_node_type(self):
        if self.matrix is None:
            return NodeType.Qubit
        return NodeType.Gate
    
    def reset_qubits(self, remapped_qubits, control_values):
        for i, q in enumerate(self.qubits):
            if remapped_qubits[q] is not None and remapped_qubits[q] != q:
                self.qubits[i] = remapped_qubits[q]
                if control_values[i] == 0 and self.controls[i] == 0:
                    self.controls[i] == 1
                elif control_values[i] == 1 and self.controls[i] == 0:
                    self.controls[i] == 0
                elif control_values[i] == 1 and self.controls[i] == 1:
                    self.controls[i] == 1
                elif control_values[i] == 0 and self.controls[i] == 1:
                    self.controls[i] == 1
                else:
                    assert False


class NodeSet:
    '''Takes a class of nodes'''
    def __init__(self, nodes=None, data={}):
        if nodes is None:
            nodes = []
            assert data is not {} and "gate data" in data
            for gateSet in data["gate data"]:
                preppedGateSet = []
                for gate in gateSet:
                    preppedGateSet.append(LabelledNode.read_in(gate))
                    preppedGateSet[-1].setName(gate["display_name"])
                    self.dimension = preppedGateSet[-1].get_matrix().shape[0]
                nodes.append(preppedGateSet)
        self.node_set = nodes
        for gateSet in nodes:
            for gate in gateSet:
                assert gate.get_node_type().is_gate()
                assert gate.getName() is not None
                self.dimension = gate.get_matrix().shape[0]
        gate_qubits = []
        for gateSet in self.node_set:
            for gate in gateSet:
                qubits = [x for x in gate.get_qubits() if x not in gate_qubits]
                gate_qubits.extend(qubits)
        self.qubits = gate_qubits
        self.check()

    def check(self):
        assert len(self.node_set) == 3
        if len(self.node_set[0]) != len(self.node_set[2]):
            for n in self.node_set:
                print(len(n))
                for i in n:
                    print(i)
                print("___________")
            assert False
        
        
    def serialize(self):
        '''Dump a set of nodes as once'''
        data = {"gate data":None}
        gate_data = []
        for gateSet in self.node_set:
            gateSetData = []
            for gate in gateSet:
                gateSetData.append({"qubits": gate.qubits, "controls": gate.controls,
                                  "matrix": {"real": np.real(gate.matrix).tolist() if gate.matrix is not None else "",
                                         "imaginary": np.imag(gate.matrix).tolist() if gate.matrix is not None else "",},
                                  "status_rule_3": gate.status_rule_3, "display_name": gate.getName(),
                                  "root":False, "sink":False})
            gate_data.append(gateSetData)
        assert len(gate_data) == 3
        assert len(gate_data[0]) == len(gate_data[2])
        data["gate data"] = gate_data
        json.dumps(data)
        self.check()
        return data

    @classmethod
    def read_in(cls, data):
        return cls(nodes=None, data=data)
    
    def getNodeSet(self):
        return self.node_set

    def getAllNodes(self):
        return [*self.node_set[0], * self.node_set[1], *self.node_set[2]]
    
    def getActiveSet(self):
        return self.node_set[1]
    
    def get_qubits(self):
        return self.qubits
    
    def realizeNode(self, qubits, if_final_circ=False, breakdown=""):
        assert not if_final_circ
        self.check()
        all_nodes = self.getAllNodes()
        if isinstance(all_nodes[0].getName(), list):
            new_name = [""]*len(self.qubits)
            for i, q in enumerate(all_nodes[0].get_qubits()):
                new_name[q] += all_nodest[0].getName()[i]
        else:
            new_display_name = all_nodes[0].getName()            
        gate_qubits = []
        for gates in all_nodes[1:]:
            if isinstance(new_display_name, str):
                new_display_name += "_" + gates.getName()
            else:
                for i, q in enumerate(gates.get_qubits()):
                    new_name[q] += "_" + all_nodes[0].getName()[i]
            nq = [x for x in gates.get_qubits() if x not in gate_qubits]
            gate_qubits.extend(nq)
        gate_qubits = [qubits[x] for x in gate_qubits]
        if breakdown == "":
            return LabelledNode.NodeSetVisualizer(new_display_name, gate_qubits, self.dimension)
        else:
            circuit = cirq.Circuit()
            for r in self.getAllNodes():
                circuit.append(r.realizeNode(qubits=qubits, breakdown=breakdown))
            return circuit
    
    def deepcopy(self):
        copied_nodes = []
        for nodelist in self.node_set:
            copied_nodes.append([x.deepcopy(use_old_name=True) for x in nodelist])
        return NodeSet(nodes=copied_nodes)

    def return_inverse(self):
        copied_nodes = []
        for nodelist in reversed(self.node_set):
            copied_nodes.append([x.return_inverse() for x in reversed(nodelist)])
        return NodeSet(nodes=copied_nodes)
    
    def get_node_type(self):
        '''Assume I am a gate'''
        return NodeType.GateSet
    
    def get_gate_type(self):
        self.check()
        oldgatetype = None
        for g in self.getActiveSet():
            if oldgatetype is None or oldgatetype == g.get_gate_type():
                gatetype = g.get_gate_type()
                oldgatetype = g.get_gate_type()
            elif (oldgatetype is None or oldgatetype.is_binary()) and g.get_gate_type().is_binary():
                gatetype = GateType.basis01Gate
                oldgatetype = GateType.basis01Gate
            elif (oldgatetype is None or oldgatetype.state_permuation_gate()) and g.get_gate_type().state_permuation_gate():
                gatetype = GateType.GeneralPermuationGate
                oldgatetype = GateType.GeneralPermuationGate
            else:
                return GateType.GeneralGate
        if len(self.getActiveSet()) == 0:
            return GateType.PhaseGate
        else:
            return gatetype
        assert False
    
    def is_phase_gate(self):
        self.check()
        for g in self.getActiveSet():
          if not g.is_phase_gate():
            return False
        return True
    
    def simplify_active_gates(self, prevSliceInfo, nextSliceInfo):
        self.check()
        prevSliceInfo = prevSliceInfo.deepcopy()
        nextSliceInfo = nextSliceInfo.deepcopy()
        for gates in self.node_set[0]:
            if prevSliceInfo.get_wire_status(gates.get_target()) == EdgeType.Basis01 and not gates.get_gate_type().is_binary():
                prevSliceInfo.set_wire_status(gates.get_target(), EdgeType.Basis012)
        for gates in reversed(self.node_set[1]):    
            if nextSliceInfo.get_wire_status(gates.get_target()) == EdgeType.Basis01 and not gates.get_gate_type().is_binary():
                nextSliceInfo.set_wire_status(gates.get_target(), EdgeType.Basis012)
        delete_gates = []
        for gates in self.getActiveSet():
            for c in gates.get_control():
                if gates.get_control_value_for_qubit(c) == 2 and \
                    (prevSliceInfo.get_wire_status(c) == EdgeType.Basis01 or \
                        nextSliceInfo.get_wire_status(c) == EdgeType.Basis01):
                    delete_gates.append(gates)
                    break
        for gates in delete_gates:
            self.node_set[1].remove(gates)
        self.check()

    def speculatively_swap_with_gate(self, gate, forward_dir=True, swap_through_gates=False, invalid_swap=False, desiredOp=None, all_gates=[], debug=False):
        self.check()
        result = ["", []]
        if gate.get_node_type().is_qubit():
            return True, result, [None, None]
        if gate.get_node_type().is_gateset():
            '''This is a modification since here the gates need to match perfectly for gate match'''
            '''Whereas otherwise all the gates in middle will modify the gates of the self object'''
            '''We can safely swap if gates in the middle can be swapped. We can merge if the bookends match'''
            res1, res2, res3 = self.resolve_gatesets(gate, forward_dir=forward_dir, desiredOp=desiredOp, debug=debug)
            self.check()
            return res1, res2, res3
            # "TODO: NEED TO FIX THIS SWAP AND PLACEMENT SINCE THE REWRITE TYPE AHS NOT BEEN DEFINED"
        if gate.get_node_type().is_gate():
            gate_list_set = self.getActiveSet()
            for gates in self.getNodeSet()[0]:
                res1, res2, res3 = gates.speculatively_swap_with_gate(gate, forward_dir=forward_dir, swap_through_gates=True,
                                                                invalid_swap=True, debug=debug)
                if not res1:
                    return res1, res2, [None, None]
            if not forward_dir:
                gate_list_set = reversed(gate_list_set)
            for gates in gate_list_set:
                res1, res2, res3 = gates.speculatively_swap_with_gate(gate, forward_dir=forward_dir, swap_through_gates=swap_through_gates,
                                                                invalid_swap=invalid_swap, all_gates=all_gates, debug=debug)
                # if len(self.getActiveSet()) == 1 and desiredOp=="PartialMerge":
                #     g1 = self.getActiveSet()[0]
                #     if len(g1.get_control()) == len(gates.get_control()) and g1.get_target() == gates.get_target()\
                #         and np.allclose(g1.get_matrix()@gates.get_matrix().conjugate().T, np.eye(g1.get_matrix().shape[0])):
                #         '''This can be partially merged here'''
                #         return True, ["PartialMerge", []], [0, 0]
                if not res1:
                    self.check()
                    return False, ["", []], [None, None]
                if not res1:
                    return res1, res2, [None, None]
            '''We are not allowed to match a NodeSet with a single gate unfortunately but we can perform a swap'''
            self.check()
            return True, ["", []], [None, None]
    
    def resolve_gatesets(self, gateset, forward_dir, desiredOp, debug=False):
        # First check if the bookends are same:
        if forward_dir:
            pos1, pos2 = (2, 0)
        else:
            pos1, pos2 = (0, 2)
        new_gates = []        
        for i, g1 in enumerate(self.getActiveSet()):
            for j, g2 in enumerate(gateset.getActiveSet()):
                res1, res2, res3 = g1.speculatively_swap_with_gate(g2, swap_through_gates=True, forward_dir=forward_dir, debug=debug)
                if desiredOp=="PartialMerge" and len(self.getMismatchedBookends(gateset, forward_dir=forward_dir)) > 0:
                    if len(g1.get_control()) == len(g2.get_control()) and g1.get_target() == g2.get_target()\
                        and np.allclose(g1.get_matrix()@g2.get_matrix().conjugate().T, np.eye(g1.get_matrix().shape[0])):
                        '''This can be partially merged here'''
                        return True, ["PartialMerge", []], [i, j]
                if not res1:
                    self.check()
                    return False, ["", []], [None, None]
                elif res1 and len(res2[0]) != 0 and self.is_mergable(gateset, forward_dir=forward_dir, debug=debug):
                    '''In the compile and rewrite figure out how this optimization happens'''
                    self.check()
                    return res1, res2, (i, j)
                elif res1 and len(res2[0]) != 0 and not self.is_mergable(gateset, forward_dir=forward_dir, debug=debug):
                    return True, ["", []], [None, None]
        return True,  ["", []], [None, None]
    
    def swap_and_modify_gates(self, swapping_gate, forward_direction=True):
        self.check()
        if forward_direction:
            for gates in self.getAllNodes():
                swapping_gate = gates.swap_and_modify_gates(swapping_gate, forward_direction=True)
        else:
            for gates in reversed(self.getAllNodes()):
                swapping_gate = gates.swap_and_modify_gates(swapping_gate, forward_direction=False)
        # if self.getActiveSet()[0].get_target() in swapping_gate.get_control():
        #     curr_gate = self.getActiveSet()[0]
        #     if curr_gate.is_control_subset(swapping_gate):
        #         if curr_gate.get_gate_type().state_permutation_gates():
        #             swapping_gate.control_update_qubit(curr_gate, curr_gate.get_target(), forward_direction=forward_direction)
        #     elif curr_gate.is_control_mismatch(swapping_gate):
        #         pass
        #     else:
        #         assert False
        # elif swapping_gate.get_target() in self.get_control():
        #     if swapping_gate.is_control_subset(self):
        #         if swapping_gate.get_gate_type().state_permutation_gates():
        #             self.control_update_qubit(swapping_gate, swapping_gate.get_target(), forward_direction=forward_direction)
        #     elif swapping_gate.is_control_mismatch(self):
        #         pass
        #     else:
        #         assert False
        # elif swapping_gate.get_target() == self.get_target():
        #     if swapping_gate.is_control_subset(self):
        #         pass
        #     elif swapping_gate.is_control_mismatch(self):
        #         pass
        #     elif np.allclose(self.get_matrix() @ swapping_gate.get_matrix(), swapping_gate.get_matrix() @ self.get_matrix()):
        #         pass
        #     else:
        #         assert False
        # elif swapping_gate.get_target() not in self.get_qubits() and self.get_target() not in swapping_gate.get_qubits():
        #     pass
        # else:
        #     assert False
        self.check()
        return swapping_gate
    
    def is_mergable(self, node, forward_dir=False, debug=False):
        self.check()
        node.check()
        matching_nodes1 = []
        matching_nodes2 = []
        if not forward_dir:
            pos = (self.getNodeSet()[0], reversed(node.getNodeSet()[2]))
        else:
            pos = (reversed(self.getNodeSet()[2]), node.getNodeSet()[0])
        for g1 in pos[0]:
            for g2 in pos[1]:
                if g1.is_mergable(g2) and g2 not in matching_nodes2:
                    matching_nodes1.append(g1)
                    matching_nodes2.append(g2)
                    break
                if not g1.speculatively_swap_with_gate(g2, forward_dir=forward_dir)[0]:
                    '''Cant safely swap and merge'''
                    if debug:
                        print("\t\t\t\t\t", g1.get_qubits(), g2.get_qubits())
                    return False
        self.check()
        return (len(matching_nodes1) == len(self.getNodeSet()[0])) and (len(matching_nodes2)==len(node.getNodeSet()[2]))
    
    def getMismatchedBookends(self, gateset, forward_dir=True):
        matching_nodes1 = []
        matching_nodes2 = []
        if gateset.get_node_type().is_gateset():
            for g1 in reversed(self.getNodeSet()[2]):
                gates_mergable = True
                for g2 in gateset.getNodeSet()[0]:
                    if (g2 not in matching_nodes2) and (g1 not in matching_nodes1) and gates_mergable:
                        if g1.is_mergable(g2):
                            matching_nodes1.append(g1)
                            matching_nodes2.append(g2)
                        elif not g1.speculatively_swap_with_gate(g2, forward_dir=forward_dir)[0]:
                            gates_not_mergable = False
            list1 = [g1 for g1 in self.getNodeSet()[2] if g1 not in matching_nodes1]       
            list1.extend([g1 for g1 in gateset.getNodeSet()[2] if g1 not in matching_nodes2])
            return list1
        return self.getNodeSet()[2]

    def mergeMismatchedBookends(self, gateset, mismatchBooks, gateids, forward_dir=True):
        def ctrl_pruner(ctrlLists):
            if ctrlLists.shape[0] > 1:
                result = []
                for i in range(ctrlLists.shape[0]):
                    if not (np.any(np.all(ctrlLists[:i, :] == ctrlLists[i, :], axis=1)) or\
                        np.any(np.all(ctrlLists[i+1:, :] == ctrlLists[i, :], axis=1))):
                        result.append(ctrlLists[i, :])
                return np.asarray(result)
            return np.asarray(ctrlLists)
        
        def ctrl_deleter(ctrlLists, qubitLists, newTarget):
            if ctrlLists.shape[1] > 1:
                non_identical_columns = [col for col in range(ctrlLists.shape[1]) if not np.all(ctrlLists[:, col] == ctrlLists[0, col])]
                # print("\t\t", non_identical_columns)
                newqubitLists = []
                for i in non_identical_columns:
                    newqubitLists.append(qubitLists[i])
                # print("\t\tNew qubits:", newqubitLists)
                for i in list(set(qubitLists) - (set(newqubitLists))):
                    newTarget.delete_qubit(i)
                return ctrlLists[:, non_identical_columns], newqubitLists, newTarget
            return ctrlLists, qubitLists, newTarget

        def recursively_add_gates(target_gate, qubitLists, ctrlLists, bookends_dir=None, direction="forward"):
            newGate = target_gate.deepcopy(use_old_name=True)
            newTarget = target_gate.deepcopy(use_old_name=True)
            # print("shape: ", ctrlLists.shape)
            # , ctrlLists[0])
            ctrlLists = ctrl_pruner(ctrlLists)
            ctrlLists, qubitLists, newTarget = ctrl_deleter(ctrlLists, qubitLists, newTarget)
            # print("ctrl del: ", ctrlLists, qubitLists, newTarget)
            if ctrlLists.shape[0] == 1:
                for i in range(ctrlLists.shape[1]):
                    newGate.set_control_value_for_qubit(qubitLists[i], ctrlLists[0, i])
                # print("new gate is- ", newGate)
                return [newGate]
            qidx = ctrlLists.shape[1]-1
            curr_qubit = qubitLists[qidx]
            counts = Counter(ctrlLists[:, qidx])
            common_control = counts.most_common(1)[0][0]
            # print("common control ", common_control)
            newGate.set_control_value_for_qubit(qubitLists[qidx], common_control)
            for c in newGate.get_control():
                if c in qubitLists[:qidx]:
                    newGate.delete_qubit(c)
            subfield = []
            for gtIdx in range(ctrlLists.shape[0]):
                if common_control != ctrlLists[gtIdx, qidx]:
                    '''This matches the control value here
                    We can skip this control beinga dded into the middle gate.'''
                    subfield.append(ctrlLists[gtIdx, :qidx])
            subfield = np.array(subfield)
            # print("subfild info ", subfield)
            newTarget.delete_last_qubit(newTarget.get_target(), qubitLists[qidx], new_matrix=np.array([[0, 1], [1, 0]], dtype=np.complex128))
            # print("new gate is: ", newGate)
            if bookends_dir is None:
                return [recursively_add_gates(newTarget, qubitLists[:qidx], subfield, bookends_dir="Forward"), [newGate], recursively_add_gates(newTarget, qubitLists[:qidx], subfield, bookends_dir="Backward")]
            elif bookends_dir=="Forward":
                return [*recursively_add_gates(newTarget, qubitLists[:qidx], subfield, bookends_dir="Forward"), newGate]
            elif bookends_dir=="Backward":
                return [newGate, *recursively_add_gates(newTarget, qubitLists[:qidx], subfield, bookends_dir="Backward")]
            else:
                raise Exception
        # print("START GATESET MERGES HERE")
        controls = [g.get_qubits() for g in mismatchBooks]
        controlList = []
        for cs in controls:
          controlList = sorted(list(set(controlList).union(set(cs))))
        # print("check control list ", controlList)
        
        newTargetGates = []
        newGates = []
        if gateset.get_node_type().is_gateset():
            _gate_list = gateset.getActiveSet()
        else:
            gateset.setName("dummy")
            _gate_list = [gateset]        
        for i, g1 in enumerate(set(_gate_list).union(self.getActiveSet())):
            controlList = sorted(list(set(controlList).union(set(g1.get_control()))))
        for i, g1 in enumerate(_gate_list):
            if i == gateids[1]:
                newTargetGates.extend(g1.broadcastQubits(controlList))    
            else:    
                newGates.extend(g1.broadcastQubits(controlList)) 
        for i in range(len(newTargetGates)):
            for gate in mismatchBooks:
                res = newTargetGates[i].speculatively_swap_with_gate(gate, forward_dir=forward_dir)
                # print("RES CHECK ", res)
        qubits = newTargetGates[0].get_control()
        for j, g2 in enumerate(self.getActiveSet()):
          if j == gateids[0]:
            newTargetGates.extend(g2.broadcastQubits(qubits))    
          else:
            newGates.append(g2)

        ctrlLists = []
        for gates in newTargetGates:
            ctrlLists.append(gates.get_control_values())
        # print("control Values after swapping are")
        # print(ctrlLists)        
        ctrlLists = np.asarray(ctrlLists)
        '''Deleted controls of interest have been prepped'''
        result = recursively_add_gates(newTargetGates[0], qubits, ctrlLists, bookends_dir=None)
        if len(newGates) == 0:
            for g in self.getAllNodes():
                assert not isinstance(g, list)
            self.node_set[0] = [*self.node_set[0], *result[0]]
            self.node_set[2] = [*result[2], *self.node_set[2]]
            self.addActiveGate(result[1], position=gateids[1])
            self.delActiveGate(position=gateids[1])
            for g in self.getAllNodes():
                assert not isinstance(g, list)
        else:
            raise NotImplementedError

    def merge(self, node):
        self.check()
        node.check()
        new_active_nodes = []
        for g1 in self.getNodeSet()[0]:
            for g2 in node.getNodeSet()[2]:
                if g1.is_mergable(g2):
                    if not np.allclose(g1.get_matrix(), g2.get_matrix().conjugate().T):
                        new_active_nodes.append(g1.set_matrix(g1.get_matrix() @ g2.get_matrix()))
                    break
        if len(new_active_nodes) == 0:
            for i, g1 in enumerate(self.getActiveSet()):
                for g2 in node.getActiveSet():
                    if g1.is_mergable(g2):
                        if not np.allclose(g1.get_matrix(), g2.get_matrix().conjugate().T):
                            new_active_nodes.append(g1.set_matrix(g1.get_matrix() @ g2.get_matrix()))
                        break
            if len(new_active_nodes) == self.getActiveSet():
                self.node_set[1] = new_active_nodes
            else:
                self.node_set[1] = [*node.getActiveSet(), * self.getActiveSet()]
        else:
            self.node_set[1] = [*node.getActiveSet(), * new_active_nodes, * self.getActiveSet()]
        self.node_set[0] = node.getNodeSet()[0]
        self.check()

    def addActiveGate(self, gates, position):
        self.node_set[1] = [*self.node_set[1][:position+1], *gates, *self.node_set[1][position+1:]]
        
    def delActiveGate(self, position=None, gates=None):
        self.check()
        if position is not None:
            self.node_set[1].pop(position)
        elif gates is not None:
            self.node_set[1].remove(gates)
        else:
            raise Exception
        self.check()
    
    def prune_if_null(self):
        allgates = []
        for i, gates in enumerate(self.getActiveSet()):
            if np.allclose(gates.get_matrix(), np.eye(gates.get_matrix().shape[0])):
                self.delActiveGate(gates=gates)
        self.pruneBookends()
        self.check()

    def pruneBookends(self):
        self.check()
        deleted_gates = []
        all_gates = self.getNodeSet()[0]
        controls_considered = []
        for gate in self.getNodeSet()[1]:
            controls_considered = list(set(controls_considered).union(gate.get_control()))
        gates_to_delete = []
        for gate in reversed(self.getNodeSet()[0]):
            if gate.get_target() not in controls_considered:
                gates_to_delete.append(gate)
            elif np.allclose(gate.get_matrix(), np.eye(gate.get_matrix().shape[0])):
                gates_to_delete.append(gate)
            controls_considered = list(set(controls_considered).union(gate.get_control()))
        for gate in gates_to_delete:
            self.node_set[0].remove(gate)
        for gate in self.getNodeSet()[1]:
            controls_considered = list(set(controls_considered).union(gate.get_control()))
        gates_to_delete = []
        for gate in self.getNodeSet()[2]:
            if gate.get_target() not in controls_considered:
                gates_to_delete.append(gate)
            elif np.allclose(gate.get_matrix(), np.eye(gate.get_matrix().shape[0])):
                gates_to_delete.append(gate)    
            controls_considered = list(set(controls_considered).union(gate.get_control()))
        for gate in gates_to_delete:
            self.node_set[2].remove(gate)
        self.check()
    
    def reset_qubits(self, remapped_qubits, control_values):
        for gates in self.getAllNode():
            gates.reset_qubits(remapped_qubits, control_values)
