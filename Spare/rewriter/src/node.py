import ast
import numpy as np
from src.toffoli import *
from src.qutrit_decompositions import *
from rewriter.utils.graph_utils import *
import json
from rewriter.utils.node_gate_types import *
from rewriter.utils.edge_types import EdgeType

class Node:
    def __init__(self, data={}, node_info={}):
        self.qubits = []
        self.controls = []
        self.matrix = None
        self.label = ""
        self.root = False
        self.status = False
        self.sink = None
        self.status_rule_3 = None
        self.rule_2_insert = False
        self._is_ancilla = False
        if "label" in node_info:
          self.label = node_info["label"]
          assert len(self.label) != 0
        elif "labels" in node_info:
          self.label = node_info["labels"]
          assert len(self.label) != 0
        self.data_initialized = False
        self.dimension = None
        if "matrix" in node_info and "qubits" in node_info:
            data = node_info
        for d in data.keys():
            self.qubits = data["qubits"] if isinstance(data["qubits"], list) else [int(x.strip("\"")) for x in data["qubits"].split(",")]
            if len(self.qubits) > 1:
                self.controls = data["controls"] if "controls" in data else [1]*len(self.qubits[:-1])
            if "matrix" in data and isinstance(data["matrix"], dict):
                if (not isinstance(data.get("matrix")["real"], str)):
                    self.matrix = np.array(data.get("matrix")["real"], dtype=np.complex128) +\
                                  1j * np.array(data.get("matrix")["imaginary"], dtype=np.complex128)
                    self.dimension = self.matrix.shape[0]
                elif data.get("matrix")["real"] != "":
                    # print("matrix confirming ", data.get("matrix")["real"])
                    raise Exception
            elif "matrix" in data and isinstance(data["matrix"], str) and "real" in data["matrix"] and "imaginary" in data["matrix"]:
                matrix_string = data["matrix"]
                matrix_dict2 = ast.literal_eval(matrix_string)
                matrix_dict = ast.literal_eval(matrix_dict2)
                # Convert the 'real' and 'imaginary_part' strings to NumPy arrays
                real_part = np.array(eval(matrix_dict['real']))
                imaginary_part = np.array(eval(matrix_dict['imaginary_part']))
                # Combine the real and imaginary parts to form a complex NumPy array
                self.matrix = real_part + 1j * imaginary_part
                self.dimension = self.matrix.shape[0]
            elif "matrix" in data and isinstance(data["matrix"], str) and data["matrix"].strip("\"")=="None":
                self.matrix = None
            elif "matrix" in data and isinstance(data["matrix"], str) and data["matrix"].strip("\"")!="None":
                matrix_data = data["matrix"].strip("\"").replace(" ", "")
                # print(matrix_data)
                row_list = matrix_data.split("],[")
                matrix_data = []
                for r in row_list:
                    row_data = []    
                    for x in r.strip("\"").strip("]").strip("[").split(","):
                        if "j" not in x:
                            row_data.append(float(x))
                        else:
                            row_data.append(complex(x.replace("(", "").replace(")", "")))
                    matrix_data.append(row_data)
                
                self.matrix = np.array(matrix_data, dtype=np.complex128)
                self.dimension = self.matrix.shape[0]
            else:
                self.matrix = None
                self.dimension = None
            if "status_rule_3" in data:
                self.status_rule_3 = data["status_rule_3"]
            else:
                self.status_rule_3 = [None]*len(self.controls)
            self._is_ancilla = (self.matrix is None) and ("ancilla" in data.keys() and (data["ancilla"] == "true" or\
                data["ancilla"] == "True" or data["ancilla"] == True))
            self.data_initialized = True
        
    def serialize(self):
        data = {
            "qubits": self.qubits,
            "controls": self.controls,
            "matrix": {
                "real": np.real(self.matrix).tolist() if self.matrix is not None else "",
                "imaginary": np.imag(self.matrix).tolist() if self.matrix is not None else "",
            },
            "status_rule_3": self.status_rule_3,
            "ancilla": self._is_ancilla,
            "root": self.is_root(),
            "sink": self.is_sink()
        }
        json.dumps(data)
        return data

    @classmethod
    def read_in(cls, data):
        nodeObj = cls(data, {})
        nodeObj.set_root() if data["root"] else None
        nodeObj.set_sink() if data["sink"] else None
        return nodeObj
    
    def is_ancilla(self):
        assert self.get_node_type().is_qubit()
        return bool(self._is_ancilla)
    
    def inserted_node(self):
        self.rule_2_insert = True

    def get_status(self):
        return self.status

    def set_status(self, status):
        self.status = status
    
    def set_rule_3_status(self, status, qubit):
        assert qubit in self.qubits[:-1]
        assert isinstance(qubit, int)
        self.status_rule_3[self.qubits.index(qubit)] = status
    
    def get_rule_3_status(self):
        res = []
        for q in self.get_control():
            res.append(self.status_rule_3[self.qubits.index(q)])
        return res
    
    def get_rule_3_status_qubit(self, qubit):
        return self.status_rule_3[self.qubits.index(qubit)]
        
    def is_root(self):
        return self.root
    
    def is_sink(self):
        return self.sink
        
    def set_sink(self):
        assert self.get_node_type().is_qubit()
        self.sink = True
    
    def set_root(self):
        assert self.get_node_type().is_qubit()
        self.root = True

    def __str__(self):
        if self.matrix is None:
            return " qubits " + str(self.qubits) + " controls " + str(self.controls) + " anc: " + str(self._is_ancilla) + " " + str(self.status_rule_3)
        return " qubits " + str(self.qubits) + " controls " + str(self.controls) + "matrix " + np.array2string(self.matrix) + " anc: " + str(self._is_ancilla) + " " + str(self.status_rule_3)
    
    def set_matrix(self, matrix, assert_disabled=False):
        if matrix is not None and not assert_disabled:
            if self.dimension is not None:
                if self.dimension == 2 and matrix.shape[0] == 3 and np.allclose(matrix[2, 2], 1):
                    self.matrix = np.eye(2, dtype=np.complex128)
                    self.matrix = matrix[:2, :2]
                    return
                elif self.dimension == 3 and matrix.shape[0] == 2:
                    self.matrix = np.eye(3, dtype=np.complex128)
                    self.matrix[:2, :2] = matrix
                    return
                else:
                    # print(matrix)
                    # print(self.dimension)
                    # print(self.matrix)
                    assert matrix.shape[0] == self.dimension
            else:
                self.dimension = matrix.shape[0]
        elif assert_disabled and self.dimension != matrix.shape[0]:
            self.dimension = matrix.shape[0]
        self.matrix = matrix
    
    def update_qubit_targets(self, control, target):
        if control in self.qubits:
            self.controls.append(1)
            self.qubits.append(target)
        else:
            assert target in self.qubits
            self.add_control(control, 1)
    
    def add_control(self, control, value):
        # Sort qubits and add the control and qubit
        control_qubits = [*self.qubits[:-1], control]
        control_values = [*self.controls, value]
        zipped = list(zip(control_qubits, control_values))
        zipped = sorted(zipped, key = lambda x: x[0])
        qubits, controls = zip(*zipped)
        self.qubits = [*qubits, self.qubits[-1]]
        self.controls = list(controls)
    
    def broadcastQubits(self, controls):
        controls = controls.copy()
        for q in self.get_qubits():
            if q in controls:
                controls.remove(q)
        '''Are we not assumiong tghese are qubits'''
        values = 2**len(controls)
        gatelists = []
        for c in range(values):
            gate_copy = self.deepcopy(use_old_name=True)
            control_values = [int(x) for x in bin(c)[2:]]
            while len(control_values) < len(controls):
                control_values = [0, *control_values] 
            for q in range(len(controls)):
                gate_copy.add_control(controls[q], control_values[q])
            assert len(gate_copy.get_qubits()) == len(controls) + len(self.get_qubits())
            gatelists.append(gate_copy)
        assert len(gatelists) == 2**len(controls)
        return gatelists
    
    def set_control(self, controls):
        if isinstance(self.qubits, list):
            assert len(self.qubits) - 1 == len(controls)
            self.controls = controls

    def set_qubits(self, qubits):
        self.qubits = qubits
        self.status_rule_3 = ["n"] * len(self.qubits[:-1])
    
    def delete_qubit(self, qubit, control_only=False):
        if qubit not in self.qubits:
            raise Exception("qbit not in node")
        index = self.qubits.index(qubit)
        if control_only and index == len(self.qubits) - 1:
            raise Exception("qbit not a control")
        if index != len(self.qubits) - 1:
            self.controls.pop(index)
            self.label = self.label[1:]
            if len(self.status_rule_3) > 0 and len(self.status_rule_3) > index:
                self.status_rule_3.pop(index)
        self.qubits.remove(qubit)

    def delete_last_qubit(self, qubit, new_target=None, new_matrix=None, dim=3):
        assert qubit == self.get_target()
        dim = self.matrix.shape[0]
        if new_matrix is None:
            res = self.is_uniform_phase_gate({"status":EdgeType.Basis01}, {"status":EdgeType.Basis01})
            assert res[0]
            if len(self.qubits) == 1:
                return
            if res[1] == 1:
                self.status = "delete"
                return
        matrix = np.eye(dim, dtype=np.complex128)
        if new_target is not None:
            assert new_matrix is not None
            matrix = new_matrix
            self.delete_qubit(new_target)
            self.delete_qubit(qubit)
            self.qubits.append(new_target)
        else:
            matrix[self.get_control_values()[-1], self.get_control_values()[-1]] = res[1]
            self.controls.pop(-1)
            self.status_rule_3.pop(-1)
            self.delete_qubit(qubit)
        self.matrix = matrix
        self.label = self.label[1:]
        
    def set_control_for_qubit(self, qubit, control):
        self.controls[self.qubits.index(qubit)] = control
    
    def get_control(self):
        if not len(self.qubits) > 1:
            return []
        qubits_copy = self.qubits[:-1]
        qubits_copy.sort()
        return qubits_copy
    
    def replace_qubit(self, qubit1, qubit2):
        assert qubit1 in self.qubits
        index = self.qubits.index(qubit1)
        self.qubits[index] = qubit2
    
    def get_target(self):
        return self.qubits[-1]
    
    def set_target(self, q):
        self.qubits = [* self.qubits[:-1], q]
        assert len(self.qubits) - 1 == len(self.controls)
    
    def get_control_values(self):
        qubits_copy = self.qubits.copy()[:-1]
        qubits_copy.sort()
        control_arr = []
        for q in qubits_copy:
            control_arr.append(self.controls[self.qubits.index(q)])
        return control_arr
    
    def set_control_value_for_qubit(self, qubit, control):
        qubit = int(qubit)
        control = int(control)
        if qubit not in self.qubits[:-1]:
            raise Exception("Qubit not in controls")
        self.controls[self.qubits.index(qubit)] = control
    
    def get_control_value_for_qubit(self, qubit):
        if qubit not in self.qubits[:-1]:
            raise Exception("Qubit not in controls ")
        return self.controls[self.qubits.index(qubit)]
    
    def get_matrix(self):
        if self.matrix is None:
            return None
        return self.matrix.copy()
    
    def is_high_phase_gate(self):
        if self.get_node_type().is_gate() and np.allclose(self.matrix[:2, :2], np.eye(2)):
            return self
    
    def is_synchronized(self, gate, check_status=False):
        assert self.get_node_type().is_gate() and gate.get_node_type().is_gate()
        if self.get_control() == gate.get_control() and self.get_target() == gate.get_target() \
            and self.get_control_values() == gate.get_control_values():
            if not check_status:
                return True
            else:
                return self.get_rule_3_status() == gate.get_rule_3_status()
        return False
    
    def is_control_mismatch(self, gate):
        for q in gate.get_control():
            if q in self.get_control():
                c1 = gate.get_control_value_for_qubit(q)
                c2 = self.get_control_value_for_qubit(q)
                if c1 != c2:
                    return True
        return False

    def control_cancellable(self, gate):
        assert self.get_node_type().is_gate() and gate.get_node_type().is_gate()
        if self.qubits == gate.get_qubits() and if_mat_preserves_binary(self.get_matrix() @ gate.get_matrix().conjugate().T):
            return True
        return False

    def is_inverse(self, gate, phase_ok=True, check_status=True):
        assert self.get_node_type().is_gate() and gate.get_node_type().is_gate()
        if self.get_control() == gate.get_control() and self.get_target() == gate.get_target()\
              and self.get_control_values() == gate.get_control_values():
            if not check_status:
                if phase_ok:
                    return if_mat_preserves_binary(self.get_matrix() @ gate.get_matrix())
                else:
                    return np.allclose(self.get_matrix())
            else:
                if self.get_rule_3_status() != gate.get_rule_3_status():
                    return False
                else:
                    return self.is_inverse(gate, phase_ok=phase_ok, check_status=False)
        return False

    def is_immutable(self, state):
        if state is None:
            return False
        assert isinstance(state, int)
        # Also need to check no gate in between changes it
        assert self.get_node_type().is_gate()
        if self.get_gate_type() == GateType.PhaseGate:
            return True 
        vec = np.zeros((self.matrix.shape[0], 1), dtype=np.complex128)
        if state >= self.matrix.shape[0]:
            return True
        vec[state, 0] = 1
        if np.allclose(vec, np.absolute(self.matrix @ vec)):
            return True
        return False
    
    def is_control_subset(self, gate, qubits_and_controls=None):
        for x in self.get_control():
            if x not in gate.get_control():
                return False
            if qubits_and_controls is not None and self.get_control_value_for_qubit(x) != qubits_and_controls[x]:
                return False
            if qubits_and_controls is None and self.get_control_value_for_qubit(x) != gate.get_control_value_for_qubit(x):
                return False
        return True
    
    def qubits_match(self, gate):
        return set(self.get_qubits()[:-1]) == set(gate.get_qubits()[:-1]) and self.get_target() == gate.get_target()

    def get_control_mismatches(self, gate):
        assert self.qubits_match(gate)
        gate_controls = [gate.get_control_value_for_qubit(c) for c in gate.get_control()]
        self_controls = [self.get_control_value_for_qubit(c) for c in self.get_control()]
        mismatch_indexes = [i for i, (x, y) in enumerate(zip(self_controls, gate_controls)) if x != y]
        return [self.get_control()[i] for i in mismatch_indexes]

    def is_mergable(self, gate):
        if gate.get_qubits() == self.get_qubits() and gate.get_control_values() == self.get_control_values():
            return True
        return False

    def speculatively_swap_with_gate(self, gate, forward_dir=True, swap_through_gates=False,
                                     invalid_swap=False, desiredOp=None, all_gates=[], debug=False):
        result = ["", []]
        if gate.get_node_type().is_qubit():
            return True, result, [None, None]
        if gate.get_node_type().is_gateset():
            '''This is a modification since here the gates need to match perfectly for gate match'''
            '''Whereas otherwise all the gates in middle will modify the gates of the self object'''
            for i, g1 in enumerate(gate.getNodeSet()[0]):
                if self.get_target() in g1.get_control():
                    return False, ["", []], [None, None]
            if len(gate.getActiveSet()) == 1 and desiredOp=="PartialMerge":
                    g1 = gate.getActiveSet()[0]
                    if len(g1.get_control()) == len(self.get_control()) and g1.get_target() == self.get_target()\
                        and np.allclose(g1.get_matrix()@self.get_matrix().conjugate().T, np.eye(g1.get_matrix().shape[0])):
                        '''This can be partially merged here'''
                        return True, ["PartialMerge", []], [0, 0]
            for i, g1 in enumerate(gate.getActiveSet()):
                res1, res2, res3 = self.speculatively_swap_with_gate(g1, forward_dir=forward_dir, swap_through_gates=swap_through_gates,
                                                                     invalid_swap=invalid_swap, all_gates=all_gates, debug=debug)                    
                if not res1:
                    return res1, res2, [None, None]
            return True, ["", []], [None, None]
        if gate.get_node_type().is_gate():
            res, details = self.gate_swap_status(gate, forward_dir=forward_dir, swap_through_gates=swap_through_gates,
                                                 invalid_swap=invalid_swap, all_gates=all_gates, debug=debug)
            if debug:
                print(f"Forward dir {forward_dir}, invalid {invalid_swap}, swap through {swap_through_gates}, res {res}")
                if res:
                    print(gate)
                    print(self)
                    print("\t", details)
            return res, details, [None, None]
        raise Exception
    
    def gate_swap_status(self, gate, forward_dir=True, swap_through_gates=False, invalid_swap=False, all_gates=[], debug=False):
        result = ["", []]
        if self.qubits_match(gate) and np.allclose(self.get_matrix(), gate.get_matrix()):
            '''Same qubits exist now we should have a match'''
            control_mismatches = self.get_control_mismatches(gate)
            if len(control_mismatches) == 0:
                return True, ["GateMatch", []]
            if len(control_mismatches) == 1:
                return True, ["GateMerge", control_mismatches]
            if len(control_mismatches) == 2 and \
                np.allclose(np.eye(self.get_matrix().shape[0]), self.get_matrix() @ gate.get_matrix().conjugate().T):
                return True, ["ControlSwapper", control_mismatches]
            elif len(control_mismatches) == 3 and \
                np.allclose(np.eye(self.get_matrix().shape[0]), self.get_matrix() @ gate.get_matrix().conjugate().T):
                return True, ["ControlSwapper", control_mismatches]
            elif np.allclose(np.eye(self.get_matrix().shape[0]), self.get_matrix() @ gate.get_matrix().conjugate().T):
                return True, ["ControlSwapper", control_mismatches]
            else:
                print(control_mismatches)
                raise Exception
                return False, ["manyMistmatches", control_mismatches]
        '''This is for all the immutable cases where controls dont change'''
        if self.get_target() in gate.get_control() and self.is_immutable(gate.get_control_value_for_qubit(self.get_target())):
            return True, result
        if gate.get_target() in self.get_control() and gate.is_immutable(self.get_control_value_for_qubit(gate.get_target())):
            return True, result
        if self.is_control_mismatch(gate):
            '''Control mismatch definetly means we can perform the swap. But dont do it if it is not enabled?'''
            return True, result
        if gate.is_control_subset(self):
            '''The currrent gate is will be applied only when the gate being pushed through is also applied'''
            if debug:
                print("subset")
            if gate.get_target() in self.get_control():
                if gate.get_gate_type().state_permutation_gates():
                    copy_temp = gate.deepcopy()
                    if forward_dir:
                        copy_temp.set_matrix(copy_temp.get_matrix().conjugate().T)
                    self.control_update_qubit(copy_temp, gate.get_target())
                    return True, result
                else:
                    return False, result
            if gate.get_target() == self.get_target():
                if forward_dir:
                    self.set_matrix(gate.get_matrix() @ self.get_matrix() @ gate.get_matrix().conjugate().T)
                else:
                    self.set_matrix(gate.get_matrix().conjugate().T @ self.get_matrix() @ gate.get_matrix())
                return True, result
            if gate.get_target() not in self.get_qubits():
                return True, result
        if gate.get_target() in self.get_control():
            '''This swap can only happen if the middle gates are afggressively modified so careful;ly used to find bookending gates'''
            if debug:
                print("target in ctrl")
            res = invalid_swap and gate.get_gate_type().state_permutation_gates()
            return res, result
        if self.is_control_subset(gate) and self.get_target() in gate.get_control():
            '''Only if swapping through the gates is allowed then we can push a gate through something that modifies this current gate'''
            if debug:
                print("subset and targt in")
            res = swap_through_gates and self.get_gate_type().state_permutation_gates()
            return res, result
        if gate.get_target() not in self.get_qubits() and self.get_target() not in gate.get_qubits():
            '''gates dont interact'''
            if debug:
                print("no interatc")
            return True, result
        if self.get_target() in gate.get_control():
            '''If I am not a control subset but can be swapped when gate is broadcasted'''
            if debug:
                print("broadcast2")
            res = invalid_swap and self.get_gate_type().state_permutation_gates()
            return invalid_swap, result
        if forward_dir and gate.get_target() == self.get_target() and\
            np.allclose(gate.get_matrix().conjugate().T @ self.get_matrix() @ gate.get_matrix(), self.get_matrix()):
            '''If the gate hasn't changed we can commute it'''
            return True, result
        if not forward_dir and gate.get_target() == self.get_target() and\
            np.allclose(gate.get_matrix() @ self.get_matrix() @ gate.get_matrix().conjugate().T, self.get_matrix()):
            return True, result
        if not gate.get_gate_type().state_permutation_gates():
            return False, result
        return False, result
        #print(gate.get_gate_type(), gate.get_gate_type().state_permutation_gates())
        # if self.get_target() in gate.get_qubits():
        #     '''Maybe if the targets match but some other thing doesnt match, maybe make it redundant'''
        #     return False, result
        #print(self)
        #print(gate)
        #assert False        
    
    def swap_and_modify_gates(self, swapping_gate, forward_direction=True):
        if self.get_target() in swapping_gate.get_control():
            if self.is_control_subset(swapping_gate):
                if self.get_gate_type().state_permutation_gates():
                    swapping_gate.control_update_qubit(self, self.get_target(), forward_direction=forward_direction)
            elif self.is_control_mismatch(swapping_gate):
                pass
            else:
                assert False
        elif swapping_gate.get_target() in self.get_control():
            if swapping_gate.is_control_subset(self):
                if swapping_gate.get_gate_type().state_permutation_gates():
                    self.control_update_qubit(swapping_gate, swapping_gate.get_target(), forward_direction=forward_direction)
            elif swapping_gate.is_control_mismatch(self):
                pass
            else:
                assert False
        elif swapping_gate.get_target() == self.get_target():
            if swapping_gate.is_control_subset(self):
                pass
            elif swapping_gate.is_control_mismatch(self):
                pass
            elif np.allclose(self.get_matrix() @ swapping_gate.get_matrix(), swapping_gate.get_matrix() @ self.get_matrix()):
                pass
            else:
                assert False
        elif swapping_gate.get_target() not in self.get_qubits() and self.get_target() not in swapping_gate.get_qubits():
            pass
        else:
            assert False
        return swapping_gate        
        
    def is_control_superset(self, gate):
        return gate.is_control_subset(self)
    
    def cannot_apply_together_with(self, gate, qubits_and_controls):
        gates_match_control_cond = False
        for x in self.get_control():
            if x not in gate.get_control():
                return False
            if qubits_and_controls is None and self.get_control_value_for_qubit(x) != gate.get_control_value_for_qubit(x):
                gates_match_control_cond = True
            if qubits_and_controls is not None and self.get_control_value_for_qubit(x) != qubits_and_controls[x]:
                gates_match_control_cond = True
        return gates_match_control_cond
    
    def swap_or_merge_gate(self, node, dim=3):
        if self.get_node_type().is_qubit():
            return node, True
        if node is None and len(self.get_qubits()) == 1:
            copy = self.deepcopy()
            self.delete_gate()
            if not np.allclose(copy.get_matrix(), np.eye(dim)):
                return copy, True
            else:
                return None, True
        elif node is None and len(self.get_qubits()) > 1:
            return None, True
        elif node is not None and node.get_target() == self.get_target() and len(self.qubits) == 1:
            self.set_matrix(self.get_matrix() @ node.get_matrix())
            assert np.allclose(self.get_matrix() @ self.get_matrix().conjugate().T, np.eye(dim, dtype=np.complex128))
            self.delete_gate()
            if np.allclose(self.get_matrix(), np.eye(dim)):
                return None, True
            else:
                copy = self.deepcopy()
                return copy, True
        elif node is not None and node.get_target() == self.get_target():
            self.set_matrix(node.get_matrix().conjugate().T @ self.get_matrix() @ node.get_matrix())
            return node, True
        elif node is not None and node.get_target() in self.get_control():
            # print("NODE FOUND IN CONTROL", node.get_matrix(), node.get_gate_type().state_permutation_gates())
            # We will set control value here
            if node.get_gate_type().state_permutation_gates():
                self.control_update_qubit(node, node.get_target())
                return node, True
            else:
                return node, False 
        else:
            assert False
    
    def control_update_qubit(self, node, qbit=None, forward_direction=True):
        assert node.get_target() == qbit
        if node.get_gate_type() == GateType.MinusGate:
            if forward_direction:
                self.set_control_value_for_qubit(qubit=qbit, control=(self.get_control_value_for_qubit(qbit)+1)%3)
            else:
                self.set_control_value_for_qubit(qubit=qbit, control=(self.get_control_value_for_qubit(qbit)-1)%3)    
        elif node.get_gate_type() == GateType.PlusGate:
            if forward_direction:
                self.set_control_value_for_qubit(qubit=qbit, control=(self.get_control_value_for_qubit(qbit)-1)%3)
            else:
                self.set_control_value_for_qubit(qubit=qbit, control=(self.get_control_value_for_qubit(qbit)+1)%3)    
        elif node.get_gate_type() == GateType.X01Gate:
            if self.get_control_value_for_qubit(qbit) != 2:
                self.set_control_value_for_qubit(qubit=qbit, control=1 - self.get_control_value_for_qubit(qbit))
        elif node.get_gate_type() == GateType.X12Gate:
            if self.get_control_value_for_qubit(qbit) != 0:
                self.set_control_value_for_qubit(qubit=qbit, control=3 - self.get_control_value_for_qubit(qbit))
        elif node.get_gate_type() == GateType.X02Gate:
            if self.get_control_value_for_qubit(qbit) != 1:
                self.set_control_value_for_qubit(qubit=qbit, control=2 - self.get_control_value_for_qubit(qbit))
        elif node.get_gate_type() == GateType.PhaseGate:
            pass
        else:
            assert False        

    def get_qubits(self):
        q_copy = self.qubits[:-1]
        q_copy.sort()
        q_copy.append(self.qubits[-1])
        return q_copy
    
    def delete_control(self, qubit):
        self.delete_qubit(qubit, control_only=True)
    
    def node_retains_binary(self):
        if np.allclose(self.matrix[2, 2] * self.matrix[2, 2].conjugate(), 1):
            return True
        else:
            return False
        
    def to_delete(self):
        if self.status == "delete":
            return True
        return False

    def delete_gate(self):
        if self.get_node_type().is_qubit():
            assert False
        self.status = "delete"
    
    def delete_top_control(self):
        assert len(self.qubits) > 1
        self.delete_qubit(self.get_qubits()[0])
    
    def is_phase_gate(self):
        return_state = True
        for i in range(self.matrix.shape[0]):
            return_state = return_state and np.allclose(self.matrix[i, i] * self.matrix[i, i].conjugate().T, 1)
        return return_state
    
    def is_uniform_phase_gate(self, edges_prev, edges_nxt):
        assert isinstance(edges_nxt["status"], EdgeType)
        assert isinstance(edges_prev["status"], EdgeType)
        if self.is_phase_gate():
            if self.matrix.shape[0] == 2 and np.allclose(self.matrix[0, 0], self.matrix[1, 1]):
                return True, self.matrix[0, 0]
            if np.allclose(self.matrix[0, 0], self.matrix[1, 1]) and \
                np.allclose(self.matrix[1, 1], self.matrix[2, 2]):
                return True, self.matrix[0, 0]
            if np.allclose(self.matrix[0, 0], self.matrix[1, 1]) and \
                (edges_prev["status"] == EdgeType.Basis01 or edges_nxt["status"] == EdgeType.Basis01):
                return True, self.matrix[0, 0]
        return False, None

    def is_high_phase_gate(self, dim=3):
        if dim == 2 and self.get_matrix().shape[0] == 2:
            return False
        if np.allclose(self.matrix[0, 0], 1):
            if np.allclose(self.matrix[1, 1], 1):
                if np.allclose(self.matrix[2, 2] * self.matrix[2, 2].conjugate().T, 1):
                    return True
        return False
    
    def get_gate_type(self):
        m1 = self.get_matrix()
        if m1.shape[0] == 2 and if_mat_preserves_binary(m1 @ cirq.unitary(cirq.X)):
          return GateType.X01Gate
        if self.is_phase_gate():
          return GateType.PhaseGate
        if m1.shape[0] == 2:
          return GateType.QubitGate
        if if_mat_preserves_binary(m1 @ cirq.unitary(QutritX01Gate())):
          return GateType.X01Gate
        if if_mat_preserves_binary(m1 @ cirq.unitary(QutritX12Gate())):
          return GateType.X12Gate
        if if_mat_preserves_binary(m1 @ cirq.unitary(QutritX02Gate())):
          return GateType.X02Gate
        if if_mat_preserves_binary(m1 @ cirq.unitary(QutritMinusGate())):
          return GateType.PlusGate
        if if_mat_preserves_binary(m1 @ cirq.unitary(QutritPlusGate())):
          return GateType.MinusGate
        if np.allclose(m1[2, 2] * m1[2, 2].conjugate().T, 1):
          return GateType.basis01Gate
        if np.allclose(m1[0, 0] * m1[0, 0].conjugate().T, 1):
          return GateType.basis12Gate
        if np.allclose(m1[1, 1] * m1[1, 1].conjugate().T, 1):
          return GateType.basis02Gate
        return GateType.GeneralGate
    
    def get_label(self):
        return self.label
    
    # def set_label(self, label):
    #     self.label = label
    
    def remove_phase(self):
        mat = self.matrix
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if mat[i, j].real < -0.9:
                    mat[i, j] = -mat[i, j]
        self.matrix = mat
