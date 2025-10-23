import collections
import json
from rewriter.utils.edge_types import EdgeType
from rewriter.utils.node_gate_types import GateOpType

def filter_nodes(gateset, q):
    res = []
    for n in gateset:
        if n.get_node_type().is_gateset():
            res.append(n)
        elif n.get_target() == q:
            res.append(n)
    return res

class EdgeCollection:
    def __init__(self, slice, status=None):
        self.wire_status = []
        self.qubit = []
        # print(slice)
        if slice is not None:
            for qubits in range(slice.get_num_qubits()):
                # print(slice.get_node(qubits), end=" ")
                if slice.get_node(qubits) is not None and slice.get_node(qubits).get_node_type().is_qubit() and slice.get_node(qubits).is_ancilla():
                    self.wire_status.append(EdgeType.Basis0)
                else:
                    self.wire_status.append(EdgeType.Basis01)
                self.qubit.append(qubits)
        # print("")
        # print(self.wire_status)
        if status is not None:
            self.wire_status = status
            self.qubit = list(range(len(status)))
    
    def __eq__(self, node):        
        if self.wire_status != node.wire_status:
            return False
        if self.qubit != node.qubit:
            return False
        return True

    def serialize(self):
        data = {}
        data["qubits"] = len(self.wire_status)
        for i in range(len(self.wire_status)):
            data["status_" + str(i)] = str(self.wire_status[i])
        json.dumps(data)
        return data
    
    @classmethod
    def read_in(cls, data):
        def convert_string_to_enum(in_str):
            return EdgeType.from_string(in_str)
        
        qubits = data["qubits"]
        arr = []
        for i in range(qubits):
            arr.append(convert_string_to_enum(data["status_" + str(i)]))
        return cls(None, arr)
    
    @classmethod
    def create_copy(cls, data):
        return cls(None, status=data[""])        
    
    def __str__(self):
        return_string = ""
        for i, q in enumerate(self.wire_status):
            if i == 0:
                return_string += str(i) + " : " + str(q)     
            else:
                return_string += ", " + str(i) + " : " + str(q) 
        return_string += "\n"
        return return_string
    
    def get_wire_qubit(self, qubit):
        assert qubit in self.qubit
        return self.qubit[qubit]
    
    def get_wire_status(self, qubit):
        assert qubit in self.qubit
        return self.wire_status[qubit]

    def get_wire_mapping(self, qubit):
        assert qubit in self.qubit
        return self.wire_mapping[qubit]
    
    def set_wire_mapping(self, qubit, matrix):
        self.wire_mapping[qubit] = matrix
    
    def set_wire_status(self, qubit, status):
        self.wire_status[qubit] = status

    def get_qubits(self):
        return self.qubit

    def get_num_qubits(self):
        assert len(self.qubit) == len(self.wire_status)
        return len(self.wire_status)

    def set_default_mapping(self, sl, debug=False):
        if debug:
            print(sl)
        for qubits in range(len(self.wire_status)):
            if debug:
                print(sl.get_node(qubits))
                print("anc: ", sl.get_node(qubits)._is_ancilla)
                print(sl.get_node(qubits) is not None, sl.get_node(qubits).get_node_type().is_qubit(), sl.get_node(qubits).is_ancilla())
                
            if sl.get_node(qubits) is not None and sl.get_node(qubits).get_node_type().is_qubit() and sl.get_node(qubits).is_ancilla():
                self.wire_status[qubits] = EdgeType.Basis0
            else:
                self.wire_status[qubits] = EdgeType.Basis01
    
    def deepcopy(self):
        new_copy = EdgeCollection(slice=None)
        new_copy.wire_status = self.wire_status.copy()
        new_copy.qubit = self.qubit.copy()
        return new_copy
    
    def merge(self, old_slice):
        assert self.qubit == old_slice.qubit
        for q in range(len(self.qubit)):
            if old_slice.get_wire_status(q) == EdgeType.Basis0 or self.get_wire_status(q) == EdgeType.Basis0:
                self.set_wire_status(q, EdgeType.Basis0)
            elif old_slice.get_wire_status(q) == EdgeType.Basis01 or self.get_wire_status(q) == EdgeType.Basis01:
                self.set_wire_status(q, EdgeType.Basis01)

class EdgeSliceNodes:
    def __init__(self, gate_slice_collection=None, list_of_edges=None):
        self.list_of_edges = []
        self.gate_slice_collection = gate_slice_collection
        self.num_slices = 0
        self.del_arr = []
        if list_of_edges is not None:
            self.list_of_edges = list_of_edges
    
    def set_gate_slice_collection(self, gate_slice_collection):
        self.gate_slice_collection = gate_slice_collection

    def serialize(self):
        data = collections.defaultdict()
        data["edges_num"] = len(self.list_of_edges)
        for i in range(len(self.list_of_edges)):
            data[i] = self.list_of_edges[i].serialize()
        json.dumps(data)
        return data

    @classmethod
    def read_in(cls, data):
        edgelist = []
        for i in range(data["edges_num"]):
            edgelist.append(EdgeCollection.read_in(data[str(i)]))
        return cls(None, edgelist)

    @classmethod
    def create_copy(cls, gateslices, list_of_edges):
        return cls(gate_slice_collection=gateslices, list_of_edges=list_of_edges)
    
    def print(self):
        for l in self.list_of_edges:
            print(l)
    
    def merge_index(self, idx2, idx1, dont_create_mappings=False):
        if dont_create_mappings:
            self.del_arr.append(idx1)
            return
        edge_slice2 = self.list_of_edges[idx2]
        edge_slice1 = self.list_of_edges[idx1]
        edge_slice2.merge(edge_slice1)
        self.del_arr.append(idx1)
        return
    
    def prune_edgeslices(self):
        new_list_edges = []
        for i in range(len(self.list_of_edges)):
            if i not in self.del_arr:
                new_list_edges.append(self.list_of_edges[i])
        self.list_of_edges = new_list_edges
        self.del_arr = []

    def add_slice(self, slice, idx):
        edge_collection = EdgeCollection(slice=slice)
        self.list_of_edges.append(edge_collection)

    def add_slice_after_idx(self, slice, idx):
        edge_collection = EdgeCollection(slice=slice)
        self.list_of_edges = [* self.list_of_edges[:idx+1], edge_collection, * self.list_of_edges[idx+1:]]
    
    def get_slice_edge(self, idx):
        return self.list_of_edges[idx]
    
    def delete_slice(self, idx):
        if idx > 0:
            idx = idx - 1
            self.list_of_edges.pop(idx)
        else:
            pass
    
    def set_edge_collection(self, idx=0):
        collection = self.list_of_edges[idx]
        if idx == 0:
            collection.set_default_mapping()
            return
        oldcollection = self.list_of_edges[idx - 1]
        gateslice = self.gate_slice_collection.get_slice(idx)
        for q in gateslice.get_qubits():
            if gateslice.get_node(q) is None:
                collection.set_wire_status(q, oldcollection.get_wire_status(q))
                continue
            if gateslice.get_node(q).get_node_type().is_gate() and (gateslice.get_node(q).get_gate_type().is_binary()):
                gate_type_binary = True
            else:
                gate_type_binary = False
            if gateslice.get_node(q).get_node_type().is_gateset():
                active_gates = gateslice.get_node(q).getActiveSet()
            else:
                active_gates = [gateslice.get_node(q)]
            curr_history = []
            for curr_gate in active_gates:
                curr_history.append(curr_gate)
                if curr_gate.get_node_type().is_qubit():
                    pass
                elif curr_gate.is_phase_gate():
                    collection.set_wire_status(q, oldcollection.get_wire_status(q))
                elif q == curr_gate.get_target() and gate_type_binary and\
                    (oldcollection.get_wire_status(q) == EdgeType.Basis0 or oldcollection.get_wire_status(q) == EdgeType.Basis1):
                    collection.set_wire_status(q, EdgeType.basis01)    
                elif q == curr_gate.get_target() and oldcollection.get_wire_status(q) == EdgeType.Basis01 and gate_type_binary:
                    collection.set_wire_status(q, oldcollection.get_wire_status(q))
                elif q == curr_gate.get_target() and curr_gate.get_gate_type().turns_state_to_t and\
                    oldcollection.get_wire_status(q) == EdgeType.Basis01:
                    print("SET TO 2")
                    print(curr_gate)
                    print(curr_gate.get_gate_type())
                    collection.set_wire_status(q, EdgeType.Basis012)
                elif q == curr_gate.get_target():
                    # Fetch all the pre and post gates
                    q_pre, q_post = self.gate_slice_collection.get_nodes_before_after_gate(gateslice.get_node(q))
                    gateset = q_pre
                    '''TODO: Modify this to capture all gates before all gates'''
                    gateset.append(curr_history)
                    for gates in gateset:
                        if gates.get_gate_type().turns_state_to_t():
                            collection.set_wire_status(q, EdgeType.Basis012)
                            break
                    if collection.get_wire_status(q) != EdgeType.Basis012:
                        '''Above set has not taken place yet'''
                        collection.set_wire_status(q, EdgeType.Basis01)
                else:
                    if (oldcollection.get_wire_status(q) == EdgeType.Basis01):
                        collection.set_wire_status(q, EdgeType.Basis01)
                    else:
                        collection.set_wire_status(q, oldcollection.get_wire_status(q))

    def update_all_mappings(self, debug=False):
        if len(self.gate_slice_collection.slice_collection) - 1 != len(self.list_of_edges):
            print(len(self.gate_slice_collection.slice_collection))
            print(len(self.list_of_edges))
            assert False
        self.deleting_pairs = {}
        for i in range(self.list_of_edges[0].get_num_qubits()):
            self.deleting_pairs[i] = []
        for s_idx in range(len(self.gate_slice_collection.slice_collection) - 1):
            self.set_mapping(s_idx, debug=debug)
        if debug:
            print("______________________________________")
        self.deleting_pairs = {}
        for i in range(self.list_of_edges[0].get_num_qubits()):
            self.deleting_pairs[i] = []
        for s_idx in reversed(range(len(self.gate_slice_collection.slice_collection) - 1)):
            self.set_mapping(s_idx, mode="reverse", debug=debug)

    def set_mapping(self, idx, mode="forward", debug=False):
        if idx >= len(self.list_of_edges):
            raise Exception()
        collection = self.list_of_edges[idx]
        if idx == 0:
            collection.set_default_mapping(self.gate_slice_collection.slice_collection[0], debug=debug)
            return
        if idx == len(self.gate_slice_collection.slice_collection) - 2:
            collection.set_default_mapping(self.gate_slice_collection.slice_collection[idx+1])
            return
        if mode == "forward":
            oldcollection = self.list_of_edges[idx - 1]
            gateslice = self.gate_slice_collection.get_slice(idx)
        elif mode == "reverse":
            oldcollection = self.list_of_edges[idx + 1]
            gateslice = self.gate_slice_collection.get_slice(idx + 1)
        else:
            assert False
        # Get current copy for the edge configuration
        origconfig = collection.deepcopy()
        if debug:
            print("old config --- ", idx-1, oldcollection)
            print("orig values are ", idx, origconfig)
            print(gateslice)
        for q in gateslice.get_qubits():
            if gateslice.get_node(q) is None:
                if oldcollection.get_wire_status(q).contains(origconfig.get_wire_status(q)):
                    collection.set_wire_status(q, origconfig.get_wire_status(q))
                elif origconfig.get_wire_status(q).contains(oldcollection.get_wire_status(q)):
                    collection.set_wire_status(q, oldcollection.get_wire_status(q))
                else:
                    collection.set_wire_status(q, oldcollection.get_wire_status(q))                    
                continue
            if gateslice.get_node(q).get_node_type().is_qubit():
                # Expand this new system for the edge marking correctly
                if gateslice.get_node(q).is_ancilla():
                    collection.set_wire_status(q, EdgeType.Basis0)
                else:
                    collection.set_wire_status(q, EdgeType.Basis01)
                continue
            gate_type_binary = gateslice.get_node(q).get_gate_type().is_binary()
            if gateslice.get_node(q).get_node_type().is_gateset():
                active_gates = [gateslice.get_node(q)]
            else:
                active_gates = [gateslice.get_node(q)]
            curr_history = []
            for gate_id, curr_gate in enumerate(active_gates):
                curr_history.append(curr_gate)
                if curr_gate.is_phase_gate():
                    collection.set_wire_status(q, oldcollection.get_wire_status(q))
                elif not curr_gate.get_node_type().is_gateset() and q == curr_gate.get_target() and \
                    origconfig.get_wire_status(q).is_ancilla_state():
                    collection.set_wire_status(q, EdgeType.Basis0)
                elif not curr_gate.get_node_type().is_gateset() and q == curr_gate.get_target() and \
                    oldcollection.get_wire_status(q).is_ancilla_state() and gate_type_binary:
                    collection.set_wire_status(q, EdgeType.Basis01)
                elif not curr_gate.get_node_type().is_gateset() and q == curr_gate.get_target() and \
                    oldcollection.get_wire_status(q).is_binary_state() and not gate_type_binary:
                    collection.set_wire_status(q, EdgeType.Basis012)
                elif curr_gate.get_node_type().is_gateset() and gate_type_binary:
                    '''Under optimized currently'''
                    collection.set_wire_status(q, EdgeType.Basis01)
                elif q == curr_gate.get_target():
                    q_pre, q_post = self.gate_slice_collection.get_nodes_before_after_gate(gateslice.get_node(q))
                    if mode == "forward":
                        gateset = q_pre
                        gateset.extend(curr_history)                    
                    else:
                        gateset = q_post
                        gateset = [*curr_history, * gateset]
                    orig_gateset = gateset.copy()
                    # remove old delete pairs
                    for g in self.deleting_pairs.keys():
                        for x in self.deleting_pairs[g]:
                            if x in gateset:
                                gateset.remove(x)
                    '''TODO: Check the get inverse gates for the update for the data for finding the inverse'''
                    if curr_gate in gateset:
                        gateset_delete = []
                        for g in gateset_delete:
                            if g in gateset:
                                gateset.remove(g)
                        self.deleting_pairs[q] = list(set(self.deleting_pairs[q]).union(set(gateset_delete)))
                        gateset = filter_nodes(gateset, q)
                    elif curr_gate not in self.deleting_pairs[q]:
                        raise Exception
                    if len(gateset) == 0:
                        '''Set default mapping for binary case'''
                        collection.set_wire_status(q, EdgeType.Basis01)
                    marked = False
                    for gates in gateset:
                        if gates.get_gate_type().turns_state_to_t():      
                            if mode == "reversed" and collection.get_wire_status(q) != EdgeType.Basis01:
                                collection.set_wire_status(q, EdgeType.Basis012)
                                marked = True
                            elif mode == "forward":
                                collection.set_wire_status(q, EdgeType.Basis012)
                                marked = True
                    if not marked:
                        collection.set_wire_status(q, EdgeType.Basis01)
                else:
                    collection.set_wire_status(q, oldcollection.get_wire_status(q))
        if debug:
            print("updated collection ", collection)
        if origconfig != collection:
            return False
        else:
            return True
        
    def update_propogate(self, idx):
        if idx > 0:
            for i in range(idx, len(self.list_of_edges)):
                gateslice = self.gate_slice_collection.get_slice(i)
                res = self.set_mapping(i)
                if res:
                    break
        return
    
    def get_edge_slice(self, idx):
        return self.list_of_edges[idx]
    
    def node_works_on_qutrit(self, sliceid, node):
        assert sliceid > 0
        if node.get_node_type().is_gateset():
            result = True
            for gates in node.getActiveSet():
                res, gatetype = self.node_works_on_qutrit(sliceid, gates)
                result = result and res
            return result, gatetype
        elif node.get_node_type().is_gate():
            if self.list_of_edges[sliceid].get_wire_status(node.get_target()) == EdgeType.Basis012 and\
                self.list_of_edges[sliceid-1].get_wire_status(node.get_target()) == EdgeType.Basis012:
                return True, GateOpType.MaintainingGate
            elif self.list_of_edges[sliceid].get_wire_status(node.get_target()) == EdgeType.Basis01 and\
                self.list_of_edges[sliceid-1].get_wire_status(node.get_target()) == EdgeType.Basis012:
                return True, GateOpType.CollapsingGate
            elif self.list_of_edges[sliceid].get_wire_status(node.get_target()) == EdgeType.Basis012 and\
                self.list_of_edges[sliceid-1].get_wire_status(node.get_target()) == EdgeType.Basis01:
                return True, GateOpType.LiftingGate
            return False, GateOpType.NoOpGate
