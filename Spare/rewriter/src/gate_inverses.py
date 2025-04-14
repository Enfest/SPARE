import numpy as np
from rewriter.utils.graph_utils import *
from rewriter.src.nodeset import LabelledNode

class DictEntry:
  def __init__(self, node, nodeid):
    assert type(node) == LabelledNode 
    self.matrix_entry = node.get_matrix()
    self.gate_entry = node
    self.gate_id = nodeid
    self.id_name = tuple(tuple(node.get_qubits()), tuple(node.get_controls()))
  
  def __eq__(self, __value: object) -> bool:
    assert type(__value) == LabelledNode 
    return self.id_name == __value.get_id_name()

  def get_id(self):
    return self.id_name
  
  def get_gate_id(self):
    return self.gate_id

class DictEntries:
  def __init__(self, new_entry=None):
    if new_entry is None:
      self.id = None
      self.entries = []
    else:
      assert type(new_entry) == DictEntry
      self.id = new_entry.get_id()
      self.entries = [new_entry]

  def __eq__(self, __value:object) -> bool:
    if self.id is None:
      return False
    return self.id == __value.get_id_name()
  
  def add_entry(self, entry):
    if self.id is None:
      self.id = entry.get_id()
      self.entries.append(entry)
    elif self.id == entry.get_id():
      self.entries.append(entry)
    else:
      raise Exception()
    
  def remove_last_entry(self):
    self.entries.pop(-1)

  def get_last_entry(self):
    assert len(self.entries) > 0
    return self.entries[-1]


class GateInverses:
  def __init__(self, graph):
    self.DictEntryCollection = []
    self.graph = graph

  def check_if_mergable(self, new_entry):
    index = self.DictEntryCollection.index(new_entry)
    intersect = self.graph.print_intersection(new_entry.get_gate_id(), self.DictEntryCollection[index].get_last_gate_id())
    subset = []
    for n in intersect:
      if is_in([self.graph.get_target(n)], self.graph.get_control(new_entry.get_gate_id())):
        if if_mat_preserves_binary(self.graph.get_matrix(n)):
          pass
        if (self.graph.check_gate_type(n) == "01" or self.graph.check_gate_type(n) == "10" or self.graph.check_gate_type(n) == "01-") and \
          self.graph.get_control_value(nodeid, self.graph.get_target(n)) == 2:
          pass
        if (self.graph.check_gate_type(n) == "12" or self.graph.check_gate_type(n) == "12-") and \
          self.graph.get_control_value(nodeid, self.graph.get_target(n)) == 0:
          pass
        if (self.graph.check_gate_type(n) == "02" or self.graph.check_gate_type(n) == "02-") and \
          self.graph.get_control_value(nodeid, self.graph.get_target(n)) == 1:
          pass
        else:
          subset.append(n)
    
    ginv_small = GateInverses(self.graph)
    for n_ in subset:
      ginv_small.add_node(n_)
    if ginv_small.if_free():
      return True
    else:
      return False

    return True

  def add_node(self, nodeid, node):
    new_entry = DictEntry(node, nodeid)
    if new_entry not in self.DictEntryCollection:
      self.DictEntryCollection.append(DictEntries(new_entry))
    else:
      res = self.check_if_mergable(nodeid, node, new_entry)
      index = self.DictEntryCollection.index(new_entry)
      if res:
        self.DictEntryCollection[index].update(new_entry)
      else:
        self.DictEntryCollection[index].add_entry(new_entry)
    return

  def get_status(self, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    if if_mat_preserves_binary(self.matrix_collection[arr]):
      return True
    return False
  
  def get_all_status(self, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    if if_mat_preserves_binary(self.matrix_collection[arr]):
      return True
    if if_mat_preserves_binary(self.matrix_collection[arr] @ cirq.unitary(QutritX01Gate())):
      return True
    if np.allclose(self.matrix_collection[arr]["matrix"][2, 2], 1):
      return True 
    return False

  def reset_node(self, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    return self.all_gates_map[arr]
  
  def hard_reset(self, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    del self.gate_collection[arr]
    del self.matrix_collection[arr]
    del self.all_gates_map[arr]
  
  def get_qubits(self, node):
    return node.get_qubits(), node.get_control_values()
  
  def exists(self, debug=False):
    assert len(self.all_gates_map.keys()) == len(self.matrix_collection.keys())
    assert len(self.gate_collection.keys()) == len(self.matrix_collection.keys())
    assert len(self.gate_collection.keys()) == len(self.all_gates_map.keys())
    if debug:
      print(self.all_gates_map)
      print(self.gate_collection)
      print(self.gate_collection)

    if len(self.matrix_collection.keys()) > 0 and len(self.gate_collection.keys()) > 0 and len(self.all_gates_map.keys()) > 0:
      return True
    for a in self.matrix_collection.keys():
      if np.allclose(self.matrix_collection[a], np.eye(3)):
        del self.matrix_collection[a]
        del self.gate_collection[a]
        del self.all_gates_map[a]
    if len(self.matrix_collection.keys()) > 0 and len(self.gate_collection.keys()) > 0 and len(self.all_gates_map.keys()) > 0:
      if debug:
        print("dict exists")
      return True
    if debug:
      print("doesnt exist")
    return False

  def remove(self, e, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    if arr in self.all_gates_map:
      # print("AA ", self.all_gates_map[arr])
      if e in self.all_gates_map[arr]:
        self.all_gates_map[arr].remove(e)
        # print("BB ", self.all_gates_map[arr])
        if len(self.all_gates_map[arr]) == 0:
          del self.all_gates_map[arr]
          del self.matrix_collection[arr]
          del self.gate_collection[arr]
          

class TrackingDict:
  def __init__(self, graph):
    # self.DictEntryCollection = 
    self.gate_collection = {}
    self.matrix_collection = {}
    self.all_gates_map = {}
    self.graph = graph

  def check_if_cancellable(self, nodeid, node):
    qubits = node.get_qubits()
    controls = node.get_control_values()
    print(node)
    print(qubits, controls)
    qubits, controls = self.get_qubits(node)
    elem = self.all_gates_map[tuple([tuple(qubits), tuple(controls)])][-1]
    intersect = self.graph.print_intersection(elem, nodeid)
    for n in intersect:
      assert isinstance(n, str)
      assert isinstance(nodeid, str)
      if is_in([self.graph.get_target(n)], self.graph.get_control(nodeid)):
        if if_mat_preserves_binary(self.graph.get_matrix(n)):
          pass
        if (self.graph.check_gate_type(n) == "01" or self.graph.check_gate_type(n) == "10" or self.graph.check_gate_type(n) == "01-") and \
          self.graph.get_control_value(nodeid, self.graph.get_target(n)) == 2:
          pass
        if (self.graph.check_gate_type(n) == "12" or self.graph.check_gate_type(n) == "12-") and \
          self.graph.get_control_value(nodeid, self.graph.get_target(n)) == 0:
          pass
        if (self.graph.check_gate_type(n) == "02" or self.graph.check_gate_type(n) == "02-") and \
          self.graph.get_control_value(nodeid, self.graph.get_target(n)) == 1:
          pass
        else:
          return False
    return True

  def add_node(self, nodeid, node):
    # new_entry = DictEntry(node, nodeid)
    qubits, controls = self.get_qubits(node)
    controls = node.get_control_values()
    # if new_entry not in self.DictEntryCollection:
    #   self.DictEntryCollection.append(new_entry)
    # else:
    #   res = self.check_if_cancellable(new_entry)
    print("q and c ", qubits, controls, self.get_qubits(node))
    if tuple([tuple(qubits), tuple(controls)]) not in self.gate_collection.keys():
      self.gate_collection[tuple([tuple(qubits), tuple(controls)])] = node
      self.matrix_collection[tuple([tuple(qubits), tuple(controls)])] = node.get_matrix()
      self.all_gates_map[tuple([tuple(qubits), tuple(controls)])] = [nodeid]
    else:
      res = self.check_if_cancellable(nodeid, node)
      if res:
        self.matrix_collection[tuple([tuple(qubits), tuple(controls)])] = node.get_matrix() @ \
                                                                          self.matrix_collection[tuple([tuple(qubits), tuple(controls)])] 
        self.all_gates_map[tuple([tuple(qubits), tuple(controls)])].append(nodeid)
      else:
        # We cant support this node should not check this further for now
        return False
    return

  def get_status(self, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    if if_mat_preserves_binary(self.matrix_collection[arr]):
      return True
    return False
  
  def get_all_status(self, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    if if_mat_preserves_binary(self.matrix_collection[arr]):
      return True
    if if_mat_preserves_binary(self.matrix_collection[arr] @ cirq.unitary(QutritX01Gate())):
      return True
    if np.allclose(self.matrix_collection[arr]["matrix"][2, 2], 1):
      return True 
    return False

  def reset_node(self, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    return self.all_gates_map[arr]
  
  def hard_reset(self, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    del self.gate_collection[arr]
    del self.matrix_collection[arr]
    del self.all_gates_map[arr]
  
  def get_qubits(self, node):
    return node.get_qubits(), node.get_control_values()
  
  def exists(self, debug=False):
    assert len(self.all_gates_map.keys()) == len(self.matrix_collection.keys())
    assert len(self.gate_collection.keys()) == len(self.matrix_collection.keys())
    assert len(self.gate_collection.keys()) == len(self.all_gates_map.keys())
    if debug:
      print(self.all_gates_map)
      print(self.gate_collection)
      print(self.gate_collection)

    if len(self.matrix_collection.keys()) > 0 and len(self.gate_collection.keys()) > 0 and len(self.all_gates_map.keys()) > 0:
      return True
    for a in self.matrix_collection.keys():
      if np.allclose(self.matrix_collection[a], np.eye(3)):
        del self.matrix_collection[a]
        del self.gate_collection[a]
        del self.all_gates_map[a]
    if len(self.matrix_collection.keys()) > 0 and len(self.gate_collection.keys()) > 0 and len(self.all_gates_map.keys()) > 0:
      if debug:
        print("dict exists")
      return True
    if debug:
      print("doesnt exist")
    return False

  def remove(self, e, node):
    qubits, controls = self.get_qubits(node)
    arr = tuple([tuple(qubits), tuple(controls)])
    if arr in self.all_gates_map:
      # print("AA ", self.all_gates_map[arr])
      if e in self.all_gates_map[arr]:
        self.all_gates_map[arr].remove(e)
        # print("BB ", self.all_gates_map[arr])
        if len(self.all_gates_map[arr]) == 0:
          del self.all_gates_map[arr]
          del self.matrix_collection[arr]
          del self.gate_collection[arr]
