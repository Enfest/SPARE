
from enum import Enum
import re
import numpy as np

class EdgeType(Enum):
    Basis0 = "0"
    Basis1 = "1"
    Basis2 = "2"
    Basis01 = "01"
    Basis02 = "02"
    Basis12 = "12"
    Basis012 = "012"

    def is_0_basis(self):
        return self == EdgeType.Basis0

    def is_1_basis(self):
        return self == EdgeType.Basis1

    def is_01_basis(self):
        return self == EdgeType.Basis01
    
    def is_02_basis(self):
        return self == EdgeType.Basis02
    
    def is_12_basis(self):
        return self == EdgeType.Basis12
    
    def is_012_basis(self):
        return self == EdgeType.Basis012
    
    def is_ancilla_state(self):
        return self.is_0_basis() or self.is_1_basis()
    
    def is_binary_state(self):
        return self.is_0_basis() or self.is_1_basis() or self.is_01_basis()
    
    def is_equal(self, state):
        states_possible = re.findall(r'[0-9]+', str(self))
        # print("EQUAL CHECK-", self, str(state), str(states_possible[-1]))
        return str(state) == str(states_possible[-1])

    def contains(self, state):
        states_possible = re.findall(r'[0-9]+', str(self))
        return str(state) in str(states_possible[-1])
    
    def decimate_gate(self, matrix):
        for i in range(matrix.shape[0]):
            assert np.allclose(matrix[i, i] * matrix[i, i].conjugate(), 1)
        if self == EdgeType.Basis0:
            if matrix.shape[0] == 2:
                return np.allclose(matrix[1, 1], 1)
            return np.allclose(matrix[1, 1], 1) and np.allclose(matrix[2, 2], 1)
        if self == EdgeType.Basis1:
            if matrix.shape[0] == 2:
                return np.allclose(matrix[0, 0], 1)
            return np.allclose(matrix[0, 0], 1) and np.allclose(matrix[2, 2], 1)
        if self == EdgeType.Basis2 and matrix.shape[0] > 2:
            return np.allclose(matrix[0, 0], 1) and np.allclose(matrix[1, 1], 1)
        if self == EdgeType.Basis01 and matrix.shape[0] > 2:
            return np.allclose(matrix[2, 2], 1)
        if self == EdgeType.Basis02 and matrix.shape[0] > 2:
            return np.allclose(matrix[1, 1], 1)
        if self == EdgeType.Basis12 and matrix.shape[0] > 2:
            return np.allclose(matrix[0, 0], 1)
        return False

    @classmethod
    def from_string(cls, enum_str):
        enum_str = enum_str.replace("EdgeType.", "")
        for member in cls:
            if member.name == enum_str:
                return member
        raise ValueError(f"No such enum value: {enum_str}")