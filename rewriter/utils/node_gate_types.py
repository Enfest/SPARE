from enum import Enum
import numpy as np

class GateType(Enum):
    PhaseGate = "01"
    X01Gate = "10"
    basis01Gate = "01-basis"
    X12Gate = "12"
    basis12Gate = "12-basis"
    X02Gate = "02"
    basis02Gate = "02-basis"
    PlusGate = "+"
    MinusGate = "-"
    GeneralPermutationGate = "Perm012"
    GeneralGate = "012"
    QubitGate = "qubit_gate"

    def is_qubit_gate(self):
        return self == GateType.QubitGate

    def is_phase_gate(self):
        return self == GateType.PhaseGate

    def is_binary(self):
        return self == GateType.PhaseGate or self == GateType.X01Gate or self == GateType.basis01Gate or self == GateType.QubitGate
    
    def turns_state_to_t(self):
        return not self.is_binary()
        
    def on_01_basis(self):
        return self == GateType.X01Gate or self == GateType.PhaseGate or self == GateType.basis01Gate
    
    def on_02_basis(self):
        return self == GateType.X02Gate or self == GateType.basis01Gate
    
    def on_12_basis(self):
        return self == GateType.X12Gate or self == GateType.basis12Gate
    
    def state_permutation_gates(self):
        return self == GateType.X01Gate or self == GateType.X12Gate or self == GateType.X02Gate \
            or self == GateType.PlusGate or self == self.MinusGate or self == self.GeneralPermutationGate

class NodeType(Enum):
    Qubit = 0
    Gate = 1
    GateSet = 2

    def is_qubit(self):
        return self == NodeType.Qubit
    
    def is_gate(self):
        return self == NodeType.Gate or self == NodeType.GateSet

    def is_gateset(self):
        return self == NodeType.GateSet
    

class GateOpType(Enum):
    LiftingGate = "lifts"
    CollapsingGate = "collapses"
    MaintainingGate = "maintains"
    NoOpGate = "None"

    def normalize_matrix(self, matrix):
        matrix_copy = matrix.copy()
        if matrix.shape[0] == 2:
            return matrix_copy
        if self == GateOpType.LiftingGate:
            matrix_copy[:, 2] = np.abs(matrix[:, 2])
        elif self == GateOpType.CollapsingGate:
            matrix_copy[2, :] = np.abs(matrix[2, :])
        assert np.allclose(matrix_copy @ matrix_copy.conj().T, np.eye(matrix.shape[0])) 
        return matrix_copy


