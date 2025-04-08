""" 2of5.py - Cirq implementation of the 2of5 function.

=== Description ===
2of5 function takes 5 inputs and 1 output. The output is 1 if and only if
exactly two of its inputs are 1.

=== Reference ===
Adapted from the Reversible Logic Synthesis Benchmarks Page:
http://webhome.cs.uvic.ca/~dmaslov/

Part of the SQUARE benchmark set: https://arxiv.org/abs/2004.08539

=== Note ===
This implementation performs uncomputation.

"""

import numpy as np
import random
import cirq
from cirq.ops import raw_types
from src.cirq_depth_counter import depth_counter
from src.toffoli import *

# Get the original Toffoli gate class from Cirq
original_toffoli = cirq.TOFFOLI

# Define a new decomposition method that prevents decomposition
def no_decomposition(self):
    # Prevent decomposition by returning NotImplemente
    return NotImplemented

# Dynamically override the _decompose_ method for the Cirq Toffoli gate

cirq.TOFFOLI._decompose_ = no_decomposition
cirq.CNOT._decompose_ = no_decomposition

class Compute(raw_types.Gate):
    """
    """
    def __init__(self, num_in, num_out, num_anc, is_inverse=False):
        self.num_in = num_in
        self.num_out = num_out
        self.num_anc = num_anc
        self.is_inverse = is_inverse
        self._num_qubits = num_in + num_out + num_anc

    def num_qubits(self):
        return self._num_qubits

    def __pow__(self, power):
        if power == -1:
            return Compute(self.num_in, self.num_out, self.num_anc, is_inverse=True)
        else:
            return NotImplemented
        
    def _decompose_(self, qubits):
        in_qubits = qubits[:self.num_in]
        out_qubits = qubits[self.num_in:self.num_in+self.num_out]
        anc_qubits = qubits[self.num_in+self.num_out:]
        self.data = [cirq.TOFFOLI(in_qubits[0], in_qubits[1], anc_qubits[0]),
        cirq.CNOT(in_qubits[0], in_qubits[1]),
        cirq.TOFFOLI(in_qubits[2], anc_qubits[0], anc_qubits[1]),
        cirq.TOFFOLI(in_qubits[1], in_qubits[2], anc_qubits[0]),
        cirq.CNOT(in_qubits[1], in_qubits[2]),
        cirq.TOFFOLI(in_qubits[3], anc_qubits[0], anc_qubits[1]),
        cirq.TOFFOLI(in_qubits[2], in_qubits[3], anc_qubits[0]),
        cirq.CNOT(in_qubits[2], in_qubits[3]),
        cirq.TOFFOLI(in_qubits[4], anc_qubits[0], anc_qubits[1]),
        cirq.TOFFOLI(in_qubits[3], in_qubits[4], anc_qubits[0]),
        cirq.CNOT(in_qubits[3], in_qubits[4]),
        cirq.CNOT(anc_qubits[1], anc_qubits[0])]
        if self.is_inverse:
            for i in range(len(self.data)-1, -1, -1):
                yield self.data[i]
        else:
            for i in range(len(self.data)):
                yield self.data[i]

def print_circ():
    n_iter = 1
    num_in = 5
    num_out = 1
    num_anc = 2 # per iter
    c = cirq.Circuit()
    in_qubits = [cirq.GridQubit(i,0) for i in range(num_in)]
    out_qubits = [cirq.GridQubit(i+num_in,0) for i in range(num_out)]
    j = 0
    anc_qubits = [cirq.GridQubit(i+j*num_anc+num_in+num_out,0) for i in range(num_anc)]
    all_qubits = in_qubits + out_qubits + anc_qubits
    compute = Compute(num_in, num_out, num_anc)
    circ = compute(*all_qubits) # forward
    c.append([circ])
    c2 = cirq.Circuit()
    c2.append([cirq.decompose(circ)])
    print(c2)
        

def main():
    n_iter = 1
    num_in = 5
    num_out = 1
    num_anc = 2 # per iter
    #num_qubit = num_in + num_out + num_anc
    # Initialize qubits and circuit
    c = cirq.Circuit()
    new_circ = cirq.Circuit()
    compute = Compute(num_in, num_out, num_anc)

    cntr = 0
    for j in range(n_iter):
        in_qubits = [cirq.GridQubit(i,0) for i in range(num_in)]
        out_qubits = [cirq.GridQubit(i+num_in,0) for i in range(num_out)]
        # Pick a random input
        rand_bits = [random.randint(0,1) for i in range(num_in)]
        anc_qubits = [cirq.GridQubit(i+j*num_anc+num_in+num_out,0) for i in range(num_anc)]
        all_qubits = in_qubits + out_qubits + anc_qubits
        circ = compute(*all_qubits)
        two_circ = cirq.Circuit()
        two_circ.append(circ)
        c.append([circ])
        c.append([cirq.CNOT(anc_qubits[0], out_qubits[0])])
        c.append([cirq.inverse(circ)])
        new_circ.append(cirq.decompose(c))

    mat = cirq.unitary(c)
    closest_to_1_indices = np.argmax(np.abs(mat), axis=1)
    mat = cirq.unitary(c)[::2**num_anc, ::2**num_anc]
    assert np.allclose(mat @ mat.T, np.eye(mat.shape[0]))
    closest_to_1_indices = np.argmax(np.abs(mat), axis=1)
    print(closest_to_1_indices)
    print("new circ here")
    print(new_circ)
    print(all_qubits)
    print(depth_counter(new_circ, all_qubits))
    return new_circ, num_in+num_out, num_anc

def encoding_circuit():
    c = cirq.Circuit()
    qubits = cirq.LineQid.range(2, dimension=3)
    c.append(QutritPlusGate().on(qubits[0]).controlled_by(qubits[1]))
    c.append(QutritX12Gate().on(qubits[0]))
    c.append(QutritX01Gate().on(qubits[1]).controlled_by(qubits[0]))
    c.append(QutritX12Gate().on(qubits[0]))
    print(c)
    print(cirq.unitary(c))

def make_unitary():
    mat = np.zeros((27, 27))
    vals_set = []
    arr = []
    for i in range(8):
        in_binary = bin(i)[2:]
        string = []
        i_copy = i
        v = 3
        while v:
            string.append(i_copy % 3)
            i_copy = i_copy // 3
            v = v - 1
        string.reverse()
        ternary_i = sum([int(j_val) * 3**(j) for j, j_val in enumerate(reversed(in_binary))])
        arr.append(ternary_i)
    print(arr)
    for i, j in enumerate(arr):
        # if i != ternary_i:
        print("i and j - ", i, j)

        mat[i, j] = 1
        if i not in arr and i > 8:
            print("j and i", j, i)
            if j not in vals_set:
                vals_set.append(j)
            mat[j, i] = 1
        if i not in vals_set:
            vals_set.append(i)

    print(vals_set)
    print(list(set(range(mat.shape[0])).difference(vals_set)))
    for i in list(set(range(mat.shape[0])).difference(vals_set)):
        mat[i, i] = 1
    # print(mat) 
    assert np.allclose(mat @ mat.T, np.eye(mat.shape[0]))
    print(mat)


if __name__ == '__main__':
    # main()
    # encoding_circuit()
    print_circ()
