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
import math
from cirq.ops import raw_types

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
        self.data = [
        cirq.CNOT(in_qubits[2], anc_qubits[3]),
        cirq.CNOT(in_qubits[1], anc_qubits[0]),
        cirq.CNOT(in_qubits[1], anc_qubits[1]),
        cirq.CNOT(in_qubits[1], anc_qubits[2]),
        cirq.CNOT(in_qubits[1], anc_qubits[3]),
        cirq.CNOT(in_qubits[0], anc_qubits[1]),
        cirq.CNOT(in_qubits[0], anc_qubits[2]),
        cirq.CNOT(in_qubits[0], anc_qubits[3])
        ]
        if self.is_inverse:
            for i in range(len(self.data)-1, -1, -1):
                yield self.data[i]
        else:
            for i in range(len(self.data)):
                yield self.data[i]

def main():
    n_iter = 1
    num_in = 3
    num_out = 0
    num_anc = 4 # per iter
    #num_qubit = num_in + num_out + num_anc

    # Initialize qubits and circuit
    c = cirq.Circuit()
    compute = Compute(num_in, num_out, num_anc)

    new_circ = cirq.Circuit()
    cntr = 0
    for j in range(n_iter):
        in_qubits = [cirq.GridQubit(i,0) for i in range(num_in)]
        # anc_qubits = [cirq.GridQubit(i+num_in,0) for i in range(num_anc)]
        # out_qubits = [cirq.GridQubit(i+num_in+num_anc,0) for i in range(num_out)]
        out_qubits = [cirq.GridQubit(i+num_in,0) for i in range(num_out)]

        # Pick a random input
        rand_bits = [random.randint(0,1) for i in range(num_in)]
        # for r in range(num_in):
        #    if rand_bits[r] == 1:
        #        c.append([cirq.X(in_qubits[r])])

        anc_qubits = [cirq.GridQubit(i+j*num_anc+num_in+num_out,0) for i in range(num_anc)]
        all_qubits = in_qubits + out_qubits + anc_qubits
        circ = compute(*all_qubits) # forward
        print(circ)
        c.append([circ])
        for i in range(num_in):
            c.append([cirq.Rz(rads=math.pi/8)(in_qubits[i])])
        for i in range(num_anc):
            c.append([cirq.Rz(rads=math.pi/8)(anc_qubits[i])])
        # c.append([cirq.CNOT(anc_qubits[0], out_qubits[0])]) # copy
        c.append([cirq.inverse(circ)]) # backward
        # c.append([cirq.measure(*out_qubits, key='result')])
        new_circ.append(cirq.decompose(c))

    print(new_circ)
    mat = cirq.unitary(c)
    closest_to_1_indices = np.argmax(np.abs(mat), axis=1)
    print("old matrix")
    print(mat)
    mat = cirq.unitary(c)[::2**num_anc, ::2**num_anc]
    print(mat)
    print(mat @ mat.conjugate().T)
    assert np.allclose(mat @ mat.conjugate().T, np.eye(mat.shape[0]))
    closest_to_1_indices = np.argmax(np.abs(mat), axis=1)
    print(closest_to_1_indices)
    simulator = cirq.Simulator()
    # result = simulator.run(c, repetitions=2)
    # print("Results:")
    # print(rand_bits, out_qubits, * out_qubits)
    # print(result)
    return new_circ, num_in+num_out, num_anc

if __name__ == '__main__':
    main()
