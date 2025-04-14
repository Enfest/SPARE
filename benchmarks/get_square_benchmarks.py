import cirq
import numpy as np
from .square_benchmarks.bench.square_cirq.application import simple, errorcorr, twooffive, sixsim, rd53, adder

def handle_on_square_benchmark(benchmark, qubits, ancilla):
    def binary_to_ternary(binary_str, debug=False):
        # Convert a binary string to its equivalent ternary string
        ternary_value = int(binary_str, 3)
        if debug:
            print("--", int(binary_str, 2), int(binary_str, 3), np.base_repr(int(binary_str), base=3))
        return ternary_value
    matrix = cirq.unitary(benchmark)[::2**ancilla, ::2**ancilla]
    ternary_matrix = np.zeros((3**qubits, 3**qubits), dtype=complex)
    other_indexes = list(range(ternary_matrix.shape[0]))
    n = qubits
    maps = {}
    for i in range(matrix.shape[0]):
        num = binary_to_ternary(format(i, f'0{n}b'))
        arr = []
        while num > 0:
            arr.append(num % 3)
            num = num // 3
        assert 2 not in arr
        for j in range(matrix.shape[1]):
            ternary_i = binary_to_ternary(format(i, f'0{n}b'))
            ternary_j = binary_to_ternary(format(j, f'0{n}b'))
            ternary_matrix[ternary_i, ternary_j] = matrix[i, j]
            maps[ternary_i] = ternary_j
            # if matrix[i, j] > 0.8:
            #     print(ternary_i, ternary_j, matrix[i, j])
        other_indexes.remove(ternary_i)
    for i in other_indexes:
        ternary_matrix[i, i] = 1
    assert np.allclose(ternary_matrix @ ternary_matrix.conjugate().T, np.eye(ternary_matrix.shape[0]))
    return ternary_matrix

def get_simple_circuits():
    benchmark, qubits, ancilla = simple.main()
    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    return ternary_matrix

def get_errorcorr_circuits(return_base=False):
    benchmark, qubits, ancilla = errorcorr.main()
    if return_base:
        return benchmark, ancilla, None
    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    return benchmark, ancilla, ternary_matrix

def get_twooffive_circuits(return_base=False):
    benchmark, qubits, ancilla = twooffive.main()
    if return_base:
        return benchmark, ancilla, None
    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    return benchmark, ancilla, ternary_matrix

def get_sixsim_circuits(return_base=False):
    benchmark, qubits, ancilla = sixsim.main()
    if return_base:
        return benchmark, ancilla, None
    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    return benchmark, ancilla, ternary_matrix

def get_rd53_circuits(return_base=False):
    benchmark, qubits, ancilla = rd53.main()
    if return_base:
        return benchmark, ancilla, None
    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    return ternary_matrix

def get_adder_circuits():
    benchmark, qubits, ancilla = adder.main()
    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    return ternary_matrix

if __name__ == '__main__':
    mat = get_sixsim_circuits()
    for i in range(mat.shape[0]):
        assert np.allclose(np.max(mat[i, :]), 1)
