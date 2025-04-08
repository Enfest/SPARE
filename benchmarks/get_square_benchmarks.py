import cirq
import numpy as np
import square_benchmarks.bench.square_cirq.application.simple as simple
import square_benchmarks.bench.square_cirq.application.errorcorr as errorcorr
import square_benchmarks.bench.square_cirq.application.twooffive as twooffive
import square_benchmarks.bench.square_cirq.application.sixsim as sixsim
import square_benchmarks.bench.square_cirq.application.rd53 as rd53
import square_benchmarks.bench.square_cirq.application.adder as adder


def handle_on_square_benchmark(benchmark, qubits, ancilla):
    def binary_to_ternary(binary_str, debug=False):
        # Convert a binary string to its equivalent ternary string
        ternary_value = int(binary_str, 3)
        if debug:
            print("--", int(binary_str, 2), int(binary_str, 3), np.base_repr(int(binary_str), base=3))
        return ternary_value
    decomposed_circuit = cirq.decompose(benchmark)
    matrix = cirq.unitary(benchmark)
    matrix = cirq.unitary(benchmark)[::2**ancilla, ::2**ancilla]
    print("mat is ", matrix)
    print(matrix)
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
    print("check the norm here ", np.linalg.norm(ternary_matrix @ ternary_matrix.conjugate().T - np.eye(ternary_matrix.shape[0]), 'fro'))
    print(ternary_matrix @ ternary_matrix.conjugate().T)
    assert np.allclose(ternary_matrix @ ternary_matrix.conjugate().T, np.eye(ternary_matrix.shape[0]))
    print("returning ternary matrix ")
    return ternary_matrix

def get_simple_circuits():
    benchmark, qubits, ancilla = simple.main()
    print(benchmark)
    print(qubits)
    print(ancilla)
    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    return ternary_matrix

def get_errorcorr_circuits(return_base=False):
    benchmark, qubits, ancilla = errorcorr.main()
    print(benchmark)
    print(qubits)
    print(ancilla)
    if return_base:
        print(benchmark)
        print(f"qubits {qubits}")
        print(f"ancilla {ancilla}")
        print(type(ancilla))
        return benchmark, ancilla, None

    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    return benchmark, ancilla, ternary_matrix

def get_twooffive_circuits(return_base=False):
    benchmark, qubits, ancilla = twooffive.main()
    print("==========================================")
    print("return base is ", return_base)
    if return_base:
        print(benchmark)
        print(f"qubits {qubits}")
        print(f"ancilla {ancilla}")
        print(type(ancilla))
        return benchmark, ancilla, None

    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    return benchmark, ancilla, ternary_matrix

def get_sixsim_circuits(return_base=False):
    benchmark, qubits, ancilla = sixsim.main()
    print("circuit received")
    if return_base:
        return benchmark, ancilla, None
    ternary_matrix = handle_on_square_benchmark(benchmark, qubits, ancilla)
    print("matrix generated")
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
    print("shape is ", mat.shape[0])
    for i in range(mat.shape[0]):
        assert np.allclose(np.max(mat[i, :]), 1)
