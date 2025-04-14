import numpy as np
import math
from scipy.stats import unitary_group


def select_rows_and_columns(matrix, dimension, considered_qubits):
    """
    Select rows and columns of a matrix that are either 0 or 1 in ternary format.
    
    Parameters:
    matrix (numpy.ndarray): Input matrix of size n x n.
    
    Returns:
    numpy.ndarray: Output sub-matrix of size m x m, where m <= n.
    """
    # convert row and column indexes to ternary format
    # Convert positions to a set for faster lookup
    positions_set = set(considered_qubits)
    # Define a function to check if a number can be represented using only '0' and '1' in the given base
    def is_binary_in_base(n, base, positions_set, size):
        size = int(size)
        repr_str = np.base_repr(n, base=base)
        if positions_set is None:
            return np.isin(list(np.base_repr(n, base=base)), ['0', '1']).all()
        else:
            repr_str = repr_str.zfill(size)  # Ensure it's long enough
            return all(repr_str[i] in ['0', '1'] for i in positions_set) and all(repr_str[i] == '0' for i in range(len(repr_str)) if i not in positions_set)

    # Compute ternary_rows and ternary_cols
    ternary_rows = np.array([
        is_binary_in_base(i, base=dimension, positions_set=positions_set, size=math.log(matrix.shape[0], dimension)) for i in range(matrix.shape[0])
    ])
    ternary_cols = np.array([
        is_binary_in_base(i, base=dimension, positions_set=positions_set, size=math.log(matrix.shape[0], dimension)) for i in range(matrix.shape[1])  # Fixed from shape[0] to shape[1] for columns
    ])    
    # ternary_rows = np.array([np.isin(list(np.base_repr(i, base=dimension)), ['0', '1']).all() for i in range(matrix.shape[0])])
    # ternary_cols = np.array([np.isin(list(np.base_repr(i, base=dimension)), ['0', '1']).all() for i in range(matrix.shape[0])])
    # select valid rows and columns
    submatrix = matrix[ternary_rows]
    submatrix = submatrix[:, ternary_cols]
    return submatrix

def check_circuit_equality(unitary1, unitary2, dimension=3, qubits_considered=None, dtype=np.complex128):
    size = 2
    qubits_considered = [x.get_qubits()[0] for x in qubits_considered]
    submat1 = select_rows_and_columns(unitary1, dimension, considered_qubits=qubits_considered)
    submat2 = select_rows_and_columns(unitary2, dimension, considered_qubits=qubits_considered)
    norm = np.sqrt(np.sum(np.abs(submat1)**2))
    submat_ = submat1 / norm
    if not np.allclose(submat1, submat2):
        for i in range(submat1.shape[0]):
            if not np.allclose(submat1[i, :], submat2[i, :]):
                print("orig ", i, submat1[i, :])
                print("newe ", i, submat2[i, :])
    assert np.allclose(submat1, submat2)
    return submat1, submat2

def return_qubit_circuit(matrix):
    return select_rows_and_columns(matrix)

if __name__ == "__main__":
    dimension = 3
    size = 2
    unitary1 = unitary_group.rvs(dimension**size)
    submat1, submat_ = check_circuit_equality(unitary1, unitary1)
    print(unitary1)
    print(submat1)
    print(submat_ @ submat_.conj().T)
