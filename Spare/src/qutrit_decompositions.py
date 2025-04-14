import cirq
import math
import numpy as np

class QutritCZ102Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3, dtype=float)
        flip = np.array([[1., 0., 0.],
                         [0., 1., 0.],
                         [0., 0., -1.]])
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, flip, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)


    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(Z01)]'


class QutritCZ101Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3, dtype=float)
        flip = np.array([[1., 0., 0.],
                         [0., -1., 0.],
                         [0., 0., 1.]])
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, flip, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)


    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(Z01)]'


class QutritCX101Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3, dtype=float)
        flip = np.array([[0., 1., 0.],
                         [1., 0., 0.],
                         [0., 0., 1.]])
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, flip, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1) 

        
    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(X01)]'



class QutritCX102Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        flip = np.array([[0, 0, 1],
                         [0, 1, 0],
                         [1, 0, 0]])
        
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, flip, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1) 

        
    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(X02)]'



class QutritCX112Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        flip = np.array([[1, 0, 0],
                         [0, 0, 1],
                         [0, 1, 0]])
        
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, flip, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1) 

        
    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(X12)]'


class QutritCX201Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        flip = np.array([[0, 1, 0],
                         [1, 0, 0],
                         [0, 0, 1]])
        
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, ones, zero), axis=0), np.concatenate((zero, zero, flip), axis=0)), axis=1) 

        
    def _circuit_diagram_info_(self, args):
        return '[(2)]', '[(X01)]'



class QutritCX202Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        flip = np.array([[0, 0, 1],
                         [0, 1, 0],
                         [1, 0, 0]])
        
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, ones, zero), axis=0), np.concatenate((zero, zero, flip), axis=0)), axis=1) 

        
    def _circuit_diagram_info_(self, args):
        return '[(2)]', '[(X02)]'



class QutritCX212Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        flip = np.array([[1, 0, 0],
                         [0, 0, 1],
                         [0, 1, 0]])
        
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, ones, zero), axis=0), np.concatenate((zero, zero, flip), axis=0)), axis=1) 

        
    def _circuit_diagram_info_(self, args):
        return '[(2)]', '[(X12)]'


class QutritX01Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[0, 1, 0],
                        [1, 0, 0],
                        [0, 0, 1]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(X01)]'


class QutritX02Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[0, 0, 1],
                        [0, 1, 0],
                        [1, 0, 0]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(X02)]'


class QutritX12Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[1, 0, 0],
                        [0, 0, 1],
                        [0, 1, 0]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(X12)]'



class QutritH01Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[1 / math.sqrt(2), 1 / math.sqrt(2), 0],
                        [1 / math.sqrt(2), -1 / math.sqrt(2), 0],
                        [0,                 0,                1]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(H01)]'

class QutritPlusGate(cirq.Gate):
    """A gate that adds one in the computational basis of a qutrit.

    This gate acts on three-level systems. In the computational basis of
    this system it enacts the transformation U|x〉 = |x + 1 mod 3〉, or
    in other words U|0〉 = |1〉, U|1〉 = |2〉, and U|2> = |0〉.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix.
        return np.array([[0, 0, 1],
                         [1, 0, 0],
                         [0, 1, 0]])

    def _circuit_diagram_info_(self, args):
        return '[+1]'

class QutritH02Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[1 / math.sqrt(2), 0, 1/ math.sqrt(2)],
                        [0, 1, 0],
                        [1/ math.sqrt(2), 0, -1/ math.sqrt(2)]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(H02)]'

class QutritMinusGate(cirq.Gate):
    """A gate that adds one in the computational basis of a qutrit.
    This gate acts on three-level systems. In the computational basis of
    this system it enacts the transformation U|x〉 = |x + 1 mod 3〉, or
    in other words U|0〉 = |1〉, U|1〉 = |2〉, and U|2> = |0〉.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        return np.array([[0, 1, 0],
                         [0, 0, 1],
                         [1, 0, 0]])

    def _circuit_diagram_info_(self, args):
        return '[-1]'

class QutritH12Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = 1/math.sqrt(2) * np.array([[math.sqrt(2), 0, 0],
                                         [0, 1 , 1],
                                         [0, 1, -1]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(H12)]'


class QutritT01Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[1, 0, 0],
                        [0, (1 + 1j) / (math.sqrt(2)), 0],
                        [0, 0, 1]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(T01)]'


class QutritT02Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, (1 + 1j)/math.sqrt(2)]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(T02)]'


class QutritT_01Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[1, 0, 0],
                        [0, (1 - 1j) / (math.sqrt(2)), 0],
                        [0, 0, 1]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(T+01)]'


class QutritT_02Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, (1 - 1j)/math.sqrt(2)]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(T+02)]'


class QutritT12Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, (1 + 1j)/math.sqrt(2)]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(T12)]'


class QutritT_12Gate(cirq.Gate):
    """A gate that implements a flip with CX01 gate.
    """
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        arr = np.array([[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, (1 - 1j)/math.sqrt(2)]])
        
        return arr
        
    def _circuit_diagram_info_(self, args):
        return '[(T+12)]'


def qutrits_decomposition(circuit, passes=1):
    for p in range(passes):
        new_circuit = cirq.Circuit()
        for m_index, moment in enumerate(circuit):
            print("__________________")
            print(m_index, " ", moment)
            for op in moment: 
                print(op)
                #print(op._gate, " ", <toffoli.QutritCC1MinusGate>)
                print(cirq.inverse(op))
                # print(cirq.decompose_multi_controlled_x(op)) # decompose(op))
                #print(op.inverse()) #default_decompose()) 

