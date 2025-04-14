import cirq
import numpy as np
from src.qutrit_decompositions import *

class QutritCC1MinusGate(cirq.Gate):
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
        return (3,3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 1, 0],
                         [0, 0, 1],
                         [1, 0, 0]])
        
        controlledC = np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, plus, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)
        return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), controlledC, np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)), axis=0)), axis=1)

    
    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(1)]', '[-1]'


class QutritControlPlusGate(cirq.Gate):
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (3,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a single qutrit.
        return (3, 3,)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix.
        return np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1]])

    def _circuit_diagram_info_(self, args):
        return '[+1]'




class QutritC1MinusGate(cirq.Gate):
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
        return (3, 3)

    def _num_qubits_(self):
        return 2

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 1, 0],
                         [0, 0, 1],
                         [1, 0, 0]])
        
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, plus, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)
    
    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[-1]'


class QutritCC2MinusGate(cirq.Gate):
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
        return (3,3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 1, 0],
                         [0, 0, 1],
                         [1, 0, 0]])
        
        controlledC = np.array([[ones, zero, zero],
                                [zero, ones, zero],
                                [zero, zero, plus]])
        # return np.array([[np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), np.eye(9), np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), np.zeros((9, 9)), controlledC]], dtype=object)
        controlledC = np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, ones, zero), axis=0), np.concatenate((zero, zero, plus), axis=0)), axis=1)
        return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.eye(9), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), controlledC), axis=0)), axis=1)

   
    def _circuit_diagram_info_(self, args):
        return '[(2)]', '[(2)]', '[-1]'


class QutritC2MinusGate(cirq.Gate):
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
        return (3, 3)

    def _num_qubits_(self):
        return 2

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 1, 0],
                         [0, 0, 1],
                         [1, 0, 0]])
        #return np.array([[ones, zero, zero],
        #                 [zero, ones, zero],
        #                 [zero, zero, plus]], dtype=object)
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, ones, zero), axis=0), np.concatenate((zero, zero, plus), axis=0)), axis=1)
        # return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0), np.concatenate((np.zeros((9, 9)), controlledC, np.zeros((9, 9))), axis=0), np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)), axis=0)), axis=1)


    def _circuit_diagram_info_(self, args):
        return '[(2)]', '[-1]'
 

class QutritCC1PlusGate(cirq.Gate):
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
        return (3,3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 0, 1],
                         [1, 0, 0],
                         [0, 1, 0]])
        
        controlledC = np.array([[ones, zero, zero],
                                [zero, plus, zero],
                                [zero, zero, ones]])
        #return np.array([[np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), controlledC, np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)]], dtype=object)
        controlledC = np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, plus, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)
        return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), controlledC, np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)), axis=0)), axis=1)



    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(1)]', '[+1]'


class QutritCC1Flip01Gate(cirq.Gate):
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
        return (3,3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 1, 0],
                         [1, 0, 0],
                         [0, 0, 1]])
        
        controlledC = np.array([[ones, zero, zero],
                                [zero, plus, zero],
                                [zero, zero, ones]])
        #return np.array([[np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), controlledC, np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)]], dtype=object)
        controlledC = np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, plus, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)
        return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), controlledC, np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)), axis=0)), axis=1)

    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(1)]', '[X01]'


class QutritCC1Flip12Gate(cirq.Gate):
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
        return (3,3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[1, 0, 0],
                         [0, 0, 1],
                         [0, 1, 0]])
        
        controlledC = np.array([[ones, zero, zero],
                                [zero, plus, zero],
                                [zero, zero, ones]])
        #return np.array([[np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), controlledC, np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)]], dtype=object)
        controlledC = np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, plus, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)
        return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), controlledC, np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)), axis=0)), axis=1)

    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(1)]', '[X12]'


class QutritCC1Flip02Gate(cirq.Gate):
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
        return (3,3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 0, 1],
                         [0, 1, 0],
                         [1, 0, 0]])
        
        controlledC = np.array([[ones, zero, zero],
                                [zero, plus, zero],
                                [zero, zero, ones]])
        #return np.array([[np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), controlledC, np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)]], dtype=object)
        controlledC = np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, plus, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)
        return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), controlledC, np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)), axis=0)), axis=1)

    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(1)]', '[X02]'


class QutritCC2Flip02Gate(cirq.Gate):
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
        return (3,3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 0, 1],
                         [0, 1, 0],
                         [1, 0, 0]])
        
        controlledC = np.array([[ones, zero, zero],
                                [zero, ones, zero],
                                [zero, zero, plus]])
        #return np.array([[np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), controlledC, np.zeros((9, 9))],
        #                 [np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)]], dtype=object)
        controlledC = np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, plus, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)
        return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.eye(9) , np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), controlledC), axis=0)), axis=1)

    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(1)]', '[X02]'


class QutritC1PlusGate(cirq.Gate):
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
        return (3, 3)

    def _num_qubits_(self):
        return 2

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 0, 1],
                         [1, 0, 0],
                         [0, 1, 0]])
        #return np.array([[ones, zero, zero],
        #                 [zero, plus, zero],
        #                 [zero, zero, ones]], dtype=object)
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, plus, zero), axis=0), np.concatenate((zero, zero, ones), axis=0)), axis=1)
        # return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0), np.concatenate((np.zeros((9, 9)), controlledC, np.zeros((9, 9))), axis=0), np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)), axis=0)), axis=1)

    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[+1]'


class QutritCC2PlusGate(cirq.Gate):
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
        return (3,3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 0, 1],
                         [1, 0, 0],
                         [0, 1, 0]])
        
        controlledC = np.array([[ones, zero, zero],
                                [zero, ones, zero],
                                [zero, zero, plus]])
        # return np.array([[np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))],
        #                  [np.zeros((9, 9)), np.eye(9), np.zeros((9, 9))],
        #                  [np.zeros((9, 9)), np.zeros((9, 9)), controlledC]], dtype=object)
        controlledC = np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, ones, zero), axis=0), np.concatenate((zero, zero, plus), axis=0)), axis=1)
        return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.eye(9), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), controlledC), axis=0)), axis=1)
    
    def _circuit_diagram_info_(self, args):
        return '[(2)]', '[(2)]', '[+1]'


class QutritCC1X01Gate(cirq.Gate):
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
        return (3,3,3)

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 1, 0],
                         [1, 0, 0],
                         [0, 0, 1]])
        
        controlledC = np.array([[ones, zero, zero],
                                [zero, plus, zero],
                                [zero, zero, ones]])
        # return np.array([[np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))],
        #                  [np.zeros((9, 9)), np.eye(9), np.zeros((9, 9))],
        #                  [np.zeros((9, 9)), np.zeros((9, 9)), controlledC]], dtype=object)
        controlledC = np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, ones, zero), axis=0), np.concatenate((zero, zero, plus), axis=0)), axis=1)
        return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), controlledC, np.zeros((9, 9))), axis=0),
                               np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)), axis=0)), axis=1)
    
    def _circuit_diagram_info_(self, args):
        return '[(1)]', '[(1)]', '[X01]'


class QutritC2PlusGate(cirq.Gate):
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
        return (3, 3)

    def _num_qubits_(self):
        return 2

    def _unitary_(self):
        # Since the gate acts on three level systems it has a unitary 
        # effect which is a three by three unitary matrix
        zero = np.zeros((3, 3))
        ones = np.eye(3)
        plus = np.array([[0, 0, 1],
                         [1, 0, 0],
                         [0, 1, 0]])
        # return np.array([[ones, zero, zero],
        #                  [zero, ones, zero],
        #                  [zero, zero, plus]], dtype=object)
        return np.concatenate((np.concatenate((ones, zero, zero), axis=0), np.concatenate((zero, ones, zero), axis=0), np.concatenate((zero, zero, plus), axis=0)), axis=1)
        # return np.concatenate((np.concatenate((np.eye(9), np.zeros((9, 9)), np.zeros((9, 9))), axis=0), np.concatenate((np.zeros((9, 9)), controlledC, np.zeros((9, 9))), axis=0), np.concatenate((np.zeros((9, 9)), np.zeros((9, 9)), np.eye(9)), axis=0)), axis=1)



    def _circuit_diagram_info_(self, args):
        return '[(2)]', '[+1]'


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
        # effect which is a three by three unitary matrix
        return np.array([[0, 0, 1],
                         [1, 0, 0],
                         [0, 1, 0]])

    def _circuit_diagram_info_(self, args):
        return '[+1]'

def generate_toffoli(qutrits, i, j, k, num):
    circuit = cirq.Circuit()
    if num == "101" or num == "110":
        print("PRINT CASE VAL", num)
        circuit.append(QutritH01Gate().on(qutrits[k]))   # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX101Gate().on(qutrits[j], qutrits[k]))

        circuit.append(QutritT_01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX101Gate().on(qutrits[i], qutrits[k]))

        circuit.append(QutritT01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX101Gate().on(qutrits[j], qutrits[k]))
        circuit.append(QutritT_01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX101Gate().on(qutrits[i], qutrits[k]))
        circuit.append(QutritT01Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX101Gate().on(qutrits[i], qutrits[j]))
        circuit.append(QutritH01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritT01Gate().on(qutrits[i]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT_01Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX101Gate().on(qutrits[i], qutrits[j]))
    if num == "102" or num == "120":

        circuit.append(QutritH02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX102Gate().on(qutrits[j], qutrits[k]))

        circuit.append(QutritT_02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX102Gate().on(qutrits[i], qutrits[k]))

        circuit.append(QutritT02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX102Gate().on(qutrits[j], qutrits[k]))
        circuit.append(QutritT_02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX102Gate().on(qutrits[i], qutrits[k]))
        circuit.append(QutritT02Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX102Gate().on(qutrits[i], qutrits[j]))
        circuit.append(QutritH02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritT02Gate().on(qutrits[i]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT_02Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX102Gate().on(qutrits[i], qutrits[j]))
    if num == "112" or num == "121":
        circuit.append(QutritH12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX112Gate().on(qutrits[j], qutrits[k]))

        circuit.append(QutritT_12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX112Gate().on(qutrits[i], qutrits[k]))

        circuit.append(QutritT12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX112Gate().on(qutrits[j], qutrits[k]))
        circuit.append(QutritT_12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX112Gate().on(qutrits[i], qutrits[k]))
        circuit.append(QutritT12Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX112Gate().on(qutrits[i], qutrits[j]))
        circuit.append(QutritH12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritT12Gate().on(qutrits[i]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT_12Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX112Gate().on(qutrits[i], qutrits[j]))
    if num == "201" or num == "210":
        circuit.append(QutritH01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX201Gate().on(qutrits[j], qutrits[k]))

        circuit.append(QutritT_01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX201Gate().on(qutrits[i], qutrits[k]))

        circuit.append(QutritT01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX201Gate().on(qutrits[j], qutrits[k]))
        circuit.append(QutritT_01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX201Gate().on(qutrits[i], qutrits[k]))
        circuit.append(QutritT01Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX201Gate().on(qutrits[i], qutrits[j]))
        circuit.append(QutritH01Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritT01Gate().on(qutrits[i]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT_01Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX201Gate().on(qutrits[i], qutrits[j]))
    if num == "202" or num == "220":
        circuit.append(QutritH02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX202Gate().on(qutrits[j], qutrits[k]))

        circuit.append(QutritT_02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX202Gate().on(qutrits[i], qutrits[k]))

        circuit.append(QutritT02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX202Gate().on(qutrits[j], qutrits[k]))
        circuit.append(QutritT_02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX202Gate().on(qutrits[i], qutrits[k]))
        circuit.append(QutritT02Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX202Gate().on(qutrits[i], qutrits[j]))
        circuit.append(QutritH02Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritT02Gate().on(qutrits[i]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT_02Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX202Gate().on(qutrits[i], qutrits[j]))
    if num == "212" or num == "221":
        circuit.append(QutritH12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX212Gate().on(qutrits[j], qutrits[k]))

        circuit.append(QutritT_12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX212Gate().on(qutrits[i], qutrits[k]))

        circuit.append(QutritT12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritCX212Gate().on(qutrits[j], qutrits[k]))
        circuit.append(QutritT_12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX212Gate().on(qutrits[i], qutrits[k]))
        circuit.append(QutritT12Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX212Gate().on(qutrits[i], qutrits[j]))
        circuit.append(QutritH12Gate().on(qutrits[k]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritT12Gate().on(qutrits[i]))  # QutritPlusGate().on(qutrits[0])
        circuit.append(QutritT_12Gate().on(qutrits[j]))  # QutritPlusGate().on(qutrits[0])

        circuit.append(QutritCX212Gate().on(qutrits[i], qutrits[j]))
    return circuit


def generate_incrementer(qutrits, i, j, k, lvl):
    #print(cirq)
    #for m_index, moment in enumerate(circ):
    #    print(m_index)

    print("#######################")
 
    circ = cirq.Circuit()

    # circ.append(QutritH01Gate().on(qutrits[k])) 
    # circ.append(QutritCX101Gate().on(qutrits[j], qutrits[k]))
    #print(cirq.unitary(QutritCX101Gate()))
    #print(cirq.unitary(QutritCX101Gate()).shape)
    # circ.append(QutritT_01Gate().on(qutrits[k]))
    # circ.append(QutritCX101Gate().on(qutrits[i], qutrits[k]))
    # circ.append(QutritT01Gate().on(qutrits[k])) # QutritPlusGate().on(qutrits[0])
    #circ.append(QutritCX101Gate().on(qutrits[j], qutrits[k]))
    #circ.append(QutritT_01Gate().on(qutrits[k])) # QutritPlusGate().on(qutrits[0])
    #circ.append(QutritCX101Gate().on(qutrits[i], qutrits[k])
    #circ.append(QutritT01Gate().on(qutrits[j])) # QutritPlusGate().on(qutrits[0])
    #circ.append(QutritT01Gate().on(qutrits[k])) # QutritPlusGate().on(qutrits[0])
    #circ.append(QutritCX101Gate().on(qutrits[i], qutrits[j]))
    #circ.append(QutritH01Gate().on(qutrits[k])) # QutritPlusGate().on(qutrits[0])
    #circ.append(QutritT01Gate().on(qutrits[j])) # QutritPlusGate().on(qutrits[0])
    #circ.append(QutritT_01Gate().on(qutrits[j])) # QutritPlusGate().on(qutrits[0])
    #circ.append(QutritCX101Gate().on(qutrits[i], qutrits[j]))
    #return circ 

    if lvl == 1:
        circ.append(generate_toffoli(qutrits, i, j, k, str(lvl) + "01"))
        # circ.append(generate_toffoli(qutrits, i, j, k, str(lvl) + "02"))
    if lvl == 2:
        circ.append(generate_toffoli(qutrits, i, j, k, str(lvl) + "01"))
        circ.append(generate_toffoli(qutrits, i, j, k, str(lvl) + "02"))
    return circ

def generate_decrementer(qutrits, i, j, k, lvl):
    circ = cirq.Circuit()
    if lvl == 1:
        circ.append(generate_toffoli(qutrits, i, j, k, str(lvl) + "12"))
        circ.append(generate_toffoli(qutrits, i, j, k, str(lvl) + "20"))
    if lvl == 2:
        circ.append(generate_toffoli(qutrits, i, j, k, str(lvl) + "12"))
        circ.append(generate_toffoli(qutrits, i, j, k, str(lvl) + "20"))
    return circ



if __name__ == "__main__":
    # Here we create a qutrit for the gate to act on. 
    q0 = cirq.LineQid(0, dimension=3)
    q1 = cirq.LineQid(1, dimension=3)
    q2 = cirq.LineQid(2, dimension=3)

    q = cirq.LineQid.range(16, dimension=3)


    # We can now enact the gate on this qutrit.
    circuit = cirq.Circuit(
            QutritPlusGate().on(q0)
    )
    circuit2 = cirq.Circuit(
                    QutritC2PlusGate().on(q0, q1)
    )
    circuit3 = cirq.Circuit(
                    QutritCC2PlusGate().on(q0, q1, q2)
    )
    circuit4 = cirq.Circuit(

                    QutritPlusGate().on(q[0]),
                    QutritPlusGate().on(q[1]),
                    QutritPlusGate().on(q[2]),
                    QutritPlusGate().on(q[3]),
                    QutritPlusGate().on(q[4]),
                    QutritPlusGate().on(q[5]),
                    QutritPlusGate().on(q[6]),
                    QutritPlusGate().on(q[7]),
                    QutritPlusGate().on(q[8]),
                    QutritPlusGate().on(q[9]),
                    QutritPlusGate().on(q[10]),
                    QutritPlusGate().on(q[11]),
                    QutritPlusGate().on(q[12]),
                    QutritPlusGate().on(q[13]),
                    QutritPlusGate().on(q[14]),
                    QutritCC1PlusGate().on(q[0], q[2], q[1]),
                    QutritCC1PlusGate().on(q[4], q[6], q[5]),
                    QutritCC1PlusGate().on(q[8], q[10], q[9]),
                    QutritCC1PlusGate().on(q[12], q[14], q[13]),
                    QutritCC2PlusGate().on(q[1], q[5], q[3]),
                    QutritCC2PlusGate().on(q[9], q[13], q[11]),
                    QutritCC2PlusGate().on(q[3], q[11], q[7]),
                    QutritC2PlusGate().on(q[7], q[15]),

                    QutritCC2MinusGate().on(q[3], q[11], q[7]), 
                    
                    QutritCC2MinusGate().on(q[1], q[5], q[3]),
                    QutritCC2MinusGate().on(q[9], q[13], q[11]),

                    QutritCC1MinusGate().on(q[0], q[2], q[1]),
                    QutritCC1MinusGate().on(q[4], q[6], q[5]),
                    QutritCC1MinusGate().on(q[8], q[10], q[9]),
                    QutritCC1MinusGate().on(q[12], q[14], q[13]),
                    
                    cirq.measure(q[7], key='m'),
    )
    #
    # When we print this out we see that the qutrit is labeled by its dimension.
    print(circuit)
    print(circuit2)
    print(circuit3)
    print(circuit4)

    simulator = cirq.Simulator()
    result = simulator.run(circuit4, repetitions=10)

    print(result)
