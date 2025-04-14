import cirq
import math, numpy as np


def apply_pseudo_C1C1PlusOne2(qutrits, i, j, k):
    circ = cirq.Circuit()
    circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    circ.append(C1_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    circ.append(C1_F12().on(qutrits[i], qutrits[k]))
    circ.append(Ry_12(math.pi/4).on(qutrits[k]))
    circ.append(C1_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12(math.pi/4).on(qutrits[k]))
    #circ.append(C1_F01().on(qutrits[j], qutrits[k]))
    #circ.append(Ry_01(math.pi/4).on(qutrits[k]))
    #circ.append(C1_F01().on(qutrits[i], qutrits[k]))
    #circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    #circ.append(C1_F01().on(qutrits[j], qutrits[k]))
    #circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    return circ

def apply_pseudo_C1C1PlusOne1(qutrits, i, j, k):
    circ = cirq.Circuit()
    # circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    # circ.append(C1_F12().on(qutrits[j], qutrits[k]))
    # circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    # circ.append(C1_F12().on(qutrits[i], qutrits[k]))
    # circ.append(Ry_12(math.pi/4).on(qutrits[k]))
    # circ.append(C1_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01(math.pi/4).on(qutrits[k]))
    circ.append(C1_F01().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01(math.pi/4).on(qutrits[k]))
    circ.append(C1_F01().on(qutrits[i], qutrits[k]))
    circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    circ.append(C1_F01().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    return circ

def apply_C1C1PlusOne(qutrits, i, j, k):
    circ = cirq.Circuit()
    circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    circ.append(C1_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    circ.append(C1_F12().on(qutrits[i], qutrits[k]))
    circ.append(Ry_12(math.pi/4).on(qutrits[k]))
    circ.append(C1_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12Pi4Ry01Pi4().on(qutrits[k]))
    circ.append(C1_F01().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01(math.pi/4).on(qutrits[k]))
    circ.append(C1_F01().on(qutrits[i], qutrits[k]))
    circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    circ.append(C1_F01().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    return circ

def apply_C2C2PlusOne(qutrits, i, j, k):
    circ = cirq.Circuit()
    circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    circ.append(C2_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    circ.append(C2_F12().on(qutrits[i], qutrits[k]))
    circ.append(Ry_12(math.pi/4).on(qutrits[k]))
    circ.append(C2_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12Pi4Ry01Pi4().on(qutrits[k]))
    circ.append(C2_F01().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01(math.pi/4).on(qutrits[k]))
    circ.append(C2_F01().on(qutrits[i], qutrits[k]))
    circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    circ.append(C2_F01().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    return circ

def apply_C1C1MinusOne(qutrits, i, j, k):
    circ = cirq.Circuit()
    circ.append(Ry_01(math.pi/4).on(qutrits[k]))
    circ.append(C1_F01().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01(math.pi/4).on(qutrits[k]))
    circ.append(C1_F01().on(qutrits[i], qutrits[k]))
    circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    circ.append(C1_F01().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01MinusPi4Ry12Pi4().on(qutrits[k]))
    circ.append(C1_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    circ.append(C1_F12().on(qutrits[i], qutrits[k]))
    circ.append(Ry_12(math.pi/4).on(qutrits[k]))
    circ.append(C1_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12(math.pi/4).on(qutrits[k]))
    return circ


def apply_C2C2MinusOne(qutrits, i, j, k):
    circ = cirq.Circuit()
    circ.append(Ry_01(math.pi/4).on(qutrits[k]))
    circ.append(C2_F01().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01(math.pi/4).on(qutrits[k]))
    circ.append(C2_F01().on(qutrits[i], qutrits[k]))
    circ.append(Ry_01(-math.pi/4).on(qutrits[k]))
    circ.append(C2_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_01MinusPi4Ry12Pi4().on(qutrits[k]))
    circ.append(C2_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12(-math.pi/4).on(qutrits[k]))
    circ.append(C2_F12().on(qutrits[i], qutrits[k]))
    circ.append(Ry_12(math.pi/4).on(qutrits[k]))
    circ.append(C2_F12().on(qutrits[j], qutrits[k]))
    circ.append(Ry_12(math.pi/4).on(qutrits[k]))
    return circ


class Ry_01(cirq.Gate):
    def _qid_shape_(self):
            return (3,)

    def __init__(self, theta):
        self.theta = theta

    def _unitary_(self):
        theta = self.theta
        return np.array([[math.cos(theta/2.0), -math.sin(theta/2.0), 0],
                         [math.sin(theta/2.0), math.cos(theta/2.0), 0],
                         [0, 0, 1]])
    
    def _circuit_diagram_info_(self, args):
        return '[(Ry01)]'



class Ry_12(cirq.Gate):
    def _qid_shape_(self):
        return (3,)

    def __init__(self, theta):
        self.theta = theta

    def _unitary_(self):
        theta = self.theta
        return np.array([[1, 0, 0],
                         [0, math.cos(theta/2.0), -math.sin(theta/2.0)],
                         [0, math.sin(theta/2.0), math.cos(theta/2.0)]])

    def _circuit_diagram_info_(self, args):
        return '[(Ry12)]'

class Ry_01MinusPi4Ry12Pi4(cirq.Gate):
    def _qid_shape_(self):
        return (3,)

    def _unitary_(self):
        return np.dot(Ry_12(-math.pi/4)._unitary_(), Ry_01(-math.pi/4)._unitary_())

    def _circuit_diagram_info_(self, args):
        return '[Ry01+Ry12+]'


class Ry_12Pi4Ry01Pi4(cirq.Gate):
    def _qid_shape_(self):
        return (3,)

    def _unitary_(self):
        return np.dot(Ry_01(math.pi/4)._unitary_(), Ry_12(math.pi/4)._unitary_())

    def _circuit_diagram_info_(self, args):
        return '[Ry12Ry01]'

class C1_F01(cirq.Gate):
    def _qid_shape_(self):
        return (3,3)

    def _unitary_(self):
        return np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1]])
    
    def _circuit_diagram_info_(self, args):
        return '[C1]', '[F01]'


class C2_F01(cirq.Gate):
    def _qid_shape_(self):
        return (3,3)

    def _unitary_(self):
        return np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1]])

    def _circuit_diagram_info_(self, args):
        return '[C2]', '[F01]'


class C1_F12(cirq.Gate):
    def _qid_shape_(self):
        return (3,3)

    def _unitary_(self):
        return np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1]])

    def _circuit_diagram_info_(self, args):
        return '[C1]', '[F12]'


class C2_F12(cirq.Gate):
    def _qid_shape_(self):
        return (3,3)

    def _unitary_(self):
        return np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0]])

    def _circuit_diagram_info_(self, args):
        return '[C2]', '[F12]'


class C1_PlusOne(cirq.Gate):
    def _qid_shape_(self):
        return (3,3)

    def _unitary_(self):
        return np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1]])
    
    def _circuit_diagram_info_(self, args):
        return '[C1]', '[+1]'


class C1_MinusOne(cirq.Gate):
    def _qid_shape_(self):
        return (3,3)

    def _unitary_(self):
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
        return '[C1]', '[-1]'
