import cirq
import re

def depth_counter(circuit, qubits):
    depth_all = [0] * len(qubits)
    depth_2 = [0] * len(qubits)
    qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
    for a in qd_map:
        if isinstance(a, cirq.LineQid):
            qd_map[a] = a.x
        elif isinstance(a, cirq.NamedQubit):
            qd_map[a] = int(re.findall(r'\d+', a.name)[-1])
        elif isinstance(a, cirq.GridQubit):
            qd_map[a] = int(a.row)    
        else:
            raise Exception("Not supported qubit type")

    for k in qd_map.keys():
        depth_all[qd_map[k]] = 0
        depth_2[qd_map[k]] = 0
    ttl_1_q = 0
    ttl_2_q = 0
    for m_index, moment in enumerate(circuit):
        for op in moment: 
            qudit_indices = [qd_map[qd] for qd in op.qubits]
            if len(qudit_indices) >=2:
                ttl_2_q += 1
            else:
                ttl_1_q += 1
            # print(qudit_indices)
            # print(depth_2)
            if len(qudit_indices) > 1:
                depth_2[qudit_indices[0]] = max(depth_2[qudit_indices[0]], depth_2[qudit_indices[1]]) + 1
                # max(depth_2[qudit_indices[0]], depth_2[qudit_indices[1]]) + 1
                depth_2[qudit_indices[1]] = depth_2[qudit_indices[0]]
                # max(depth_2[qudit_indices[0]], depth_2[qudit_indices[1]]) + 1  # 1
                depth_all[qudit_indices[0]] = max(depth_all[qudit_indices[0]], depth_all[qudit_indices[1]]) + 1
                depth_all[qudit_indices[1]] = depth_all[qudit_indices[0]]
                # (depth_all[qudit_indices[0]], depth_all[qudit_indices[1]]) + 1
            else:
                depth_all[qudit_indices[0]] += 1
        # print(moment)
        # print("depth is ", depth_all, depth_2)
    print("TOTAL 1 qubit gates ", ttl_1_q , "total 2 qubit gates", ttl_2_q)
    print(depth_all, depth_2)
    # single_gate_depth, double_gate_depth, final_single_gate_cnt, final_single_gate_cnt
    return max(depth_all), depth_2, ttl_1_q, ttl_2_q

def get_depth_of_gate(num_of_indices, idx, mode):
    depth = -1
    
    class PythonSwitch:
        def get(self, num_qubits):
            default = []
            if num_qubits < 8:
                return getattr(self, "case_" + str(num_qubits), lambda: default)()
            else:
                return self.case_n(num_qubits)
        
        def case_1(self):
            return [1]

        def case_2(self):
            if mode == "qutrits":
                return [1, 1]
            else:
                return [11, 11]
        
        def case_3(self):
            if mode == "qutrits":
                return [3, 6, 9]
            else:
                return [11, 11, 10]
        
        def case_4(self):
            return [6, 9, 18, 36]

        def case_5(self):
            return [9, 18, 36, 40, 72]
        
        def case_6(self):
            return [3, 6, 9, 18, 36, 72]
        
        def case_7(self):
            return [3, 6, 9, 18, 36, 72, 144]
        
        def case_8(self):
            return [2, 4, 8, 16, 32, 64, 64, 190]
        
        def case_n(self, num):
            if num == 3:
                return [11, 11, 10]
            result = [0]*num
            arr1 = self.case_n(num-1)
            for i in range(1, len(arr1)+1):
                result[i] += arr1[i-1]
            result[0] = 1 + max(arr1[0], arr1[-1])
            result[-1] = 1 + max(arr1[0], arr1[-1])
            for i in range(1, len(arr1)+1):
                result[i] += arr1[i-1]
            result[-1] = 1 + max(arr1[0], arr1[-1])
            return result
    
    my_switch = PythonSwitch()
    final_estimate = my_switch.get(num_of_indices)
    # print("final estimate ", idx, num_of_indices, final_estimate)
    return final_estimate[idx]
    # final_estimate = case.get(num_of_indices, [])
    # if len(final_estimate) == 0:
    #     raise Exception("The estimate of the gate depth is not known")
    # return final_estimate[idx]

def depth_estimater(circuit, qubits, mode="qutrits"):
    if "qubit" not in mode:
        mode = "qutrits"
    depth_all = [0] * len(qubits)
    qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
    for a in qd_map.keys():
        qd_map[a] = a.x    
    for m_index, moment in enumerate(circuit):
        for op in moment:
            # print(op)
            qudit_indices = [qd_map[qd] for qd in op.qubits]
            # print("qudit indexes are ", qudit_indices)
            if len(qudit_indices) > 1:
                max_depth_now = -1
                min_idx = 10000
                for idx in qudit_indices:
                    max_depth_now = max(depth_all[idx], max_depth_now)
                    min_idx = min(min_idx, idx)
                
                for i in range(len(qudit_indices)):
                    idx = qudit_indices[i]
                    depth_all[idx] = max_depth_now + get_depth_of_gate(len(qudit_indices), i, mode)
            else:
                depth_all[qudit_indices[0]] += 1
    return max(depth_all)

def toffoli_estimator(circuit, qubits, mode="qutrits"):
    if "qubit" not in mode:
        mode = "qutrits"
    depth_all = [0] * len(qubits)
    qd_map = {qd: i for i, qd in enumerate(circuit.all_qubits())}
    for a in qd_map.keys():
        qd_map[a] = a.x
    gate_counter = 0
    for m_index, moment in enumerate(circuit):
        for op in moment:
            qudit_indices = len([qd_map[qd] for qd in op.qubits])
            if qudit_indices == 2:
                gate_counter += 1
            elif qudit_indices > 2:
                gate_counter += (2*qudit_indices - 5)*6
    return gate_counter