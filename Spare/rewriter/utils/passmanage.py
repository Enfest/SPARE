from src.cirq_depth_counter import *
from rewriter.src.compile_basic import *
import rewriter.src.compile_basic as basic_passes

class passPair:
    def __init__(self, finderpass, rewritePass, target="all", if_spec=False, times=1):
        self._finder_pass = finderpass
        self._rewrite_pass = rewritePass
        self._times = times
        self._target = target
        self._if_spec = if_spec

    def applyPass(self, circuit_obj):
        if not self._if_spec:
            circuit_obj = self.applyDirectPass(circuit_obj)
        else:
            circuit_obj = self.applyPassSpeculatively(circuit_obj)
        return circuit_obj 
    
    def applyDirectPass(self, circuit_obj):
        counter = 0
        if self._times == float("inf"):            
            while True:
                circuit_obj_copy = circuit_obj.deepcopy()
                _, result, _ = self._applyDirectPass(circuit_obj)
                if not result:
                    return circuit_obj_copy
                counter += 1
        elif isinstance(self._times, int):
            for t in range(self._times):
                circuit_obj_copy = circuit_obj.deepcopy()
                _, result, _ = self._applyDirectPass(circuit_obj)
                if not result:
                    circuit_obj = circuit_obj_copy
                counter += 1
        return circuit_obj
    
    def applyPassSpeculatively(self, circuit_obj):
        speculative_history = []
        counter = 0
        if isinstance(self._times, int):
          for t in range(self._times):
            circuit_obj, speculative_history, _ = self._applyPassSpeculatively(circuit_obj, speculative_history)
            counter += 1
        else:
            raise NotImplementedError
        return circuit_obj
      
    def _applyDirectPass(self, circuit_obj, speculative_history=[]):
        try:
            target_opt_func = getattr(circuit_obj, self._rewrite_pass)
            if self._finder_pass is not None:
                target_finding_func = getattr(circuit_obj, self._finder_pass)
                output_targets = target_finding_func()
                if len(output_targets) == 0:
                    return circuit_obj, False, speculative_history                 
                if self._target == "all":
                    for selected_target in output_targets:
                        res = target_opt_func(selected_target)    
                else:
                    selected_target = self._target(output_targets)
                    res = target_opt_func(selected_target)
                    if output_targets.append(selected_target) in speculative_history:
                        return circuit_obj, "Already_in_hist", speculative_history                    
                    speculative_history.append(output_targets.append(selected_target))
            else:
                '''Since no parent nothing added to speculative_history'''
                res = target_opt_func()
        except:
            target_opt_func = getattr(basic_passes, self._rewrite_pass)
            circuit_obj = target_opt_func(circuit_obj)
            return circuit_obj, True, speculative_history
        return circuit_obj, res, speculative_history
    
    def _applyPassSpeculatively(self, circuit_obj, speculative_history):
        circuit_copy = circuit_obj.deepcopy()
        _, result, speculative_history = self._applyDirectPass(circuit_copy, speculative_history)
        if result == "Already_in_hist":
            # Found a target already explored
            return circuit_obj, speculative_history, True
        elif not result:
            return circuit_obj, speculative_history, False
        circuit_copy = gateslices_simplify(circuit_copy)
        orig_estimate_value = depth_estimater(circuit_obj.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
        new_estimate_value = depth_estimater(circuit_copy.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
        old2gates = toffoli_estimator(circuit_obj.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
        new2gates = toffoli_estimator(circuit_copy.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
        if new_estimate_value < orig_estimate_value or new2gates < old2gates:
            return circuit_copy, [], True
        else:
            return circuit_obj, speculative_history, True


class rewritePassManage:
    def __init__(self, passesLists, if_spec=False, times=1):
        # self._target = target
        assert isinstance(passesLists, list)
        self._pass_list = passesLists
        self._times = times
        self._style = if_spec
        
    def applyPass(self, circuit_obj):
        if not self._style:
            circuit_obj = self.applyDirectPass(circuit_obj)
        else:
            circuit_obj = self.applyPassSpeculatively(circuit_obj)
        return circuit_obj 
    
    def applyDirectPass(self, circuit_obj):
        if self._times == float("inf"):
            while True:
                orig_count = toffoli_estimator(circuit_obj.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
                for _curr_pass in self._pass_list:
                    '''Enumerate all passes here in the direct pass list'''
                    circuit_obj = _curr_pass.applyPass(circuit_obj)
                new_count = toffoli_estimator(circuit_obj.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
                if new_count == orig_count:
                    break
        elif isinstance(self._times, int):
            for t in range(self._times):
                for _curr_pass in self._pass_list:
                    '''Enumerate all passes here in the direct pass list'''
                    # print("\tInstance is ", t)
                    circuit_obj = _curr_pass.applyPass(circuit_obj)
        return circuit_obj
    
    def applyPassSpeculatively(self, circuit_obj):
        if isinstance(self._times, int):
          for t in range(self._times):
            circuit_obj, _ = self._applyPassSpeculatively(circuit_obj)
        else:
            raise NotImplementedError
        return circuit_obj
          
    def _applyPassSpeculatively(self, circuit_obj):
        circuit_copy = circuit_obj.deepcopy()
        c, _ = circuit_obj.build_circuit(breakdown="vvv")
        for _curr_pass in self._pass_list:
            circuit_copy = _curr_pass.applyPass(circuit_copy)
        circuit_copy = gateslices_simplify(circuit_copy)
        orig_estimate_value = depth_estimater(circuit_obj.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
        new_estimate_value = depth_estimater(circuit_copy.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
        old2gates = toffoli_estimator(circuit_obj.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
        new2gates = toffoli_estimator(circuit_copy.build_circuit(breakdown="vvv")[0], circuit_obj.build_circuit(breakdown="vvv")[1], mode="None")
        if new_estimate_value < orig_estimate_value or new2gates < old2gates:
            return circuit_copy, True
        else:
            return circuit_obj, True        


class rewritePassesManager:
    def __init__(self, all_passes_list):
        self._all_pass_List = all_passes_list
        
    def addPass(self, new_pass):
        assert isinstance(new_pass, rewritePassManage)
        self._all_pass_List.append(new_pass)
    
    def applyAllPasses(self, circuit_obj):
        for p in self._all_pass_List:
            circuit_obj = p.applyPass(circuit_obj)
        return circuit_obj