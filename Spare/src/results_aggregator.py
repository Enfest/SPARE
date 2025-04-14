import numpy as np
import json
import os
import sys

# print(len(sys.argv))
home = sys.argv[1]

#results_order = ["adder", "interger", "wa", "mult"]
#results_order = ["length_simplified", "remove", "push_back", "sum_orig", "sum"]
results_order = ["is_", "pop", "remove", "push_back", "sum_orig", "sum", "simpl_"]

def print_init(f):
    f.write("\\begin{table*}[t]\n")
    f.write("\\centering\n")
    f.write("\\footnotesize\n")
    f.write("\\begin{tabular}{|l|lll|lll|lll|}\n")
    f.write("\\hline\n")
    f.write("& \\multicolumn{3}{c|}{\\manualbl} & " +
            " \\multicolumn{3}{c|}{\toolbl} & " +
            " \\multicolumn{3}{c|}{\\textbf{relative change}} \\\\\n")
    f.write("\\hline\n")
    f.write("benchmark & 2-qutrit gates & 1-qutrit gates & depth  & " +
            " 2-qutrit gates & 1-qutrit gates &  depth & 2-qutrit $\Delta$ &" +
            " 1-qutrit $\Delta$ & depth $\\Delta$ \\\\\n")
    f.write("\\hline\n")

write_path = home + "/log_files/manual_and_depth_spire.txt"
results_folder = "TOOL_RESULTS_ANCILLA2"
manual_folder = "MANUAL_RESULTS_IMPROVED" # "MANUAL_RESULTS_IMPROVED"
folder_path = home + "/OUTPUT_JSON/" + results_folder # TOOL_RESULTS/"   # Replace this with the path to your folder
files = os.listdir(folder_path)
print(files)

with open(write_path, 'w') as f:
    print_init(f)
    files_list = []
    for tests in files:
        files_list.append(tests)
    files_list.sort()
    print(files_list)
    print(results_order)
    for file_types in results_order:
        filtered_list = list(filter(lambda s: file_types in s, files_list))
        filtered_list.sort()
        for tests in filtered_list:
            print(tests)
            file_path = os.path.join(folder_path, tests)
            if os.path.isfile(file_path) and ".json" in file_path:
                try:
                    with open(file_path, 'r') as file:
                        print("curr file ", file_path)
                        json_data = json.load(file)
                        if "random_compilation_itr_final" in json_data:
                            depth = json_data.get("random_compilation_itr_final").get("depth")
                            q2gates = json_data.get("random_compilation_itr_final").get("2_q_gates")
                            q1gates = json_data.get("random_compilation_itr_final").get("1_q_gates")
                            print("\twe found:", depth, q2gates)
                            print("\t", results_folder in file_path)
                            print("\tcheck ", file_path.replace(results_folder, manual_folder).replace("output_", "output_"))
                            print("___________________")
                            if os.path.isfile(file_path.replace(results_folder, manual_folder).replace("output_", "output_")):
                                with open(file_path.replace(results_folder, manual_folder).replace("output_", "output_")) as file2:
                                    json_data = json.load(file2)
                                    depth_m = json_data.get("manual_rewrites").get("depth")
                                    q2gates_m = json_data.get("manual_rewrites").get("2_q_gates")
                                    q1gates_m = json_data.get("manual_rewrites").get("1_q_gates")
                                    f.write("\\bmark{" + tests.replace("output_", "").replace(".json", "").replace("_", "-") +
                                            "} & " + str(q2gates_m) + " & " + str(q1gates_m) + " & " +
                                            str(depth_m) + "  &  " + str(q2gates) + " & " + str(q1gates) +
                                            " & " + str(depth) + " & " + str(int(q2gates) - int(q2gates_m)) +
                                            " & " + str(int(q1gates) - int(q1gates_m)) +
                                            " & " + str(int(depth) - int(depth_m)) +"\\\\\n")
                            else:
                                f.write("\\bmark{"+ tests.replace("output_", "").replace(".json", "").replace("_", "-") +
                                        "} & - & - & -  &  " + str(q2gates) + " & " + str(q1gates) + " & " + str(depth) +
                                        " & 0 & 0 & 0\\\\\n")
                except:
                    pass
    f.write("\end{tabular}\n")
    f.write("\end{table*}")
f.close()

# Now writing the manual qubit vs qutrits
write_path = home + "/log_files/manual_toffoli_and_depth4.txt"
folder_path = home + "/OUTPUT_JSON/TOOL_RESULTS/"   # Replace this with the path to your folder
files = os.listdir(folder_path)

with open(write_path, 'w') as f:
    print_init(f)
    files_list = []
    for tests in files:
        files_list.append(tests)
    
    for file_types in results_order:
        filtered_list = list(filter(lambda s: file_types in s, files_list))
        filtered_list.sort()
        for tests in filtered_list:
            file_path = os.path.join(folder_path, tests)
            if os.path.isfile(file_path) and ".json" in file_path:
                try:
                    with open(file_path, 'r') as file:
                        json_data = json.load(file)
                        if "random_compilation_itr_final" in json_data:
                            depth = json_data.get("random_compilation_itr_final").get("depth")
                            q2gates = json_data.get("random_compilation_itr_final").get("2_q_gates")
                            q1gates = json_data.get("random_compilation_itr_final").get("1_q_gates")
                            if os.path.isfile(file_path.replace("TOOL", "MANUAL_TOFFOLI").replace("output_", "output_manual")):
                                with open(file_path.replace("TOOL", "MANUAL_TOFFOLI").replace("output_", "output_manual")) as file2:
                                    json_data = json.load(file2)
                                    depth_m = json_data.get("manual_rewrites").get("depth")
                                    q2gates_m = json_data.get("manual_rewrites").get("2_q_gates")
                                    q1gates_m = json_data.get("manual_rewrites").get("1_q_gates")
                                    f.write("\\bmark{" + tests.replace("output_", "").replace(".json", "").replace("_", "-") +
                                            "} & " + str(q2gates_m) + " & " + str(q1gates_m) + " & " +
                                            str(depth_m) + "  &  " + str(q2gates) + " & " + str(q1gates) +
                                            " & " + str(depth) + " & " + str(int(q2gates) - int(q2gates_m)) +
                                            " & " + str(int(q1gates) - int(q1gates_m)) +
                                            " & " + str(int(depth) - int(depth_m)) +"\\\\\n")
                            else:
                                pass
                except:
                    pass
    f.write("\end{tabular}\n")
    f.write("\end{table*}")

# Simulation results
if True:
    write_path = home + "/log_files/simulation_results4.txt"
    folder_path = home + "/OUTPUT_JSON/TOOL_RESULTS/"   # Replace this with the path to your folder
    files = os.listdir(folder_path)

    with open(write_path, 'w') as f:
        print_init(f)
        files_list = []
        for tests in files:
            files_list.append(tests)
        
        for file_types in results_order:
            filtered_list = list(filter(lambda s: file_types in s, files_list))
            filtered_list.sort()
            for tests in filtered_list:
                file_path = os.path.join(folder_path, tests)
                if os.path.isfile(file_path) and ".json" in file_path:
                    # print("file path is ", file_path)
                    try:
                        with open(file_path, 'r') as file:
                            json_data = json.load(file)
                            if "random_compilation_itr_final" in json_data:
                                f.write("\\bmark{" + tests.replace("output_", "").replace(".json", "").replace("_", "-") + "} & ") 

                                fidelity_sc = json_data.get("random_compilation_itr_final").get("fidelity").get("current_sc")
                                fidelity_ti = json_data.get("random_compilation_itr_final").get("fidelity").get("dressed_ti")
                                fidelity_stab_2 = json_data.get("random_compilation_itr_final").get("fidelity").get("removed_2_scaling")
                                fidelity_future = json_data.get("random_compilation_itr_final").get("fidelity").get("all_scaled")
                                if os.path.isfile(file_path.replace("TOOL", "MANUAL_TOFFOLI").replace("output_", "output_manual")):
                                    with open(file_path.replace("TOOL", "MANUAL_TOFFOLI").replace("output_", "output_manual")) as file2:
                                        qmanual_fidelity_sc = json_data.get("random_compilation_itr_final").get("")
                                        qmanual_fidelity_ti = json_data.get("random_compilation_itr_final").get("")
                                        qmanual_fidelity_future = json_data.get("random_compilation_itr_final").get("")
                                        f.write(str(qmanual_fidelity_sc) + "& " + str(qmanual_fidelity_ti) + " & " + str(qmanual_fidelity_future) + "& ") 
                                else:
                                    f.write(" - & - & - & ")
                                if os.path.isfile(file_path.replace("TOOL", "MANUAL").replace("output_", "output_manual_")):
                                    with open(file_path.replace("TOOL", "MANUAL").replace("output_", "output_manual_")) as file2:
                                        tmanual_fidelity_sc = json_data.get("random_compilation_itr_final").get("fidelity").get("current_sc")
                                        tmanual_fidelity_ti = json_data.get("random_compilation_itr_final").get("fidelity").get("dressed_ti")
                                        tmanual_fidelity_stab_2 = json_data.get("random_compilation_itr_final").get("fidelity").get("removed_2_scaling")
                                        tmanual_fidelity_future = json_data.get("random_compilation_itr_final").get("fidelity").get("all_scaled")
                                        f.write(str(tmanual_fidelity_sc) + "& " + str(tmanual_fidelity_ti) + \
                                                " & " + str(tmanual_fidelity_stab_2) + " & " + str(tmanual_fidelity_future) + " & ")
                                else:
                                    f.write(" - & - & - & - & ")
                                f.write(str(fidelity_sc) + "& " + str(fidelity_ti) + " & " + str(fidelity_stab_2) + " & " + str(fidelity_future) + "\\\\\n")
                                print(str(tmanual_fidelity_sc) + "& " + str(tmanual_fidelity_ti) + " & " + str(tmanual_fidelity_stab_2) + " & " + str(tmanual_fidelity_future) + " & " + str(fidelity_sc) + "& " + str(fidelity_ti) + " & " + str(fidelity_stab_2) + " & " + str(fidelity_future))
                    except:
                        pass
        f.write("\end{tabular}\n")
        f.write("\end{table*}")
