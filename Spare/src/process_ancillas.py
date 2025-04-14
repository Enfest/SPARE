import os
import ast

def get_last_three_lines(file_path):
    """Read the last three lines from the file."""
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()[-3:]
            return lines
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return []

def parse_and_get_list_lengths(line, file):
    """Parse the first line for two lists and return their lengths."""
    try:
        # Assume the line contains two lists in a recognizable format, e.g., '[...] [...]'
        lists = line.split("] [")
        lists[0] = ast.literal_eval(lists[0] + "]")
        lists[1] = ast.literal_eval("[" + lists[1])
        if isinstance(lists, (list, tuple)) and len(lists) == 2:
            return len([x for x in lists[0] if x != 0]), len([x for x in lists[1] if x != 0])
    except Exception as e:
        print(file.split("/")[-1], end=" ")
        print(f"Error parsing line: {e}")
    return None, None

def process_directory(directory_path):
    """Process all files in all subdirectories."""
    for root, _, files in os.walk(directory_path):
        if ("new_" in root or "integer" in root or "wa" in root or "adder" in root) and "12" in root and "orig" not in root:  
            for file_name in files:
                if ("check2" in file_name or "manual_" in file_name or "log2" in file_name): 
                    # print(root, file_name)
                    file_path = os.path.join(root, file_name)
                    last_three_lines = get_last_three_lines(file_path)
                    if last_three_lines:
                        if "manual" in file_name:
                            first_line = last_three_lines[1]
                        else:
                            first_line = last_three_lines[0]    
                        length1, length2 = parse_and_get_list_lengths(first_line, root)
                        if length1 is not None and length2 is not None:
                            print(root.split("/")[-1], file_name, end=" ")
                            print(f"Lengths of lists: {length1}, {length2}", end="; ")
            print("")
# Example usage
directory_path = "/nobackupkiwi/rsharma3/QC_work/sara_work/qutrits/cirq_direct_clone/dare_check/log_files/"
process_directory(directory_path)
