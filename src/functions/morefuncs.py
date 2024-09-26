import stim
import re
import os

def extract_qubit_numbers_and_p_value(contents, file_path):
    # Extracting 'p' value from filename
    p_value_match = re.search(r'p=([\d\.]+)', file_path)
    p_value = p_value_match.group(1) if p_value_match else '0.001'  # Default value if not found

    # Extracting qubit numbers
    qubit_numbers = set()
    for line in contents:
        if line.startswith("QUBIT_COORDS"):
            parts = line.split()
            qubit_number = int(parts[-1])
            qubit_numbers.add(qubit_number)
        elif line.startswith("TICK"):
            break  # Stop reading qubit numbers after the first TICK

    return qubit_numbers, p_value

def process_file(file_path):
    with open(file_path, 'r') as file:
        contents = file.readlines()

    qubit_numbers, p_value = extract_qubit_numbers_and_p_value(contents, file_path)

    modified_contents = []
    current_section = []
    in_tick_block = False  # Flag to indicate if we are inside a TICK block

    for line in contents:
        if "TICK" in line:
            if in_tick_block and current_section:  # End of a TICK block
                # Process current section
                present_qubits = set()
                for sect_line in current_section:
                    parts = sect_line.split()
                    for part in parts[1:]:
                        if part.isdigit():
                            present_qubits.add(int(part))

                absent_qubits = qubit_numbers - present_qubits
                if absent_qubits:
                    current_section.append(f"I {' '.join(map(str, absent_qubits))}\n")
                    current_section.append(f"DEPOLARIZE1({p_value}) {' '.join(map(str, absent_qubits))}\n")
                
                modified_contents.extend(current_section)
                current_section = []
                in_tick_block = False  # Exiting TICK block

            modified_contents.append(line)
            in_tick_block = True  # Starting a new TICK block
        elif in_tick_block:
            current_section.append(line)
        else:
            modified_contents.append(line)  # Lines outside TICK blocks are added directly

    if current_section:  # Add any remaining lines after the last TICK
        modified_contents.extend(current_section)

    return modified_contents


def fix_positions_of_idling_errors(input_file_path, output_file_path):  # corrected idling errors (cut two lines directly after REPEAT and pasted them before REPEAT and at the end of the REPEAT block) <-- note this relies on there being a repeat block
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    # Find the line numbers for "REPEAT"
    repeat_line_number = None
    for i, line in enumerate(lines):
        if "REPEAT" in line:
            repeat_line_number = i
            break

    if repeat_line_number is None:
        raise ValueError("The word 'REPEAT' was not found in the file.")

    # Lines to be moved
    lines_to_move = [lines[repeat_line_number + 1], lines[repeat_line_number + 2]]

    # Remove the original lines from their place
    del lines[repeat_line_number + 1:repeat_line_number + 3]

    # Find the line number with the closing curly bracket after deletion
    closing_bracket_line_number = None
    for i, line in enumerate(lines):
        if '}' in line:
            closing_bracket_line_number = i
            break

    if closing_bracket_line_number is None:
        raise ValueError("A closing curly bracket '}' was not found in the file.")

    # Insert the lines above "REPEAT"
    lines = lines[:repeat_line_number] + lines_to_move + lines[repeat_line_number:]

    # Adjust the line number for the closing curly bracket after the first insertion
    closing_bracket_line_number += 2

    # Insert the lines above the closing curly bracket
    lines = lines[:closing_bracket_line_number] + lines_to_move + lines[closing_bracket_line_number:]

    # Save the modified contents to the output file
    with open(output_file_path, 'w') as file:
        file.writelines(lines)


def extract_qubit_numbers_from_top(file_path):
    qubit_numbers = set()
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('QUBIT_COORDS'):
                parts = line.strip().split(' ')
                qubit_number = parts[-1]
                qubit_numbers.add(qubit_number)
            else:
                break  # Stop reading after the qubit numbers section
    return qubit_numbers



def modify_error_probabilities(file_path, new_file_path):
    # Read the content of the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    modified_lines = []
    for i, line in enumerate(lines):
        # DEPOLARIZE2 errors are not modified, hence not checked here
        
        # Modify DEPOLARIZE1 errors
        if 'DEPOLARIZE1' in line:
            if i < len(lines) - 1 and (lines[i+1].strip() == 'REPEAT' or lines[i+1].strip() == '}'):
                modified_line = re.sub(r'DEPOLARIZE1\((.*?)\)', lambda m: f'DEPOLARIZE1({float(m.group(1))*2})', line)
            else:
                modified_line = re.sub(r'DEPOLARIZE1\((.*?)\)', lambda m: f'DEPOLARIZE1({float(m.group(1))/10})', line)
        else:
            modified_line = line
        
        # Modify X_ERROR errors based on position relative to 'M' and 'MR'
        if 'X_ERROR' in line:
            if i > 0 and lines[i-1].strip().startswith('MR'):
                modified_line = re.sub(r'X_ERROR\((.*?)\)', lambda m: f'X_ERROR({float(m.group(1))*5})', line)
            elif i < len(lines) - 1 and (lines[i+1].strip().startswith('M') or lines[i+1].strip().startswith('MR')):
                modified_line = re.sub(r'X_ERROR\((.*?)\)', lambda m: f'X_ERROR({float(m.group(1))*2})', line)
            else:
                modified_line = re.sub(r'X_ERROR\((.*?)\)', lambda m: f'X_ERROR({float(m.group(1))*2})', line)
        
        # Modify Z_ERROR errors based on specific conditions
        if 'Z_ERROR' in line:
            if i > 0 and (lines[i-1].strip().startswith('RX') or lines[i-1].strip().startswith('R')):
                modified_line = re.sub(r'Z_ERROR\((.*?)\)', lambda m: f'Z_ERROR({float(m.group(1))*2})', line)
            elif i < len(lines) - 1 and lines[i+1].strip().startswith('MX'):
                modified_line = re.sub(r'Z_ERROR\((.*?)\)', lambda m: f'Z_ERROR({float(m.group(1))*5})', line)
        
        modified_lines.append(modified_line)
    
    # Write the modified content back to a new file
    with open(new_file_path, 'w') as file:
        file.writelines(modified_lines)


# def modify_x_error_arguments(file_path,new_file_path):
#    # Extract p_value from the file name
#    match = re.search(r'p=([0-9.]+)', file_path)
#    if match:
#        p_value = float(match.group(1))
#    else:
#        print("p_value could not be found in the file name.")
#        return
#    # Read the lines from the file
#    with open(file_path, 'r') as file:
#        lines = file.readlines()
#    # Function to replace the argument of X_ERROR
#    def replace_argument(line, new_value):
#        return re.sub(r'X_ERROR\((.*?)\)', f'X_ERROR({new_value})', line)
#    # Find lines containing "MR" and modify the argument of X_ERROR
#    # in the lines immediately before and after
#    mr_lines = [index for index, line in enumerate(lines) if "MR" in line]
#    for line_number in mr_lines:
#        if line_number > 0:  # Check for preceding line
#            new_value = 5 * p_value
#            lines[line_number - 1] = replace_argument(lines[line_number - 1], f'{new_value:.3f}')
#        if line_number < len(lines) - 1:  # Check for succeeding line
#            new_value = 2 * p_value
#            lines[line_number + 1] = replace_argument(lines[line_number + 1], f'{new_value:.3f}')
   # Write the modified lines back to the file or a new file

def modify_x_error_arguments(file_path,new_file_path):
   # Extract p_value from the file name
   match = re.search(r'p=([0-9.]+)', file_path)
   if match:
       p_value = float(match.group(1))
   else:
       print("p_value could not be found in the file name.")
       return
   # Read the lines from the file
   with open(file_path, 'r') as file:
       lines = file.readlines()
   # Function to replace the argument of X_ERROR
   def replace_argument(line, new_value):
       return re.sub(r'X_ERROR\((.*?)\)', f'X_ERROR({new_value})', line)
   # Find lines containing "MR" or "M" and modify the argument of X_ERROR accordingly
   for index, line in enumerate(lines):
       if "MR" in line:
           if index > 0:  # Preceding line for "MR"
               lines[index - 1] = replace_argument(lines[index - 1], f'{5 * p_value:.3f}')
           if index < len(lines) - 1:  # Succeeding line for "MR"
               lines[index + 1] = replace_argument(lines[index + 1], f'{2 * p_value:.3f}')
       elif "M" in line and not "MR" in line:  # Exclusively for lines with just "M"
           if index > 0:
               lines[index - 1] = replace_argument(lines[index - 1], f'{5 * p_value:.3f}')

   with open(new_file_path, 'w') as file:
       file.writelines(lines)



def modify_depol1_argument(file_path, new_file_path):
   # Extract p_value from the file name
   match = re.search(r'p=([0-9.]+)', file_path)
   if match:
       p_value = float(match.group(1))
   else:
       raise ValueError("p_value could not be found in the file name.")
   # Read the lines from the file
   with open(file_path, 'r') as file:
       lines = file.readlines()
   # Find the line number just before the line that starts with 'REPEAT'
   for index, line in enumerate(lines):
       if line.startswith('REPEAT') and index > 0:
           # Replace the argument of DEPOLARIZE1 in the preceding line
           new_argument = 2 * p_value
           preceding_line = lines[index - 1]
           lines[index - 1] = re.sub(r'DEPOLARIZE1\((.*?)\)', f'DEPOLARIZE1({new_argument})', preceding_line)
           break

   # Write the modified lines back to a new file
   with open(new_file_path, 'w') as new_file:
       new_file.writelines(lines)
   return


def add_idling_errors_and_save_circuit(thecircuit, b, d, p, ro, x, z):
    ## Add idling errors and save circuit:
    path = f"circuits/temp_circuit.stim" # circuit without idling errors
    thecircuit.to_file(path) #saves circuit as a file

    modified_file_contents = process_file(path) # adds some idling errors 
    
    newpath = f"circuits/SD/{b}/d={d},p={p},b={b},r=3d,ro={ro},o={x[0]}{x[1]}{x[2]}{x[3]}{z[0]}{z[1]}{z[2]}{z[3]},idl=y.stim"

    with open(newpath, 'w') as file:
        file.writelines(modified_file_contents)

    fix_positions_of_idling_errors(newpath,newpath) # completes addition of idling errors 

    if os.path.exists(path):
        os.remove(path)


def make_CXSI_circuit(b, d, p, ro, x, z):
    
    newpath = f"circuits/SD/{b}/d={d},p={p},b={b},r=3d,ro={ro},o={x[0]}{x[1]}{x[2]}{x[3]}{z[0]}{z[1]}{z[2]}{z[3]},idl=y.stim" # needs to be the same as 'newpath' in add_idling_errors_and_save_circuit
    
    
    CXSI_file_path = f"circuits/CXSI/{b}/d={d},p={p},noise=CXSI,b={b},r=3d,ro={ro},o={x[0]}{x[1]}{x[2]}{x[3]}{z[0]}{z[1]}{z[2]}{z[3]},idl=y.stim"
    modify_error_probabilities(newpath,CXSI_file_path)
    modify_x_error_arguments(CXSI_file_path,CXSI_file_path)
    modify_depol1_argument(CXSI_file_path,CXSI_file_path)
