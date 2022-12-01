import openfdem as fdem
import numpy as np

# Load data
data = fdem.Model("/hdd/home/aly/Desktop/Dropbox/Python_Codes/OpenFDEM-Post-Processing/example_outputs/Irazu_UCS")

def load_data_on_DFN(cell_idx):
    # Extract cell information based on Cell ID => What you want to extract
    data_along_DFN = data.extract_cell_info(cell_idx, ['mass', 'principal stresses'])

    # For each value in the 1d stress
    mag_list = []
    strain_eng = []
    for i in data_along_DFN['principal stresses']:
        # Convert the stress array from 1d to 3x3
        stress_array = np.array(i).reshape(3,3)
        # Calculates the Magnitude
        magnitude = np.linalg.norm(stress_array)
        # Add magnitude to list
        mag_list.append(magnitude)
        stress_array = 0

    for idx, i in enumerate(data_along_DFN['mass']):
        # Calculate strain energy
        strain_eng_cal = ( (i / 2800) * ((mag_list[idx]**2) / 50000000000) ) * 0.5
        strain_eng.append(strain_eng_cal)

    # Add the list to the dataframe
    data_along_DFN['pstress_mag'] = mag_list
    data_along_DFN['strain_energy_release'] = strain_eng

    return  data_along_DFN

# # Get information for a particular TimeStep
# print(data_along_DFN.loc[10])

# Cell ID along a line.
start = [10, 10, 0]
end_point = [15, 15, 0]

# Get the Cell ID based on XYZ coordinate
cell_on_line = []
for idx, i in enumerate(np.arange(start[0], end_point[0], 0.1)):
    newx = start[0] + (0.1 * (idx +1))
    newy = start[0] + (0.1 * (idx +1))
    cellid = data.find_cell([newx,newy,0])
    cell_on_line.append(cellid)

# print(cell_on_line)
# Get ID of unique cells on the line
cell_on_line = list(set(cell_on_line))
# print(cell_on_line)

# Extract data of cell along the line
DFN_line1 = {}
for cell in cell_on_line:
    print("Extracting data along Cell ID - %s" % cell)
    data_along_DFN = load_data_on_DFN(cell)
    DFN_line1[cell] = data_along_DFN

# To print/get certain information 
# DICT [CELL ID] [DICT KEY/HEADING] [TIMESTEP]
print(DFN_line1[970]['strain_energy_release'][10])