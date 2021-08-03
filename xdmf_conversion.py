import glob
import os
import concurrent.futures

import meshio as mio
from tqdm import tqdm

def write_xdmf_data(file,writer):
    mesh = mio.read(file)
    timestep = os.path.splitext(os.path.basename(file))[0].split('_')[-1]
    writer.write_data(timestep,cell_data=mesh.cell_data,point_data=mesh.point_data)

    
    return True

if __name__ == '__main__':
    files = glob.glob(r"X:\20210622_irazu_example\UCS_example\*basic*.vtu")
#     print(len(files))
#     with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
#         result = list(tqdm(executor.map(convert2bin,files,chunksize=len(files)//16),total=len(files)))
    xdmf_file = r"X:\20210622_irazu_example\UCS_example\aggregated.xdmf"
    with mio.xdmf.TimeSeriesWriter(xdmf_file) as writer:
        mesh0 = mio.read(files[0])

        writer.write_points_cells(mesh0.points,mesh0.cells)
        result = [write_xdmf_data(file,writer) for file in files]