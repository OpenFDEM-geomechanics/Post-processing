import glob
import os
import concurrent.futures

import meshio as mio
from tqdm import tqdm

def convert2bin(file):
    mesh = mio.read(file)
    newdir = os.path.join(os.path.dirname(file),"result")
    if not os.path.exists(newdir):
        os.mkdir(newdir)
    newpath = os.path.join(newdir,os.path.basename(file))
    mesh.write(newpath,binary=True)
    return True

if __name__ == '__main__':
    files = glob.glob(r"X:\20210622_irazu_example\UCS_example\*basic.vtu")
#     print(len(files))
#     with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
#         result = list(tqdm(executor.map(convert2bin,files,chunksize=len(files)//16),total=len(files)))

    result = [convert2bin(file) for file in files]