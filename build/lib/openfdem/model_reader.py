import joblib as jl
import glob
import timeit
from multiprocessing import Process, Queue, Pool,cpu_count


import pyvista as pv
import openfdem as fd

def timed(func):
    def wrapper(*args,**kwargs):
        start = timeit.default_timer()
        func(*args,**kwargs)
        time = timeit.default_timer()-start
        print("%.3g s"%(time))
    return wrapper

@timed
def mp_read(list_of_files):
    return jl.Parallel(n_jobs=-2,prefer="threads")(jl.delayed(pv.read)(file) for file in list_of_files)

@timed
def normal_read(list_of_files):
    data = []
    for file in list_of_files:
        data.append(pv.read(file))
    return data

@timed
def pv_read(list_of_files):
    return pv.read(list_of_files)

def pv_read_queue(list_of_files,q):
    for file in list_of_files:
        q.put(pv.read(file))

@timed
def multiprocess_lib_read(list_of_files):
    # https://stackoverflow.com/questions/20727375/multiprocessing-pool-slower-than-just-using-ordinary-functions
    n_processes = 4
    if __name__ == '__main__':
        q = Queue()

        split_list = [list_of_files[x:x+size] for x in range(0,len(list_of_files),len(list_of_files)//n_processes)]
        for chunk in split_list:
            p = Process(target=pv_read_queue,args=(chunk,q))
            p.daemon = True
            p.start()
        for chunk in split_list:
            p.join()

    return q.get()

@timed
def multiprocess_async(list_of_files):
    p = Pool(processes=cpu_count()-2)
    preresults = [p.apply_async(pv.read,args=(x,)) for x in list_of_files]
    results = [x.get() for x in preresults]
    return results

if __name__ == '__main__':

    file_string = r'C:\Users\emags\Documents\20210622_irazu_example\25Red/*.vtu'

    #benchmark_sizes = [13,25,50,100,200,400,800]
    benchmark_sizes = [50]

    for size in benchmark_sizes:
        list_of_files = glob.glob(file_string)[:size]

        # mp_read(list_of_files)
        

        # normal_read(list_of_files)

        # pv_read(list_of_files)

        # multiprocess_lib_read(list_of_files)

        multiprocess_async(list_of_files)
