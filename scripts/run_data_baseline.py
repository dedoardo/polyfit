#
# Runs 'polyfit.exe' on all images in 'data' writing the outputs to 'data-baseline'
# No need to create any of the directories.
# The script is assumed to be run from the root directory (i.e. python scripts/run_data_baseline.py)
#
import os
import multiprocessing

POLYFIT_EXE = 'polyfit.exe'
DATA_DIR = 'data'
DATA_BASELINE_DIR = 'data-baseline'

def run_polyfit_multiple_inputs(inputs):
    while True:
        i = inputs.get()
        if i == -1:
            return
        cmd = f'{POLYFIT_EXE} {i[0]} {i[1]}'
        os.system(cmd)

def run_data_baseline(polyfit_exe, data_dir, data_baseline_dir):
    images = [os.path.join(dp, f) for dp, dn, filenames in os.walk(data_dir) for f in filenames if os.path.splitext(f)[1] == '.png']
    print(f'Found {len(images)} images')

    if not os.path.exists(data_baseline_dir):
        os.makedirs(data_baseline_dir)

    inputs = []
    for image in images:
        output_file = os.path.dirname(image).replace('\\', '-').replace('/', '-') + '.svg'
        output_file = os.path.join(data_baseline_dir, output_file)
        inputs.append([image, output_file])

    print(f'Launching {len(inputs)} jobs')
    manager = multiprocessing.Manager()
    tasks = manager.Queue()
    N_PROC = multiprocessing.cpu_count()        
    pool = multiprocessing.Pool(N_PROC)

    for input_dict in inputs:
        tasks.put(input_dict)
    for i in range(N_PROC):
        tasks.put(-1)

    processes = []
    for i in range(N_PROC):
        p = multiprocessing.Process(target = run_polyfit_multiple_inputs, args=(tasks,))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()

if __name__ == '__main__':
    run_data_baseline(POLYFIT_EXE, DATA_DIR, DATA_BASELINE_DIR)