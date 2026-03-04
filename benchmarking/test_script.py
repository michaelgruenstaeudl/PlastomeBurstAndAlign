#!/usr/bin/env python3

import glob
import os
import shutil
import sys
import time
import urllib.request
import tarfile
import subprocess
import pandas as pd
import psutil
from datetime import datetime
import resource
import json

# def get_memory_usage():
#     """Get current memory usage in MB"""
#     process = psutil.Process(os.getpid())
#     return process.memory_info().rss / 1024 / 1024  # MB

def run_single_test(dataset, nCPUs, mode):
    """Run a single test case and return performance metrics"""

    script_dir = os.path.dirname(os.path.abspath(__file__))
    mod_rel = ".."
    mod_path = os.path.join(script_dir, mod_rel)

    main_folder = os.path.join(script_dir, dataset)

    ts = datetime.now().strftime("%Y%m%dT%H%M")
    folder = os.path.join(main_folder, f"{os.path.basename(dataset)}_output_{mode}_{ts}")

    # Record initial memory
    # start_memory = get_memory_usage()
    start_time = time.time()

    try:
        # ----------------------------
        # Step 1: Download dataset (if not already present)
        # ----------------------------
        dataset_file = os.path.join(script_dir, f"{dataset}.tar.gz")
        print(f"Preparing dataset...{os.path.exists(main_folder) }")
        if not os.path.exists(main_folder) and not os.path.exists(dataset_file):
            print(f"Downloading dataset {dataset}...")
            url = f"https://zenodo.org/records/17137737/files/{dataset}.tar.gz?download=1"
            urllib.request.urlretrieve(url, dataset_file)

        # ----------------------------
        # Step 1b: Extract dataset (if not already extracted)
        # ----------------------------
        if not os.path.isdir(main_folder) and os.path.exists(dataset_file):
            print("Extracting dataset ...")
            with tarfile.open(dataset_file, "r:gz") as tar:
                tar.extractall(path=script_dir)

        # Ensure main_folder exists
        if not os.path.isdir(main_folder):
            extracted_dirs = [d for d in os.listdir(script_dir) if os.path.isdir(os.path.join(script_dir, d)) and dataset in d]
            if extracted_dirs:
                main_folder = os.path.join(script_dir, extracted_dirs[0])
        print(f"Using dataset folder: {main_folder}")

        # ----------------------------
        # Step 2: Count GenBank files
        # ----------------------------
        gb_files = glob.glob(os.path.join(main_folder, "*.gb"))
        print(f"Size of input dataset: {len(gb_files)} GenBank files\n")

        # ----------------------------
        # Step 3: Create necessary directories
        # ----------------------------
        dirs_to_create = [
            folder,
            os.path.join(folder, "01_unalign"),
            os.path.join(folder, "02_aligned"),
            os.path.join(folder, "02_aligned", "fasta"),
            os.path.join(folder, "02_aligned", "nexus"),
        ]

        for d in dirs_to_create:
            os.makedirs(d, exist_ok=True)

        # ----------------------------
        # Step 4: Run the module with resource monitoring
        # ----------------------------
        print(f"Running plastburstalign with {nCPUs} threads, mode: {mode}...")
        os.chdir(mod_path)

        # Set environment for thread control
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(nCPUs)
        env['OPENBLAS_NUM_THREADS'] = str(nCPUs)
        env['MKL_NUM_THREADS'] = str(nCPUs)

        run_start = time.time()
        process = subprocess.run([
            sys.executable, "-m", "plastburstalign",
            "-i", main_folder, "-o", folder,
            "-s", mode,  "-t", "5", "-l", "9","-n", str(nCPUs)
        ], check=True)
        run_end = time.time()

        # ----------------------------
        # Step 5: Organize output files
        # ----------------------------
        def move_files(pattern, dest_folder):
            for f in glob.glob(pattern):
                shutil.move(f, dest_folder)

        move_files(os.path.join(folder, "*.unalign.fasta"), os.path.join(folder, "01_unalign"))
        move_files(os.path.join(folder, "*.aligned.fasta"), os.path.join(folder, "02_aligned", "fasta"))
        move_files(os.path.join(folder, "*.aligned.nexus"), os.path.join(folder, "02_aligned", "nexus"))

        # Record final memory and time
        # end_memory = get_memory_usage()
        end_time = time.time()

        # Calculate metrics
        total_time = run_end - run_start
        
        # max_memory_used = max(get_memory_usage() for _ in range(10))  # Sample memory multiple times
        # memory_increase = end_memory - start_memory

        # Return results
        return {
            'dataset': dataset,
            'threads': nCPUs,
            'mode': mode,
            'genbank_files': len(gb_files),
            'total_time_seconds': total_time,
            # 'max_memory_mb': max_memory_used,
            # 'memory_increase_mb': memory_increase,
            'throughput_files_per_sec': len(gb_files) / total_time if total_time > 0 else 0,
            'success': True,
            'error': None,
            'timestamp': ts
        }

    except Exception as e:
        end_time = time.time()
        return {
            'dataset': dataset,
            'threads': nCPUs,
            'mode': mode,
            'total_time_seconds': time.time() - start_time,
            
            # 'max_memory_mb': get_memory_usage(),
            # 'memory_increase_mb': get_memory_usage() - start_memory,
            'throughput_files_per_sec': 0,
            'success': False,
            'error': str(e),
            'timestamp': ts
        }

def test(dataset=None, nCPUs=None, mode=None):
    """Main benchmark runner"""

    # Test configurations
   

    # Results storage
    all_results = []

    
    
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_time")
    os.makedirs(results_dir, exist_ok=True)

    
    result = run_single_test(dataset, nCPUs, mode)
    all_results.append(result)


    df = pd.DataFrame(all_results)
    results_file = os.path.join(results_dir, f"plastburst_benchmark.csv")


    if os.path.exists(results_file):
        df.to_csv(results_file, mode='a', header=False, index=False)
    else:
        df.to_csv(results_file, index=False)

                
    



if __name__ == "__main__":
    dataset = sys.argv[1]
    nCPUs = sys.argv[2]
    mode = sys.argv[3]
    
    test(dataset, nCPUs, mode)
