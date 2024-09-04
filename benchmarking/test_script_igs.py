#!/usr/bin/env python3

import glob
import os
import subprocess
import shutil
import sys
import time

# Get the current directory of the test script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define the absolute path to the module
mod_rel = ".."
mod_path = os.path.join(script_dir, mod_rel)

# Define the benchmarking dataset
dataset = sys.argv[1]
nCPUs = sys.argv[2]

# Define the path to the input and output directory
main_folder = os.path.join(script_dir, dataset)
folder_IGS = os.path.join(main_folder, dataset+"_output_IGS")

# Step 1: Obtain the benchmarking dataset from Zenodo and decompress it
subprocess.call('echo "Obtaining benchmarking dataset"', shell=True)
subprocess.run(["wget", "-q", "-c", "-O", ''.join([dataset, '.tar.gz']), ''.join(['https://zenodo.org/records/13403721/files/', dataset, '.tar.gz?download=1'])])
subprocess.run(["tar", "-xf", os.path.join(script_dir, dataset+".tar.gz"), "-C", script_dir])

# Step 2: Change the current directory to the selected dataset
os.chdir(main_folder)
subprocess.call('echo "Size of input dataset: $(ls *.gb 2>/dev/null | wc -l) GenBank files"', shell=True)
subprocess.call('echo ""', shell=True)

# Step 3: Create necessary directories
subprocess.run(["mkdir", "-p", folder_IGS])
subprocess.run(["mkdir", "-p", f"{folder_IGS}/01_unalign"])
subprocess.run(["mkdir", "-p", f"{folder_IGS}/02_aligned"])
subprocess.run(["mkdir", "-p", f"{folder_IGS}/02_aligned/fasta"])
subprocess.run(["mkdir", "-p", f"{folder_IGS}/02_aligned/nexus"])

# Step 4: Change the directory to the module and run it as a script
os.chdir(mod_path)
run_start = time.time()
subprocess.run(["python3", "-m", "plastburstalign",
                "-i", main_folder, "-o", folder_IGS, "-s", "igs", "-t", "5", "-l", "6", "-n", nCPUs])
run_end = time.time()

# Step 5: Organize the output files
subprocess.run(["mv"] + glob.glob(f"{folder_IGS}/*.unalign.fasta") + [f"{folder_IGS}/01_unalign/"])
subprocess.run(["mv"] + glob.glob(f"{folder_IGS}/*.aligned.fasta") + [f"{folder_IGS}/02_aligned/fasta/"])
subprocess.run(["mv"] + glob.glob(f"{folder_IGS}/*.aligned.nexus") + [f"{folder_IGS}/02_aligned/nexus/"])

# Step 6: Delete the dataset directory and its contents
os.chdir(script_dir)  # Move back to the directory containing the test script
#shutil.rmtree(dataset)
#print("Directory '" + dataset + "' and its contents have been deleted")
print("Time required for analysis: %.2f seconds" % (run_end - run_start))

