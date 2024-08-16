#!/usr/bin/env python3

import subprocess
import os
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

# Define the path to the input and output directory
main_folder = os.path.join(script_dir, dataset)
folder_IGS = os.path.join(main_folder, "output_IGS")

# Step 1: Extract the tar.gz file
subprocess.run(["tar", "-xf", os.path.join(script_dir, dataset+".tar.gz"), "-C", script_dir])

# Step 2: Change the current directory to the selected dataset
os.chdir(main_folder)
subprocess.call('echo "Size of input dataset: $(ls *.gb 2>\\dev\\null | wc -l) GenBank files"', shell=True)
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
                "-i", main_folder, "-o", folder_IGS, "-s", "igs", "-t", "5", "-l", "6", "-n", "10"])
run_end = time.time()

# run this to remove the folder, if not can comment out
# Step 5: Delete the dataset directory and its contents
os.chdir(script_dir)  # Move back to the directory containing the test script
shutil.rmtree(dataset)
print("Directory '" + dataset + "' and its contents have been deleted")
print("Time required for analysis: %.2f seconds" % (run_end - run_start))

