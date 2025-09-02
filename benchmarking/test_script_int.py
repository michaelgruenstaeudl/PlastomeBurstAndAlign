#!/usr/bin/env python3

import glob
import os
import shutil
import sys
import time
import urllib.request
import tarfile
import subprocess  # Only used for running the module

# ----------------------------
# Setup paths and inputs
# ----------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
mod_rel = ".."
mod_path = os.path.join(script_dir, mod_rel)

dataset = sys.argv[1]
nCPUs = sys.argv[2]

main_folder = os.path.join(script_dir, dataset)
folder_INT = os.path.join(main_folder, dataset + "_output_INT")

# ----------------------------
# Step 1: Download dataset
# ----------------------------
url = f"https://zenodo.org/records/13403721/files/{dataset}.tar.gz?download=1"
out_file = os.path.join(script_dir, f"{dataset}.tar.gz")
print("Obtaining benchmarking dataset ...")
urllib.request.urlretrieve(url, out_file)

# ----------------------------
# Step 1b: Extract dataset
# ----------------------------
print("Extracting dataset ...")
with tarfile.open(out_file, "r:gz") as tar:
    tar.extractall(path=script_dir)

# Ensure main_folder exists
if not os.path.isdir(main_folder):
    # Try to detect extracted folder automatically
    extracted_dirs = [d for d in os.listdir(script_dir) if os.path.isdir(os.path.join(script_dir, d))]
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
    folder_INT,
    os.path.join(folder_INT, "01_unalign"),
    os.path.join(folder_INT, "02_aligned"),
    os.path.join(folder_INT, "02_aligned", "fasta"),
    os.path.join(folder_INT, "02_aligned", "nexus"),
]

for d in dirs_to_create:
    os.makedirs(d, exist_ok=True)

# ----------------------------
# Step 4: Run the module
# ----------------------------
print("Running plastburstalign module ...")
os.chdir(mod_path)
run_start = time.time()
subprocess.run([
    sys.executable, "-m", "plastburstalign",
    "-i", main_folder, "-o", folder_INT,
    "-s", "int", "-t", "5", "-l", "6", "-n", nCPUs
], check=True)
run_end = time.time()

# ----------------------------
# Step 5: Organize output files
# ----------------------------
def move_files(pattern, dest_folder):
    for f in glob.glob(pattern):
        shutil.move(f, dest_folder)

move_files(os.path.join(folder_INT, "*.unalign.fasta"), os.path.join(folder_INT, "01_unalign"))
move_files(os.path.join(folder_INT, "*.aligned.fasta"), os.path.join(folder_INT, "02_aligned", "fasta"))
move_files(os.path.join(folder_INT, "*.aligned.nexus"), os.path.join(folder_INT, "02_aligned", "nexus"))

# ----------------------------
# Step 6: Cleanup (optional)
# ----------------------------
os.chdir(script_dir)
# shutil.rmtree(dataset)
# print(f"Directory '{dataset}' and its contents have been deleted")

# ----------------------------
# Done
# ----------------------------
print("Time required for analysis: %.2f seconds" % (run_end - run_start))
