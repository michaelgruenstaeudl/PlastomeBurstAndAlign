# plastburstalign
A Python tool to extract and align genes, introns, and intergenic spacers across hundreds of plastid genomes using associative arrays

### Installation on Linux (Debian)
```bash
# Alignment software
apt install mafft

# Other dependencies
apt install python3-biopython
apt install python3-coloredlogs
```

### Overview of process
![Depiction of plastomes being split according to specified feature type; these features and then aligned by region](docs/PlastomeBurstAndAlign_ProcessOverview.png)

### Main features
- Extracts all genes, introns, or intergenic spacers from a set of plastid genomes in GenBank flatfile format, aligns the corresponding regions, and then saves the alignments
- Saves both the marker-wise alignments and the concatenation of all alignments

### Additional features
- The sequence extraction and alignment of the genes/introns/intergenic spacers is parallelized and, thus, sped up by utilizing multiple CPUs.
- The order of concatenation of the aligned genes/introns/intergenic spacers can be determined (via command-line argument `--order/-r`): either corresponding to the natural sequence order of the first input genome (option `seq`) or in alphabetic order (option `alpha`).
- Even though GenBank flatfiles list trans-spliced genes (e.g., rps12) out of order (i.e., out of their natural order along the genome sequence), trans-spliced genes are repositioned and their exons automatically merged.
- The software automatically adjusts for the majority of letter-case differences between gene annotations of different plastid genome records.

### Usage
There are two main ways this package can be used, depending on one's use case.

#### Package as a script
With the package installed or with the current working directory at the same level as `plastburstalign` (relative path will not work) you can execute the package with something like
```bash
python -m plastburstalign
```
which will execute the pipeline with default arguments.

#### Package as a module
You are also able to use the objects defined in the package to perform the pipeline within Python such as

```python
from plastburstalign import PlastomeRegionBurstAndAlign

burst = PlastomeRegionBurstAndAlign()
burst.execute()
```
which will perform the same action as the above command line use. That is, an execution of the pipeline with default arguments.

You can also use individual components for your own use. For example, the `MAFFT` class can be used by itself like

```python
from plastburstalign import MAFFT

mafft_1 = MAFFT()
mafft_10 = MAFFT({"num_threads": 10})
```
which will instantiate a configuration of MAFFT that will execute with 1 thread and another that will execute with 10 threads.

### Testing
```bash
cd benchmarking
# CDS
python test_script_cds.py benchmarking1
python test_script_cds.py benchmarking2
# INT
python test_script_int.py benchmarking1
python test_script_int.py benchmarking2
# IGS
python test_script_igs.py benchmarking1
python test_script_igs.py benchmarking2
```
- Dataset `benchmarking1.tar.gz`: all Asteraceae (n=155) listed in [Yang et al. 2022](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2022.808156)
- Dataset `benchmarking2.tar.gz`: all monocots (n=733) listed in [Yang et al. 2022](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2022.808156)

### Exemplary usage
See [this document](https://github.com/michaelgruenstaeudl/PlastomeBurstAndAlign/blob/main/docs/exemplary_usage.md)


### Generating more test data
See [this document](https://github.com/michaelgruenstaeudl/PlastomeBurstAndAlign/blob/main/docs/generating_test_data.md)

