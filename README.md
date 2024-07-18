# PlastomeBurstAndAlign
A Python tool to extract and align genes, introns, and intergenic spacers across hundreds of plastid genomes using associative arrays

#### Installation on Linux (Debian)
```
# Alignment software
apt install mafft

# Other dependencies
apt install python3-biopython
apt install python3-coloredlogs
```

#### Overview of process
![](docs/PlastomeBurstAndAlign_ProcessOverview.png)

#### Main features
- Extracts all genes, introns, or intergenic spacers from a set of plastid genomes in GenBank flatfile format, aligns the corresponding regions, and then saves the alignments
- Saves both the marker-wise alignments and the concatenation of all alignments

#### Additional features
- The sequence extraction and alignment of the genes/introns/intergenic spacers is parallelized and, thus, sped up by utilizing multiple CPUs.
- The order of concatenation of the aligned genes/introns/intergenic spacers can be determined (via command-line argument `--order/-r`): either corresponding to the natural sequence order of the first input genome (option `seq`) or in alphabetic order (option `alpha`).
- Even though GenBank flatfiles list trans-spliced genes (e.g., rps12) out of order (i.e., out of their natural order along the genome sequence), trans-spliced genes are repositioned and their exons automatically merged.
- The sofwtare automatically adjusts for the majority of letter-case differences between gene annotations of different plastid genome records.

#### Testing
```
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

#### Exemplary usage
See [this document](docs/exemplary_usage.md)


#### Generating more test data
See [this document](docs/generating_test_data.md)

