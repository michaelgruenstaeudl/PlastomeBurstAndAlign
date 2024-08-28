# plastburstalign
A Python tool to extract and align genes, introns, and intergenic spacers across hundreds of plastid genomes using associative arrays

### Background
The multiple sequence alignment (MSA) of sets of complete plastid genomes is challenging. At least five factors are responsible for this challenge:
- First, the plastid genome is a mosaic of individual genome regions, and a MSA procedure must extract, group, and align those regions across the genomes that are homologous to each other.
- Second, many plastid genomes exhibit sequence annotation errors regarding gene position and/or gene name, and a MSA procedure must automatically remove the regions of those genomes that exhibit such errors.
- Third, plastid genomes comprise both coding and noncoding genome regions, and a MSA procedure must automatically employ the best employment mechanisms for each (i.e., amino acid-based alignment for coding regions, nucleotide-based alignment for noncoding regions).
- Fourth, contemporary plastid genome investigations comprise hundreds of complete plastid genomes, and a MSA procedure must perform sequence alignment within reasonable time frames (i.e., hours instead of days)
- Fifth, any user-specified exclusion of a genome region from the alignment would be prohibitively complex ater the alignment and must, conseuently, be part of the MSA procedure

plastburstalign accommodates these five (and more) challenges and extracts and aligns genes, introns, and intergenic spacers across hundreds or thousands of input plastid genomes.

### Overview of process
![Depiction of plastomes being split according to specified marker type; the extracted sequences are then aligned and concatenated](docs/PlastomeBurstAndAlign_ProcessOverview.png)

### Main features
- Extraction of all genome regions representing one of three different marker types (i.e., genes, introns, or intergenic spacers) from set of input plastid genomes, followed by grouping and alignment of the extracted regions
- Automatic exon splicing:
  - Automatic merging of all exons of any cis-spliced gene
  - Automatic grouping of all exons of any trans-spliced gene (e.g., _rps12_), followed by merging of adjacent exons [see `ExonSpliceHandler` for both]
- Automatic quality control to evaluate if extracted genes are complete (i.e., valid start and stop codon present)
- Automatic removal of
  - any duplicate regions (i.e., relevant for regions duplicated through the IRs)
  - regions that do not fulfill a minimum, user-specified sequence length
  - regions that do not fulfill a minimum, user-specified number of taxa of the dataset that the region must be found in [see `DataCleaning` for both]
  - any user-specified genome region (i.e., gene, intron, or intergenic spacer)
- For genes: DNA sequence alignment based on amino acid sequences of genes rather than nucleotide sequences
- Rapid DNA sequence extraction and alignment due to process parallelization using multiple CPUs [see `_nuc_MSA()`]

### Additional features
- Choice of
  - the order of concatenation of the aligned genes/introns/intergenic spacers to either the natural order of the first input genome (commandline option `seq`) or an alphabetic order (commandline option `alpha`)
  - automatic case standardization of gene names to adjust for letter-case differences between gene annotations of different genome records (which is especially relevant for anticodon and amino acid abbreviations of tRNAs); includes the option to remove anticodon and amino acid abbreviations from tRNA gene names altogether [see `clean_gene()`]
- If a gene/intron/intergenic spacer cannot be extracted from a GenBank record, provision of explanation why the extraction failed
- Availability of two log levels:
  - default (suitable for regular software execution), and 
  - verbose (suitable for debugging)
- Package works out of the box on Unix-like systems due to inclusion of the alignment software executable (MAFFT) into the package.

### Input/output
**Input**: Complete plastid genomes in GenBank flatfile format
**Output**:  Marker-wise alignments and concatenation of all alignments, both in FASTA format

### Installation on Linux (Debian)
```bash
# Alignment software
apt install mafft

# Other dependencies
apt install python3-biopython
apt install python3-coloredlogs
apt install python3-requests
```

### Usage

#### Option 1: As a script
If current working directory within `plastburstalign`, execute the package via:
```bash
python -m plastburstalign
```

#### Option 2: As a module
From within Python, execute the package functions via:
```python
from plastburstalign import PlastomeRegionBurstAndAlign
burst = PlastomeRegionBurstAndAlign()
burst.execute()
```

#### Usage of individual package components
Individual components can be used as well. For example, to use the class `MAFFT` by itself (e.g., instantiate a configuration of MAFFT that will execute with 1 thread; institute another that will execute with 10 threads), type:

```python
from plastburstalign import MAFFT

mafft_1 = MAFFT()
mafft_10 = MAFFT({"num_threads": 10})
```


### Explanation of exon splicing
As the gene list produced through parsing all input genomes is iterated over, genes that comprise multiple exons are automatically flagged and treated according to the distance between their exons. Cis-spliced genes only comprise exons that are adjacent to each other, trans-spliced genes comprise one or more exons that are not adjacent to each other. This software merges the exons of any cis-spliced gene in place (i.e., according to the location specified by the source GenBank file; no repositioning of the exons necessary). The exons of any trans-spliced gene (e.g., _rps12_), by contrast, undergo a repositioning before being merged. Specifically, the software accommodates the fact that GenBank flatfiles list trans-spliced genes (e.g., _rps12_) out of their natural order along the genome sequence and additionally repositions the exons of trans-spliced genes by converting them to adjacent exons and then merges these exons.

For the repositioning of trans-spliced gene, all annotations of that gene are first moved from the main gene list to a separate list. Then, the annotations are split into simple location features for each contiguous group of exons. Third, the expected location of each of these simple gene features is determined by comparing its end location with the end locations of the gene features in the main gene list: if the expected location has no overlap with either the proceeding and succeeding genes and the feature is different in name from either, it is directly inserted into that location. Alternatively, if the expected location of the feature results in a flanking gene (strictly adjacent or overlapping) with the same name, the annotations are merged; the merging is true for both the proceeding and the succeeding gene.


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

