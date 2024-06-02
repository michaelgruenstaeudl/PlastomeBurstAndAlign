# PlastomeBurstAndAlign
A Python tool to extract and align genes, introns, and intergenic spacers across plastid genomes using associative arrays

#### Installation on Linux (Debian)
```
# Alignment software
apt install mafft

# Other dependencies
apt install python3-biopython
apt install python3-coloredlogs
apt install python3-ipdb
```

#### Overview of process
![](docs/PlastomeBurstAndAlign_ProcessOverview.png)

#### Testing
```
cd benchmarking
# CDS
python test_script_cds.py
# INT
python test_script_int.py
# IGS
python test_script_igs.py
```

#### Exemplary usage
See [this document](docs/exemplary_usage.md)


#### Generating more test data
See [this document](docs/generating_test_data.md)

