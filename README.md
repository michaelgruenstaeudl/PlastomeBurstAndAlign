# PlastomeBurstAndAlign
Extract and align coding regions, introns and intergenic spacers across a set of plastomes

#### Installation
```
# Alignment software
apt install mafft  # Debian

# Other dependencies
pip install biopython
pip install coloredlogs
pip install ipdb
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

