#### How to generate test data

##### Setting up benchmark
```bash
apt install ncbi-entrez-direct  # on Debian
cd benchmarking1
while read accn; do
  efetch -db nucleotide -id $accn -format gb > ${accn}.gb 2>>log.txt;
done <Accession_numbers_Asteraceae_YangEtAl2022.txt  # doi: 10.3389/fpls.2022.808156
```
