#!/usr/bin/env bash

# Benchmarking 1
#echo "ANALYSIS: benchmarking1_cds_n1" > benchmarking1_cds_n1.txt && python3 test_script_cds.py benchmarking1 1 1>>benchmarking1_cds_n1.txt 2>&1
#echo "ANALYSIS: benchmarking1_int_n1" > benchmarking1_int_n1.txt && python3 test_script_int.py benchmarking1 1 1>>benchmarking1_int_n1.txt 2>&1
#echo "ANALYSIS: benchmarking1_igs_n1" > benchmarking1_igs_n1.txt && python3 test_script_igs.py benchmarking1 1 1>>benchmarking1_igs_n1.txt 2>&1

#echo "ANALYSIS: benchmarking1_cds_n5" > benchmarking1_cds_n5.txt && python3 test_script_cds.py benchmarking1 5 1>>benchmarking1_cds_n5.txt 2>&1
#echo "ANALYSIS: benchmarking1_int_n5" > benchmarking1_int_n5.txt && python3 test_script_int.py benchmarking1 5 1>>benchmarking1_int_n5.txt 2>&1
#echo "ANALYSIS: benchmarking1_igs_n5" > benchmarking1_igs_n5.txt && python3 test_script_igs.py benchmarking1 5 1>>benchmarking1_igs_n5.txt 2>&1

echo "ANALYSIS: benchmarking1_cds_n10" > benchmarking1_cds_n10.txt && python3 test_script_cds.py benchmarking1 10 1>>benchmarking1_cds_n10.txt 2>&1
echo "ANALYSIS: benchmarking1_int_n10" > benchmarking1_int_n10.txt && python3 test_script_int.py benchmarking1 10 1>>benchmarking1_int_n10.txt 2>&1
echo "ANALYSIS: benchmarking1_igs_n10" > benchmarking1_igs_n10.txt && python3 test_script_igs.py benchmarking1 10 1>>benchmarking1_igs_n10.txt 2>&1

# Benchmarking 2
#echo "ANALYSIS: benchmarking2_cds_n1" > benchmarking2_cds_n1.txt && python3 test_script_cds.py benchmarking2 1 1>>benchmarking2_cds_n1.txt 2>&1
#echo "ANALYSIS: benchmarking2_int_n1" > benchmarking2_int_n1.txt && python3 test_script_int.py benchmarking2 1 1>>benchmarking2_int_n1.txt 2>&1
#echo "ANALYSIS: benchmarking2_igs_n1" > benchmarking2_igs_n1.txt && python3 test_script_igs.py benchmarking2 1 1>>benchmarking2_igs_n1.txt 2>&1

#echo "ANALYSIS: benchmarking2_cds_n5" > benchmarking2_cds_n5.txt && python3 test_script_cds.py benchmarking2 5 1>>benchmarking2_cds_n5.txt 2>&1
#echo "ANALYSIS: benchmarking2_int_n5" > benchmarking2_int_n5.txt && python3 test_script_int.py benchmarking2 5 1>>benchmarking2_int_n5.txt 2>&1
#echo "ANALYSIS: benchmarking2_igs_n5" > benchmarking2_igs_n5.txt && python3 test_script_igs.py benchmarking2 5 1>>benchmarking2_igs_n5.txt 2>&1

echo "ANALYSIS: benchmarking2_cds_n10" > benchmarking2_cds_n10.txt && python3 test_script_cds.py benchmarking2 10 1>>benchmarking2_cds_n10.txt 2>&1
echo "ANALYSIS: benchmarking2_int_n10" > benchmarking2_int_n10.txt && python3 test_script_int.py benchmarking2 10 1>>benchmarking2_int_n10.txt 2>&1
echo "ANALYSIS: benchmarking2_igs_n10" > benchmarking2_igs_n10.txt && python3 test_script_igs.py benchmarking2 10 1>>benchmarking2_igs_n10.txt 2>&1

# Benchmarking 3
#echo "ANALYSIS: benchmarking3_cds_n1" > benchmarking3_cds_n1.txt && python3 test_script_cds.py benchmarking3 1 1>>benchmarking3_cds_n1.txt 2>&1
#echo "ANALYSIS: benchmarking3_int_n1" > benchmarking3_int_n1.txt && python3 test_script_int.py benchmarking3 1 1>>benchmarking3_int_n1.txt 2>&1
#echo "ANALYSIS: benchmarking3_igs_n1" > benchmarking3_igs_n1.txt && python3 test_script_igs.py benchmarking3 1 1>>benchmarking3_igs_n1.txt 2>&1

#echo "ANALYSIS: benchmarking3_cds_n5" > benchmarking3_cds_n5.txt && python3 test_script_cds.py benchmarking3 5 1>>benchmarking3_cds_n5.txt 2>&1
#echo "ANALYSIS: benchmarking3_int_n5" > benchmarking3_int_n5.txt && python3 test_script_int.py benchmarking3 5 1>>benchmarking3_int_n5.txt 2>&1
#echo "ANALYSIS: benchmarking3_igs_n5" > benchmarking3_igs_n5.txt && python3 test_script_igs.py benchmarking3 5 1>>benchmarking3_igs_n5.txt 2>&1

echo "ANALYSIS: benchmarking3_cds_n10" > benchmarking3_cds_n10.txt && python3 test_script_cds.py benchmarking3 10 1>>benchmarking3_cds_n10.txt 2>&1
echo "ANALYSIS: benchmarking3_int_n10" > benchmarking3_int_n10.txt && python3 test_script_int.py benchmarking3 10 1>>benchmarking3_int_n10.txt 2>&1
echo "ANALYSIS: benchmarking3_igs_n10" > benchmarking3_igs_n10.txt && python3 test_script_igs.py benchmarking3 10 1>>benchmarking3_igs_n10.txt 2>&1


# Parse data and generate output
OUTF="Figure_SpeedImprov_InputData.csv"
echo "DATASET,MARKER,CPUS,NTAX,TIME" > $OUTF
TMPF1=$(mktemp)
echo "" > $TMPF1
TMPF2=$(mktemp)
echo "" > $TMPF2

for i in $(ls benchmarking?_*_n*.txt); do
  grep "ANALYSIS: " $i | tr -s '_' ',' >> $TMPF1
  grep "Size of input dataset: " $i | sed 's/Size of input dataset: //g' | awk '{print $1}' >> $TMPF1
  grep "Time required for analysis: " $i | sed 's/Time required for analysis: //g' | awk '{print $1}' >> $TMPF1;
done

cat $TMPF1 | awk -v ORS="," '1' | sed 's/,ANALYSIS: /\n/g' >> $TMPF2
cat $TMPF2 | sed '/^[[:space:]]*$/d' | sed 's/,$/$/g' >> $OUTF
rm $TMPF1 $TMPF2
