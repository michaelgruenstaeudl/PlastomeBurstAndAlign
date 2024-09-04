#!/usr/bin/env bash

# Benchmarking 1
echo "ANALYSIS: benchmarking1_cds_n1" > benchmarking1_cds_n1.txt && python3 test_script_cds.py benchmarking1 1 1>>benchmarking1_cds_n1.txt 2>&1
echo "ANALYSIS: benchmarking1_int_n1" > benchmarking1_int_n1.txt && python3 test_script_int.py benchmarking1 1 1>>benchmarking1_int_n1.txt 2>&1
echo "ANALYSIS: benchmarking1_igs_n1" > benchmarking1_igs_n1.txt && python3 test_script_igs.py benchmarking1 1 1>>benchmarking1_igs_n1.txt 2>&1

echo "ANALYSIS: benchmarking1_cds_n5" > benchmarking1_cds_n5.txt && python3 test_script_cds.py benchmarking1 5 1>>benchmarking1_cds_n5.txt 2>&1
echo "ANALYSIS: benchmarking1_int_n5" > benchmarking1_int_n5.txt && python3 test_script_int.py benchmarking1 5 1>>benchmarking1_int_n5.txt 2>&1
echo "ANALYSIS: benchmarking1_igs_n5" > benchmarking1_igs_n5.txt && python3 test_script_igs.py benchmarking1 5 1>>benchmarking1_igs_n5.txt 2>&1

echo "ANALYSIS: benchmarking1_cds_n10" > benchmarking1_cds_n10.txt && python3 test_script_cds.py benchmarking1 10 1>>benchmarking1_cds_n10.txt 2>&1
echo "ANALYSIS: benchmarking1_int_n10" > benchmarking1_int_n10.txt && python3 test_script_int.py benchmarking1 10 1>>benchmarking1_int_n10.txt 2>&1
echo "ANALYSIS: benchmarking1_igs_n10" > benchmarking1_igs_n10.txt && python3 test_script_igs.py benchmarking1 10 1>>benchmarking1_igs_n10.txt 2>&1

# Benchmarking 2
echo "ANALYSIS: benchmarking2_cds_n1" > benchmarking2_cds_n1.txt && python3 test_script_cds.py benchmarking2 1 1>>benchmarking2_cds_n1.txt 2>&1
echo "ANALYSIS: benchmarking2_int_n1" > benchmarking2_int_n1.txt && python3 test_script_int.py benchmarking2 1 1>>benchmarking2_int_n1.txt 2>&1
echo "ANALYSIS: benchmarking2_igs_n1" > benchmarking2_igs_n1.txt && python3 test_script_igs.py benchmarking2 1 1>>benchmarking2_igs_n1.txt 2>&1

echo "ANALYSIS: benchmarking2_cds_n5" > benchmarking2_cds_n5.txt && python3 test_script_cds.py benchmarking2 5 1>>benchmarking2_cds_n5.txt 2>&1
echo "ANALYSIS: benchmarking2_int_n5" > benchmarking2_int_n5.txt && python3 test_script_int.py benchmarking2 5 1>>benchmarking2_int_n5.txt 2>&1
echo "ANALYSIS: benchmarking2_igs_n5" > benchmarking2_igs_n5.txt && python3 test_script_igs.py benchmarking2 5 1>>benchmarking2_igs_n5.txt 2>&1

echo "ANALYSIS: benchmarking2_cds_n10" > benchmarking2_cds_n10.txt && python3 test_script_cds.py benchmarking2 10 1>>benchmarking2_cds_n10.txt 2>&1
echo "ANALYSIS: benchmarking2_int_n10" > benchmarking2_int_n10.txt && python3 test_script_int.py benchmarking2 10 1>>benchmarking2_int_n10.txt 2>&1
echo "ANALYSIS: benchmarking2_igs_n10" > benchmarking2_igs_n10.txt && python3 test_script_igs.py benchmarking2 10 1>>benchmarking2_igs_n10.txt 2>&1

# Benchmarking 3
echo "ANALYSIS: benchmarking3_cds_n1" > benchmarking3_cds_n1.txt && python3 test_script_cds.py benchmarking3 1 1>>benchmarking3_cds_n1.txt 2>&1
echo "ANALYSIS: benchmarking3_int_n1" > benchmarking3_int_n1.txt && python3 test_script_int.py benchmarking3 1 1>>benchmarking3_int_n1.txt 2>&1
echo "ANALYSIS: benchmarking3_igs_n1" > benchmarking3_igs_n1.txt && python3 test_script_igs.py benchmarking3 1 1>>benchmarking3_igs_n1.txt 2>&1

echo "ANALYSIS: benchmarking3_cds_n5" > benchmarking3_cds_n5.txt && python3 test_script_cds.py benchmarking3 5 1>>benchmarking3_cds_n5.txt 2>&1
echo "ANALYSIS: benchmarking3_int_n5" > benchmarking3_int_n5.txt && python3 test_script_int.py benchmarking3 5 1>>benchmarking3_int_n5.txt 2>&1
echo "ANALYSIS: benchmarking3_igs_n5" > benchmarking3_igs_n5.txt && python3 test_script_igs.py benchmarking3 5 1>>benchmarking3_igs_n5.txt 2>&1

echo "ANALYSIS: benchmarking3_cds_n10" > benchmarking3_cds_n10.txt && python3 test_script_cds.py benchmarking3 10 1>>benchmarking3_cds_n10.txt 2>&1
echo "ANALYSIS: benchmarking3_int_n10" > benchmarking3_int_n10.txt && python3 test_script_int.py benchmarking3 10 1>>benchmarking3_int_n10.txt 2>&1
echo "ANALYSIS: benchmarking3_igs_n10" > benchmarking3_igs_n10.txt && python3 test_script_igs.py benchmarking3 10 1>>benchmarking3_igs_n10.txt 2>&1


# Generate output
OUTF1="input_data_DATEHERE.tmp"
OUTF2="input_data_DATEHERE.csv"
echo "" > touch $OUTF1
echo "DATASET,MARKER,CPUS,TIME" > touch $OUTF2

for i in $(ls benchmarking?_*_n*.txt); do
  grep "ANALYSIS: " $i >> $OUTF1
  grep "Size of input dataset: " $i | sed 's/Size of input dataset: //g' >> $OUTF1
  grep "Time required for analysis: " $i | sed 's/Time required for analysis: //g' >> $OUTF1;
done

tr -s '\n' ','<$OUTF1 >tmp.txt
tr -s ',ANALYSIS: ' '\n'<tmp.txt >$OUTF2

