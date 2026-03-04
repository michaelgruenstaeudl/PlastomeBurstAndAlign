#!/usr/bin/env python3

import re
import sys
from typing import List, Tuple, Dict, Set
from collections import defaultdict

def extract_accession_fasta(header: str) -> str:
    """
    Extract accession number from FASTA header.
    For headers like: >accD_NC_034851
    """
    if header.startswith('>'):
        header = header[1:].strip()
    
    # Extract the NC_ number part
    match = re.search(r'(accD_NC_\d+)', header)
    if match:
        return match.group(1)
    
    # Fallback: use the entire header (without >)
    return header

def read_fasta_sequences(filename: str) -> Dict[str, str]:
    """
    Read FASTA file and return dictionary of {accession: sequence}
    """
    sequences = {}
    current_header = None
    current_sequence = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header and current_sequence:
                    accession = extract_accession_fasta(current_header)
                    sequences[accession] = ''.join(current_sequence)
                
                # Start new sequence
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Don't forget the last sequence
        if current_header and current_sequence:
            accession = extract_accession_fasta(current_header)
            sequences[accession] = ''.join(current_sequence)
    
    return sequences

def compare_sequences(seq1: str, seq2: str) -> Dict[str, any]:
    """
    Compare two sequences and return detailed differences.
    """
    len1, len2 = len(seq1), len(seq2)
    
    # Count mismatches and gaps
    mismatches = 0
    gaps_seq1 = seq1.count('-')
    gaps_seq2 = seq2.count('-')
    
    # Compare position by position (up to the shorter length)
    min_len = min(len1, len2)
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            mismatches += 1
    
    # Count length difference as additional mismatches
    mismatches += abs(len1 - len2)
    
    return {
        'length1': len1,
        'length2': len2,
        'length_difference': abs(len1 - len2),
        'mismatches': mismatches,
        'similarity_percent': (1 - mismatches / max(len1, len2)) * 100 if max(len1, len2) > 0 else 0,
        'gaps_seq1': gaps_seq1,
        'gaps_seq2': gaps_seq2
    }

def compare_fasta_files(file1: str, file2: str, output_file: str = None, detailed: bool = False):
    """
    Compare two FASTA files by sorting and analyzing sequences based on accession numbers.
    """
    print(f"Reading {file1}...")
    sequences1 = read_fasta_sequences(file1)
    
    print(f"Reading {file2}...")
    sequences2 = read_fasta_sequences(file2)
    
    # Get all unique accessions (sorted)
    all_accessions = sorted(set(sequences1.keys()) | set(sequences2.keys()))
    
    # Categorize sequences
    only_in_file1 = []
    only_in_file2 = []
    common_sequences = []
    
    for accession in all_accessions:
        in_file1 = accession in sequences1
        in_file2 = accession in sequences2
        
        if in_file1 and not in_file2:
            only_in_file1.append(accession)
        elif in_file2 and not in_file1:
            only_in_file2.append(accession)
        else:
            comparison = compare_sequences(sequences1[accession], sequences2[accession])
            common_sequences.append((accession, sequences1[accession], sequences2[accession], comparison))
    
    # Print summary
    print("\n" + "="*80)
    print("FASTA FILE COMPARISON RESULTS")
    print("="*80)
    print(f"File 1: {file1}")
    print(f"File 2: {file2}")
    print(f"Total unique accessions: {len(all_accessions)}")
    print(f"Only in {file1}: {len(only_in_file1)}")
    print(f"Only in {file2}: {len(only_in_file2)}")
    print(f"Common accessions: {len(common_sequences)}")
    
    # Analyze common sequences
    if common_sequences:
        identical_count = sum(1 for _, seq1, seq2, _ in common_sequences if seq1 == seq2)
        different_count = len(common_sequences) - identical_count
        
        print(f"  - Identical sequences: {identical_count}")
        print(f"  - Different sequences: {different_count}")
        
        if different_count > 0:
            avg_similarity = sum(comp['similarity_percent'] for _, _, _, comp in common_sequences) / len(common_sequences)
            max_mismatches = max(comp['mismatches'] for _, _, _, comp in common_sequences)
            print(f"  - Average similarity: {avg_similarity:.2f}%")
            print(f"  - Maximum mismatches: {max_mismatches}")
    
    # Detailed output for sequences only in one file
    if only_in_file1:
        print(f"\n--- Sequences only in {file1} ({len(only_in_file1)}) ---")
        for i, acc in enumerate(only_in_file1[:5]):
            seq = sequences1[acc]
            print(f"  {i+1}. {acc} (length: {len(seq)})")
        if len(only_in_file1) > 5:
            print(f"  ... and {len(only_in_file1) - 5} more")
    
    if only_in_file2:
        print(f"\n--- Sequences only in {file2} ({len(only_in_file2)}) ---")
        for i, acc in enumerate(only_in_file2[:5]):
            seq = sequences2[acc]
            print(f"  {i+1}. {acc} (length: {len(seq)})")
        if len(only_in_file2) > 5:
            print(f"  ... and {len(only_in_file2) - 5} more")
    
    # Show most different sequences
    if common_sequences and different_count > 0:
        print(f"\n--- Most Different Common Sequences (Top 5) ---")
        # Sort by number of mismatches (descending)
        sorted_different = sorted(common_sequences, key=lambda x: x[3]['mismatches'], reverse=True)[:5]
        
        for acc, seq1, seq2, comp in sorted_different:
            if seq1 != seq2:  # Only show actually different ones
                print(f"  {acc}:")
                print(f"    Similarity: {comp['similarity_percent']:.1f}%")
                print(f"    Mismatches: {comp['mismatches']}")
                print(f"    Lengths: {comp['length1']} vs {comp['length2']} (diff: {comp['length_difference']})")
    
    # Write detailed report if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write("FASTA COMPARISON REPORT\n")
            f.write("="*60 + "\n")
            f.write(f"File 1: {file1}\n")
            f.write(f"File 2: {file2}\n\n")
            
            f.write("SUMMARY\n")
            f.write(f"Total unique accessions: {len(all_accessions)}\n")
            f.write(f"Only in {file1}: {len(only_in_file1)}\n")
            f.write(f"Only in {file2}: {len(only_in_file2)}\n")
            f.write(f"Common accessions: {len(common_sequences)}\n")
            
            if common_sequences:
                identical_count = sum(1 for _, seq1, seq2, _ in common_sequences if seq1 == seq2)
                f.write(f"  - Identical: {identical_count}\n")
                f.write(f"  - Different: {len(common_sequences) - identical_count}\n")
            
            # Write sequences only in file1
            f.write("\nSEQUENCES ONLY IN FILE 1\n")
            for acc in sorted(only_in_file1):
                seq = sequences1[acc]
                f.write(f">{acc}\n")
                f.write(f"{seq}\n\n")
            
            # Write sequences only in file2
            f.write("\nSEQUENCES ONLY IN FILE 2\n")
            for acc in sorted(only_in_file2):
                seq = sequences2[acc]
                f.write(f">{acc}\n")
                f.write(f"{seq}\n\n")
            
            # Write different common sequences
            f.write("\nCOMMON SEQUENCES WITH DIFFERENCES\n")
            for acc, seq1, seq2, comp in common_sequences:
                if seq1 != seq2:
                    f.write(f">{acc} | Similarity: {comp['similarity_percent']:.1f}% | Mismatches: {comp['mismatches']}\n")
                    f.write(f"FILE1: {seq1}\n")
                    f.write(f"FILE2: {seq2}\n\n")
        
        print(f"\nDetailed report written to: {output_file}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_fasta.py file1.fasta file2.fasta [output_report.txt]")
        print("Example: python compare_fasta.py alignment1.fasta alignment2.fasta comparison.txt")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    compare_fasta_files(file1, file2, output_file)

if __name__ == "__main__":
    main()
