"""
bioseq.py - A lightweight DNA sequence analysis toolkit.

Features:
    - Sequence validation
    - GC content
    - Reverse complement
    - Transcription (DNA -> mRNA)
    - Translation (mRNA/DNA -> protein)
    - ORF (open reading frame) finder on both strands
    - Primer melting temperature (Wallace rule and basic Tm formula)

No external dependencies required (Python 3.8+).
"""

import argparse
import sys
from typing import List, Tuple

VALID_BASES = set("ATGC")

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def clean_sequence(seq: str) -> str:
    """Uppercase, strip whitespace, and validate a DNA sequence."""
    seq = "".join(seq.split()).upper()
    invalid = set(seq) - VALID_BASES
    if invalid:
        raise ValueError(f"Invalid DNA characters found: {sorted(invalid)}")
    if not seq:
        raise ValueError("Sequence is empty.")
    return seq


def gc_content(seq: str) -> float:
    """Return GC content as a percentage."""
    seq = clean_sequence(seq)
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    seq = clean_sequence(seq)
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def transcribe(seq: str) -> str:
    """Transcribe the coding DNA strand to mRNA (T -> U)."""
    return clean_sequence(seq).replace("T", "U")


def translate(seq: str, to_stop: bool = False) -> str:
    """Translate a DNA sequence into a protein sequence (frame 1)."""
    seq = clean_sequence(seq)
    protein = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON_TABLE[seq[i:i + 3]]
        if aa == "*" and to_stop:
            break
        protein.append(aa)
    return "".join(protein)


def find_orfs(seq: str, min_length: int = 30) -> List[Tuple[str, int, int, int, str]]:
    """
    Find ORFs (ATG ... stop) on both strands.

    Returns a list of tuples:
        (strand, frame, start, end, protein)
    Coordinates are 1-based, relative to the strand searched.
    min_length is the minimum ORF length in amino acids.
    """
    seq = clean_sequence(seq)
    results = []
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                if s[i:i + 3] == "ATG":
                    j = i
                    protein = []
                    while j < len(s) - 2:
                        aa = CODON_TABLE[s[j:j + 3]]
                        if aa == "*":
                            if len(protein) >= min_length:
                                results.append(
                                    (strand, frame + 1, i + 1, j + 3, "".join(protein))
                                )
                            break
                        protein.append(aa)
                        j += 3
                    i = j + 3  # continue after this ORF / stop codon
                else:
                    i += 3
    return sorted(results, key=lambda r: len(r[4]), reverse=True)


def melting_temp(primer: str) -> float:
    """
    Estimate primer melting temperature in Celsius.

    Uses the Wallace rule (2*(A+T) + 4*(G+C)) for primers < 14 nt,
    and a basic GC-based formula for longer primers:
        Tm = 64.9 + 41 * (G + C - 16.4) / N
    """
    primer = clean_sequence(primer)
    n = len(primer)
    gc = primer.count("G") + primer.count("C")
    at = primer.count("A") + primer.count("T")
    if n < 14:
        return 2.0 * at + 4.0 * gc
    return 64.9 + 41.0 * (gc - 16.4) / n


def read_fasta(path: str) -> str:
    """Read the first record of a FASTA file and return its sequence."""
    seq_lines = []
    started = False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if started:
                    break
                started = True
                continue
            seq_lines.append(line)
    return "".join(seq_lines)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="DNA sequence analysis toolkit (GC content, ORFs, translation, primer Tm)."
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("-s", "--sequence", help="DNA sequence string")
    src.add_argument("-f", "--fasta", help="Path to a FASTA file (first record used)")
    parser.add_argument("--min-orf", type=int, default=30,
                        help="Minimum ORF length in amino acids (default: 30)")
    args = parser.parse_args(argv)

    try:
        seq = clean_sequence(read_fasta(args.fasta) if args.fasta else args.sequence)
    except (ValueError, OSError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"Length            : {len(seq)} bp")
    print(f"GC content        : {gc_content(seq):.2f}%")
    print(f"Reverse complement: {reverse_complement(seq)}")
    print(f"mRNA              : {transcribe(seq)}")
    print(f"Protein (frame 1) : {translate(seq)}")
    print(f"Melting temp (Tm) : {melting_temp(seq):.1f} C")

    orfs = find_orfs(seq, args.min_orf)
    print(f"\nORFs found (>= {args.min_orf} aa): {len(orfs)}")
    for strand, frame, start, end, protein in orfs[:10]:
        print(f"  strand {strand} frame {frame} {start}-{end} ({len(protein)} aa): {protein}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
