from Bio import pairwise2, SeqIO
import pandas as pd


def align_ref_to_variants(
    per_base_df: pd.DataFrame, reference_sequence: str | None
) -> pd.DataFrame:
    if per_base_df.empty:
        return pd.DataFrame()

    if reference_sequence is None:
        return per_base_df

    df_ref_sequence = "".join(per_base_df["ref"])

    alignments = pairwise2.align.globalxx(
        reference_sequence, df_ref_sequence, one_alignment_only=True
    )

    if alignments:
        best_alignment = alignments[0]  # Best alignment result
        aligned_seq = list(best_alignment.seqA)  # Extract aligned reference

        # Handle insertions/deletions (gap handling)
        adjusted_aligned_seq = []
        mismatch_count = 0

        for i, (ref_base, df_base) in enumerate(zip(aligned_seq, df_ref_sequence)):
            if ref_base == "-":  # Reference has an insertion
                adjusted_aligned_seq.append("-")  # Keep the gap
            elif ref_base != df_base:
                adjusted_aligned_seq.append(f"[{ref_base}]")  # Mark mismatch
                mismatch_count += 1
            else:
                adjusted_aligned_seq.append(ref_base)

        per_base_df["aligned_ref"] = adjusted_aligned_seq

        per_base_df["alignment_mismatch"] = [
            1 if ref_base != df_base else 0
            for ref_base, df_base in zip(aligned_seq, df_ref_sequence)
        ]

    return per_base_df


def process_reference_fasta(file):
    try:
        fasta_record = list(SeqIO.parse(file[0]["datapath"], "fasta"))
        if len(fasta_record) != 1:
            return None  # Only handle single-sequence FASTA files
        return str(fasta_record[0].seq)
    except Exception:
        return None


def process_reference_genbank(file):
    try:
        genbank_record = list(SeqIO.parse(file[0]["datapath"], "genbank"))
        if len(genbank_record) != 1:
            return None  # Only handle single-sequence GenBank files
        return str(genbank_record[0].seq)
    except Exception:
        return None
