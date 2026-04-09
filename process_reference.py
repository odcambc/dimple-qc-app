from Bio import SeqIO, Align
import pandas as pd


def align_ref_to_variants(
    per_base_df: pd.DataFrame, reference_sequence: str | None
) -> pd.DataFrame:
    if per_base_df.empty:
        return pd.DataFrame()

    if reference_sequence is None:
        return per_base_df

    df_ref_sequence = "".join(per_base_df["ref"])

    # Use Bio.Align.PairwiseAligner (replaces deprecated pairwise2)
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0

    alignments = aligner.align(reference_sequence, df_ref_sequence)

    if alignments:
        best_alignment = alignments[0]
        ref_aligned = str(best_alignment.seqA)
        df_aligned = str(best_alignment.seqB)

        # Build aligned_ref and alignment_mismatch by iterating through full alignment
        # Track position in df separately to handle gaps correctly
        adjusted_aligned_seq = []
        mismatch_flags = []
        df_pos = 0

        for ref_base, df_base in zip(ref_aligned, df_aligned):
            if df_base == "-":
                # Insertion in reference (deletion in data) — skip, don't add to aligned_ref
                continue
            elif ref_base == "-":
                # Deletion in reference (insertion in data) — mark as gap
                adjusted_aligned_seq.append("-")
                mismatch_flags.append(0)
                df_pos += 1
            elif ref_base != df_base:
                # Mismatch — mark with brackets
                adjusted_aligned_seq.append(f"[{ref_base}]")
                mismatch_flags.append(1)
                df_pos += 1
            else:
                # Match
                adjusted_aligned_seq.append(ref_base)
                mismatch_flags.append(0)
                df_pos += 1

        # Ensure we have the right number of aligned positions
        if len(adjusted_aligned_seq) == len(per_base_df):
            per_base_df["aligned_ref"] = adjusted_aligned_seq
            per_base_df["alignment_mismatch"] = mismatch_flags
        # If lengths don't match, leave the placeholder values in place

    return per_base_df


def process_reference_fasta(file) -> dict[str, str | list | None] | None:

    try:
        fasta_record = list(SeqIO.parse(file[0]["datapath"], "fasta"))
    except Exception:
        return None

    if len(fasta_record) == 1:
        sequence = str(fasta_record[0].seq)
    else:
        return None

    # There are no features to parse from a fasta file
    return {"sequence": sequence, "features": None}


def process_reference_genbank(file) -> dict[str, str | list | None] | None:
    try:
        genbank_record = list(SeqIO.parse(file[0]["datapath"], "genbank"))
    except Exception:
        return None

    # Only parse single-record files
    if len(genbank_record) == 1:
        genbank_record = list(genbank_record)[0]
        sequence = str(genbank_record.seq)

        # TODO: does features always return a list? Check
        features = genbank_record.features
    else:
        return None

    return {"sequence": sequence, "features": features}
