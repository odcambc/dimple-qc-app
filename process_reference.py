from Bio import Align, SeqIO
import numpy as np
import pandas as pd

from process_data import _extract_ref_base_from_aligned, compute_max_non_ref_base


def align_ref_to_variants(
    per_base_df: pd.DataFrame, reference_sequence: str | None
) -> pd.DataFrame:
    if per_base_df.empty:
        return pd.DataFrame()

    if reference_sequence is None:
        return per_base_df

    df_ref_sequence = "".join(per_base_df["ref"])

    try:
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        alignments = aligner.align(reference_sequence, df_ref_sequence)

        try:
            best_alignment = alignments[0]
        except (IndexError, OverflowError):
            return per_base_df

        ref_aligned = str(best_alignment[0])
        df_aligned = str(best_alignment[1])

        # Build aligned_ref and alignment_mismatch by iterating the full alignment.
        # Track df position separately to handle gaps without zip truncation.
        adjusted_aligned_seq = []
        mismatch_flags = []
        df_pos = 0

        for ref_base, df_base in zip(ref_aligned, df_aligned):
            if df_base == "-":
                # Insertion in reference (deletion in data) — skip this position
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

        # Only write back if lengths match; otherwise leave placeholder values
        if len(adjusted_aligned_seq) == len(per_base_df):
            per_base_df["aligned_ref"] = adjusted_aligned_seq
            per_base_df["alignment_mismatch"] = mismatch_flags

            # Recompute max_variant_base from the aligned reference
            bases = per_base_df[["A", "C", "G", "T"]].to_numpy()
            ref_bases = np.array([
                _extract_ref_base_from_aligned(v) or ""
                for v in adjusted_aligned_seq
            ])
            per_base_df["max_variant_base"] = compute_max_non_ref_base(
                per_base_df, bases=bases, ref_bases=ref_bases
            )
    except Exception:
        pass

    return per_base_df


def process_reference_fasta(file) -> dict[str, dict | None] | None:

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


def process_reference_genbank(file) -> dict[str, dict | None] | None:
    try:
        genbank_record = list(SeqIO.parse(file[0]["datapath"], "genbank"))
    except Exception:
        return None

    # Only parse single-record files
    if len(genbank_record) == 1:
        genbank_record = list(genbank_record)[0]
        sequence = str(genbank_record.seq)

        features = genbank_record.features
    else:
        return None

    features_dict = {}
    for feature in features:
        try:
            if "label" in feature.qualifiers:
                feature_name = feature.qualifiers["label"][0]
            else:
                feature_name = feature.type

            features_dict[feature_name] = feature
        except AttributeError:
            continue

    return {"sequence": sequence, "features": features_dict}
