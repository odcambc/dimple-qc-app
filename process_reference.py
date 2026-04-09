from Bio import Align, SeqIO
import pandas as pd


def align_ref_to_variants(
    per_base_df: pd.DataFrame, reference_sequence: str | None
) -> pd.DataFrame:
    if per_base_df.empty:
        return pd.DataFrame()

    if reference_sequence is None:
        return per_base_df

    df_ref_sequence = "".join(per_base_df["ref"])

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0

    alignments = aligner.align(reference_sequence, df_ref_sequence)

    if alignments:
        best_alignment = alignments[0]
        aligned_ref = list(best_alignment[0])    # Reference with gaps inserted
        aligned_query = list(best_alignment[1])  # Data sequence with gaps inserted

        adjusted_aligned_seq = []
        mismatch_flags = []

        for ref_char, query_char in zip(aligned_ref, aligned_query):
            if query_char == "-":
                # Reference has an extra base not present in the data — skip
                continue
            # TODO: implement the per-base annotation logic here.
            # For each data position (query_char != '-'), decide:
            #   - ref_char == '-': data has an insertion not in the reference
            #   - ref_char != query_char: substitution mismatch
            #   - ref_char == query_char: match
            # Append an entry to adjusted_aligned_seq and mismatch_flags (0 or 1).

        per_base_df["aligned_ref"] = adjusted_aligned_seq
        per_base_df["alignment_mismatch"] = mismatch_flags

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
