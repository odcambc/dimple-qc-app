from pathlib import Path

import pandas as pd

app_dir = Path(__file__).parent
# example_df = pd.read_csv(app_dir / "FKYSRV_1_PXR2.tsv")
example_per_base_tsv = app_dir / "FKYSRV_1_PXR2.tsv"
example_reference_fasta = app_dir / "FKYSRV_1_PXR2.fasta"
