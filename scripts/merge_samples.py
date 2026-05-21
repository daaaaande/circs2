"""
merge_samples.py
----------------
Snakemake script: merges per-sample consensus TSVs into a single
normalised count matrix for downstream differential expression analysis.

Output has two sheets / TSVs:
  all_samples_normalized.tsv  –  wide matrix (coord × sample, values = mean_RPM)
  all_samples_raw_BSJ.tsv     –  same but with mean BSJ read counts
"""

import logging
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

in_files = snakemake.input    # noqa: F821
out_file = snakemake.output.matrix  # noqa: F821

frames = []
for f in in_files:
    df = pd.read_csv(f, sep="\t")
    frames.append(df)

if not frames:
    raise RuntimeError("No input files found for merge_samples.py")

all_data = pd.concat(frames, ignore_index=True)
log.info(f"Total circRNA × sample rows: {len(all_data)}")

# annotation columns (take from any sample that has them)
annot_cols = ["coord", "chrom", "start", "end", "strand", "geneName", "circType"]
annot = (
    all_data[annot_cols]
    .drop_duplicates(subset="coord")
    .set_index("coord")
)

# ── mean_RPM wide matrix ──────────────────────────────────────────────────────
rpm_wide = (
    all_data.pivot_table(
        index="coord", columns="sample", values="mean_RPM", aggfunc="first"
    )
    .fillna(0)
)

# ── mean BSJ wide matrix ──────────────────────────────────────────────────────
all_data["mean_BSJ"] = (all_data["cx2_BSJ"] + all_data["ciri2_BSJ"]) / 2
bsj_wide = (
    all_data.pivot_table(
        index="coord", columns="sample", values="mean_BSJ", aggfunc="first"
    )
    .fillna(0)
)

# merge annotation back
rpm_out  = annot.join(rpm_wide, how="right")
bsj_out  = annot.join(bsj_wide, how="right")

out_path = Path(out_file)
out_path.parent.mkdir(parents=True, exist_ok=True)

rpm_out.to_csv(out_path, sep="\t")
bsj_out.to_csv(out_path.parent / "all_samples_raw_BSJ.tsv", sep="\t")

log.info(f"RPM matrix:  {out_path}  ({rpm_out.shape[0]} circRNAs × {rpm_wide.shape[1]} samples)")
log.info(f"BSJ matrix:  {out_path.parent / 'all_samples_raw_BSJ.tsv'}")
