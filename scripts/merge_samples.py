"""
merge_samples.py
----------------
Snakemake script: merges per-sample consensus TSVs into four wide matrices.

All four matrices contain exactly the same circRNAs (consensus coordinates,
detected by both CIRCexplorer2 AND CIRI2 in at least one sample with >= min_bsj reads).

Output files (in results/final/):
  cx2_raw.tsv        CIRCexplorer2 BSJ read counts  (raw)
  cx2_rpm.tsv        CIRCexplorer2 RPM              (normalised)
  ciri2_raw.tsv      CIRI2 BSJ read counts          (raw)
  ciri2_rpm.tsv      CIRI2 RPM                      (normalised)

Annotation columns present in all four files:
  coord      chr:start-end:strand  (0-based BED coordinates)
  chrom      chromosome
  start      0-based start
  end        end
  strand     + / -
  gene_name  from CIRCexplorer2 annotation (which uses the pipeline GTF)
  circ_type  circRNA / ciRNA
"""

import logging
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

in_files  = snakemake.input    # noqa: F821
out_dir   = Path(snakemake.output.cx2_raw).parent  # noqa: F821
min_bsj   = snakemake.params.get("min_bsj", 2)     # noqa: F821

# ── load all per-sample consensus TSVs ───────────────────────────────────────
frames = []
for f in in_files:
    df = pd.read_csv(f, sep="\t")
    frames.append(df)

if not frames:
    raise RuntimeError("No input files found for merge_samples.py")

all_data = pd.concat(frames, ignore_index=True)
log.info(f"Loaded {len(all_data)} circRNA × sample rows from {len(frames)} samples")

# ── annotation: one row per unique circRNA coord ──────────────────────────────
# CIRCexplorer2 provides gene_name and circ_type; take first occurrence per coord
annot = (
    all_data[["coord", "chrom", "start", "end", "strand", "geneName", "circType"]]
    .rename(columns={"geneName": "gene_name", "circType": "circ_type"})
    .drop_duplicates(subset="coord")
    .set_index("coord")
)

# ── build four wide matrices ──────────────────────────────────────────────────
def make_wide(data: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """Pivot to coord × sample, fill missing with 0 (integer-safe)."""
    wide = (
        data.pivot_table(
            index="coord",
            columns="sample",
            values=value_col,
            aggfunc="first"
        )
        .fillna(0)
    )
    wide.columns.name = None   # drop the "sample" axis label
    return wide

cx2_raw_wide   = make_wide(all_data, "cx2_BSJ")
ciri2_raw_wide = make_wide(all_data, "ciri2_BSJ")
cx2_rpm_wide   = make_wide(all_data, "cx2_RPM")
ciri2_rpm_wide = make_wide(all_data, "ciri2_RPM")

# ── global filter: at least min_bsj reads in >= 1 sample in EITHER tool ──────
# (per-tool min was already applied upstream, but after merging across samples
#  a circRNA with 1 read in every sample would still slip through without this)
keep_mask = (
    (cx2_raw_wide.max(axis=1)   >= min_bsj) |
    (ciri2_raw_wide.max(axis=1) >= min_bsj)
)
n_before = len(cx2_raw_wide)
cx2_raw_wide   = cx2_raw_wide[keep_mask]
ciri2_raw_wide = ciri2_raw_wide[keep_mask]
cx2_rpm_wide   = cx2_rpm_wide[keep_mask]
ciri2_rpm_wide = ciri2_rpm_wide[keep_mask]
log.info(f"Global BSJ filter (>= {min_bsj} in >= 1 sample): {n_before} → {len(cx2_raw_wide)} circRNAs")

# ── join annotation and reset coord as a plain column ────────────────────────
def annotate(wide: pd.DataFrame) -> pd.DataFrame:
    df = annot.join(wide, how="right")
    df = df.reset_index().rename(columns={"index": "coord"})
    # ensure counts are integers (RPM stays float)
    return df

cx2_raw_out   = annotate(cx2_raw_wide)
ciri2_raw_out = annotate(ciri2_raw_wide)
cx2_rpm_out   = annotate(cx2_rpm_wide)
ciri2_rpm_out = annotate(ciri2_rpm_wide)

# ── write output ──────────────────────────────────────────────────────────────
out_dir.mkdir(parents=True, exist_ok=True)

def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)
    log.info(f"Written: {path}  ({df.shape[0]} circRNAs × {df.shape[1]} columns)")

write_tsv(cx2_raw_out,   Path(snakemake.output.cx2_raw))    # noqa: F821
write_tsv(cx2_rpm_out,   Path(snakemake.output.cx2_rpm))    # noqa: F821
write_tsv(ciri2_raw_out, Path(snakemake.output.ciri2_raw))  # noqa: F821
write_tsv(ciri2_rpm_out, Path(snakemake.output.ciri2_rpm))  # noqa: F821