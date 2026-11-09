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

Zero circRNAs is a valid result: header-only matrices are written, and every
sample still gets a column (all zeros). For that, pass params.samples.
"""

import logging
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

in_files      = list(snakemake.input)                       # noqa: F821
out_dir       = Path(snakemake.output.cx2_raw).parent       # noqa: F821
min_bsj       = snakemake.params.get("min_bsj", 2)          # noqa: F821
samples_param = snakemake.params.get("samples", None)       # noqa: F821

if not in_files:
    raise RuntimeError("No input files found for merge_samples.py")

ANNOT_COLS = ["coord", "chrom", "start", "end", "strand", "geneName", "circType"]
VALUE_COLS = ["cx2_BSJ", "ciri2_BSJ", "cx2_RPM", "ciri2_RPM"]

# ── load all per-sample consensus TSVs ───────────────────────────────────────
# chrom/strand as strings (numeric contigs!); empty (header-only) files are skipped
# for the concat so they cannot upcast the numeric dtypes.
frames = []
for f in in_files:
    df = pd.read_csv(f, sep="\t", dtype={"chrom": str, "strand": str},
                     keep_default_na=False)
    if not df.empty:
        frames.append(df)

if frames:
    all_data = pd.concat(frames, ignore_index=True)
else:
    all_data = pd.DataFrame({c: pd.Series(dtype=object)
                             for c in ANNOT_COLS + ["sample"] + VALUE_COLS})
log.info(f"Loaded {len(all_data)} circRNA × sample rows from {len(in_files)} samples "
         f"({len(in_files) - len(frames)} with 0 consensus circRNAs)")

# ── sample list: every sample gets a column, even with 0 circRNAs ────────────
samples_in_data = list(dict.fromkeys(all_data["sample"]))
if samples_param:
    samples = sorted(samples_param)
    unknown = set(samples_in_data) - set(samples)
    if unknown:
        raise ValueError(f"Samples in data but not in params.samples: {sorted(unknown)}")
else:
    samples = sorted(samples_in_data)
    if len(samples) < len(in_files):
        log.warning(f"{len(in_files) - len(samples)} sample(s) have 0 circRNAs and will be "
                    f"MISSING from the matrices - pass params.samples to keep them")

# ── annotation: one row per unique circRNA coord ──────────────────────────────
annot = (
    all_data[ANNOT_COLS]
    .rename(columns={"geneName": "gene_name", "circType": "circ_type"})
    .drop_duplicates(subset="coord")
    .set_index("coord")
)

# ── build four wide matrices ──────────────────────────────────────────────────
def make_wide(data: pd.DataFrame, value_col: str, as_int: bool = False) -> pd.DataFrame:
    """Pivot to coord × sample, fill missing with 0."""
    if data.empty:
        wide = pd.DataFrame(0 if as_int else 0.0,
                            index=pd.Index([], name="coord"), columns=samples)
    else:
        wide = (
            data.pivot_table(index="coord", columns="sample",
                             values=value_col, aggfunc="first")
            .reindex(columns=samples)
            .fillna(0)
        )
        if as_int:
            wide = wide.astype("int64")
    wide.columns.name = None
    return wide

cx2_raw_wide   = make_wide(all_data, "cx2_BSJ",   as_int=True)
ciri2_raw_wide = make_wide(all_data, "ciri2_BSJ", as_int=True)
cx2_rpm_wide   = make_wide(all_data, "cx2_RPM")
ciri2_rpm_wide = make_wide(all_data, "ciri2_RPM")

# ── global filter: at least min_bsj reads in >= 1 sample in EITHER tool ──────
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
    df.index.name = "coord"
    return df.reset_index()

cx2_raw_out   = annotate(cx2_raw_wide)
ciri2_raw_out = annotate(ciri2_raw_wide)
cx2_rpm_out   = annotate(cx2_rpm_wide)
ciri2_rpm_out = annotate(ciri2_rpm_wide)

# ── write output ──────────────────────────────────────────────────────────────
out_dir.mkdir(parents=True, exist_ok=True)

def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)      # header-only if 0 rows
    log.info(f"Written: {path}  ({df.shape[0]} circRNAs × {df.shape[1]} columns)")

write_tsv(cx2_raw_out,   Path(snakemake.output.cx2_raw))    # noqa: F821
write_tsv(cx2_rpm_out,   Path(snakemake.output.cx2_rpm))    # noqa: F821
write_tsv(ciri2_raw_out, Path(snakemake.output.ciri2_raw))  # noqa: F821
write_tsv(ciri2_rpm_out, Path(snakemake.output.ciri2_rpm))  # noqa: F821