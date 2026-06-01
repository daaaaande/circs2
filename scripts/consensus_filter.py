"""
consensus_filter.py
-------------------
Snakemake script: takes per-sample CIRCexplorer2 and CIRI2 output,
applies a 2-of-2 consensus filter, and normalises to RPM.

RPM = (BSJ_reads / uniquely_mapped_reads) × 1_000_000

The normalisation denominator is the number of uniquely mapped reads from the
STAR Log.final.out (same alignment used by CIRCexplorer2, and a reliable proxy
for library size).  CIRI2 BSJ counts are kept as a separate column.
"""

import re
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


# ── helpers ───────────────────────────────────────────────────────────────────

def parse_star_log(path: str) -> int:
    """Return number of uniquely mapped reads from STAR Log.final.out."""
    pattern = re.compile(r"Uniquely mapped reads number \|\s+([\d]+)")
    with open(path) as fh:
        for line in fh:
            m = pattern.search(line)
            if m:
                return int(m.group(1))
    raise ValueError(f"Could not find 'Uniquely mapped reads number' in {path}")


def parse_circexplorer2(path: str) -> pd.DataFrame:
    """
    Parse CIRCexplorer2 annotate output.
    Columns (0-based):
      0  chrom | 1  start | 2  end | 3  name | 4  score | 5  strand |
      6  thickStart | 7  thickEnd | 8  itemRgb | 9  exonCount |
      10 exonSizes | 11 exonOffsets | 12 readNumber | 13 circType |
      14 geneName | 15 isoformName | 16 index | 17 flankIntron
    """
    col_names = [
        "chrom", "start", "end", "name", "score", "strand",
        "thickStart", "thickEnd", "itemRgb", "exonCount",
        "exonSizes", "exonOffsets", "readNumber", "circType",
        "geneName", "isoformName", "index", "flankIntron"
    ]
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=col_names,
        comment="#"
    )
    # keep only annotated circRNAs (ciRNAs and exonic circRNAs)
    df = df[df["circType"].isin(["circRNA", "ciRNA"])].copy()
    df["cx2_BSJ"] = df["readNumber"].astype(int)
    df["coord"] = df["chrom"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str) + ":" + df["strand"]
    return df[["coord", "chrom", "start", "end", "strand", "geneName", "circType", "cx2_BSJ"]]

def parse_ciri_full(path: str) -> pd.DataFrame:
    """
    Parse CIRI2 / CIRI-full TSV output.
    Columns (0-based):
      0  circRNA_ID | 1  chr | 2  circRNA_start | 3  circRNA_end |
      4  #junction_reads | 5  SM_MS_SMS | 6  #non_junction_reads |
      7  junction_reads_ratio | 8  circRNA_type | 9  gene_id |
      10 strand | 11 junction_reads_ID (may contain commas/tabs → ragged)
    Coordinates are 1-based; converted to 0-based here.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        usecols=[1, 2, 3, 4, 10],          # select only what we need by position
        header=0,
        names=["chrom", "start", "end", "ciri2_BSJ", "strand"],
        skiprows=1,                          # skip the real header row (we're supplying names)
        engine="python",                     # more tolerant of ragged rows
        on_bad_lines="warn",                 # skip/warn on malformed lines instead of crashing
    )

    log.info(f"CIRI2 raw shape after read: {df.shape}")

    # drop any rows where key fields didn't parse cleanly
    df = df.dropna(subset=["start", "end", "ciri2_BSJ", "strand"])

    df["start"]     = pd.to_numeric(df["start"],     errors="coerce").astype("Int64") - 1
    df["end"]       = pd.to_numeric(df["end"],        errors="coerce").astype("Int64")
    df["ciri2_BSJ"] = pd.to_numeric(df["ciri2_BSJ"], errors="coerce").astype("Int64")

    df = df.dropna(subset=["start", "end", "ciri2_BSJ"])  # drop any coerced NaNs

    df["coord"] = (df["chrom"] + ":" + df["start"].astype(str) + "-" +
                   df["end"].astype(str) + ":" + df["strand"])

    return df[["coord", "ciri2_BSJ"]]

# ── main ──────────────────────────────────────────────────────────────────────

sample         = snakemake.wildcards.sample      # noqa: F821
cx2_path       = snakemake.input.cx2             # noqa: F821
ciri_full_path = snakemake.input.ciri2           # noqa: F821  (Snakefile input key is still 'ciri2')
star_log       = snakemake.input.star_log        # noqa: F821
out_path       = snakemake.output.consensus      # noqa: F821





#
log.info(f"[{sample}] Parsing STAR log: {star_log}")
total_mapped = parse_star_log(star_log)
log.info(f"[{sample}] Uniquely mapped reads: {total_mapped:,}")

log.info(f"[{sample}] Parsing CIRCexplorer2: {cx2_path}")
cx2 = parse_circexplorer2(cx2_path)
log.info(f"[{sample}] CIRCexplorer2 circRNAs: {len(cx2)}")

# 
log.info(f"[{sample}] Parsing CIRI-full: {ciri_full_path}")
ciri2 = parse_ciri_full(ciri_full_path)      # was: parse_ciri2(ciri2_path)
log.info(f"[{sample}] CIRI2 circRNAs: {len(ciri2)}")

# 2-of-2 consensus: inner join on coordinate key
consensus = pd.merge(cx2, ciri2, on="coord", how="inner")
log.info(f"[{sample}] Consensus (2-of-2) circRNAs: {len(consensus)}")

# RPM normalisation (using CIRCexplorer2 BSJ count as primary; CIRI2 kept for reference)
consensus["cx2_RPM"]   = (consensus["cx2_BSJ"]   / total_mapped) * 1e6
consensus["ciri2_RPM"] = (consensus["ciri2_BSJ"]  / total_mapped) * 1e6
consensus["mean_RPM"]  = (consensus["cx2_RPM"] + consensus["ciri2_RPM"]) / 2

consensus["sample"]          = sample
consensus["total_mapped"]    = total_mapped

# reorder columns
out_cols = [
    "sample", "coord", "chrom", "start", "end", "strand",
    "geneName", "circType",
    "cx2_BSJ", "ciri2_BSJ",
    "total_mapped",
    "cx2_RPM", "ciri2_RPM", "mean_RPM"
]
consensus = consensus[out_cols].sort_values(["chrom", "start", "end"])
consensus.to_csv(out_path, sep="\t", index=False)
log.info(f"[{sample}] Written to {out_path}")
