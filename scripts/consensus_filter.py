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

# top: replace parse_ciri2() with this
def parse_ciri_full(path: str) -> pd.DataFrame:
    """
    Parse CIRI-full .ciri output (same BSJ columns as CIRI2, produced
    by the internal CIRI2 detection step inside the CIRI-full pipeline).
    Relevant columns: circRNA_ID | chr | circRNA_start | circRNA_end |
                      #junction_reads | strand
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    df.columns = df.columns.str.strip()
    df = df.rename(columns={
        "chr":             "chrom",
        "circRNA_start":   "start",
        "circRNA_end":     "end",
        "#junction_reads": "ciri2_BSJ",
    })
    df["start"]     = df["start"].astype(int) - 1   # 1-based → 0-based
    df["end"]       = df["end"].astype(int)
    df["coord"]     = df["chrom"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str) + ":" + df["strand"]
    df["ciri2_BSJ"] = df["ciri2_BSJ"].astype(int)
    return df[["coord", "ciri2_BSJ"]]




# ── main ──────────────────────────────────────────────────────────────────────

sample     = snakemake.wildcards.sample          # noqa: F821
cx2_path   = snakemake.input.cx2                 # noqa: F821
ciri_full_path = snakemake.input.ciri_full   # was: snakemake.input.ciri2star_log   = snakemake.input.star_log            # noqa: F821
out_path   = snakemake.output.consensus          # noqa: F821

log.info(f"[{sample}] Parsing STAR log: {star_log}")
total_mapped = parse_star_log(star_log)
log.info(f"[{sample}] Uniquely mapped reads: {total_mapped:,}")

log.info(f"[{sample}] Parsing CIRCexplorer2: {cx2_path}")
cx2 = parse_circexplorer2(cx2_path)
log.info(f"[{sample}] CIRCexplorer2 circRNAs: {len(cx2)}")

log.info(f"[{sample}] Parsing CIRI2: {ciri2_path}")
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
