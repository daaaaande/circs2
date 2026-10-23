"""
consensus_filter.py
-------------------
Snakemake script: takes per-sample CIRCexplorer2 and CIRI2 output,
applies a 2-of-2 consensus filter, and normalises to RPM.

RPM = (BSJ_reads / uniquely_mapped_reads) × 1_000_000

The normalisation denominator is the number of uniquely mapped reads from the
STAR Log.final.out (same alignment used by CIRCexplorer2, and a reliable proxy
for library size).  CIRI2 BSJ counts are kept as a separate column.

An empty result (0 consensus circRNAs) is a valid output: a header-only TSV
is written.
"""

import re
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

OUT_COLS = [
    "sample", "coord", "chrom", "start", "end", "strand",
    "geneName", "circType",
    "cx2_BSJ", "ciri2_BSJ",
    "total_mapped",
    "cx2_RPM", "ciri2_RPM", "mean_RPM",
]

CX2_COLS = ["coord", "chrom", "start", "end", "strand", "geneName", "circType", "cx2_BSJ"]
CIRI_COLS = ["coord", "ciri2_BSJ"]


# ── helpers ───────────────────────────────────────────────────────────────────

def empty_cx2() -> pd.DataFrame:
    return pd.DataFrame({
        "coord": pd.Series(dtype="object"),
        "chrom": pd.Series(dtype="object"),
        "start": pd.Series(dtype="int64"),
        "end": pd.Series(dtype="int64"),
        "strand": pd.Series(dtype="object"),
        "geneName": pd.Series(dtype="object"),
        "circType": pd.Series(dtype="object"),
        "cx2_BSJ": pd.Series(dtype="int64"),
    })


def empty_ciri() -> pd.DataFrame:
    return pd.DataFrame({
        "coord": pd.Series(dtype="object"),
        "ciri2_BSJ": pd.Series(dtype="int64"),
    })


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
    Parse CIRCexplorer2 annotate output (headerless BED12+).
    Columns (0-based):
      0  chrom | 1  start | 2  end | 3  name | 4  score | 5  strand |
      6  thickStart | 7  thickEnd | 8  itemRgb | 9  exonCount |
      10 exonSizes | 11 exonOffsets | 12 readNumber | 13 circType |
      14 geneName | 15 isoformName | 16 index | 17 flankIntron
    Returns an empty (typed) DataFrame if the file is empty / has no records.
    """
    col_names = [
        "chrom", "start", "end", "name", "score", "strand",
        "thickStart", "thickEnd", "itemRgb", "exonCount",
        "exonSizes", "exonOffsets", "readNumber", "circType",
        "geneName", "isoformName", "index", "flankIntron",
    ]
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=col_names,
            comment="#",
            dtype=str,              # everything as str; cast numerics explicitly
            keep_default_na=False,  # don't turn contig "NA"/"N" etc. into NaN
        )
    except pd.errors.EmptyDataError:
        log.warning(f"CIRCexplorer2 file is empty: {path}")
        return empty_cx2()

    # keep only annotated circRNAs (ciRNAs and exonic circRNAs)
    df = df[df["circType"].isin(["circRNA", "ciRNA"])].copy()
    if df.empty:
        return empty_cx2()

    df["start"] = df["start"].astype("int64")
    df["end"] = df["end"].astype("int64")
    df["cx2_BSJ"] = df["readNumber"].astype("int64")
    df["coord"] = (df["chrom"] + ":" + df["start"].astype(str) + "-" +
                   df["end"].astype(str) + ":" + df["strand"])
    return df[CX2_COLS]


def parse_ciri_full(path: str) -> pd.DataFrame:
    """
    Parse CIRI2 / CIRI-full TSV output (line-based, tolerant of ragged rows).
    Columns (0-based):
      0  circRNA_ID | 1  chr | 2  circRNA_start | 3  circRNA_end |
      4  #junction_reads | 5  SM_MS_SMS | 6  #non_junction_reads |
      7  junction_reads_ratio | 8  circRNA_type | 9  gene_id |
      10 strand | 11 junction_reads_ID
    Coordinates are 1-based; start is converted to 0-based here.
    Returns an empty (typed) DataFrame if the file has no data rows.
    """
    rows = []
    n_bad = 0
    with open(path) as fh:
        header = fh.readline()          # skip header (may be absent if file is empty)
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            f = line.split("\t")
            try:
                chrom, start, end = f[1], int(f[2]) - 1, int(f[3])
                bsj, strand = int(f[4]), f[10]
            except (IndexError, ValueError):
                n_bad += 1
                continue
            rows.append((f"{chrom}:{start}-{end}:{strand}", bsj))

    if n_bad:
        log.warning(f"CIRI2: skipped {n_bad} malformed line(s) in {path}")
    if not rows:
        log.warning(f"CIRI2 file has no data rows: {path}")
        return empty_ciri()

    return pd.DataFrame(rows, columns=CIRI_COLS)


# ── main ──────────────────────────────────────────────────────────────────────

sample         = snakemake.wildcards.sample      # noqa: F821
cx2_path       = snakemake.input.cx2             # noqa: F821
ciri_full_path = snakemake.input.ciri2           # noqa: F821  (Snakefile input key is still 'ciri2')
star_log       = snakemake.input.star_log        # noqa: F821
out_path       = snakemake.output.consensus      # noqa: F821

log.info(f"[{sample}] Parsing STAR log: {star_log}")
total_mapped = parse_star_log(star_log)
log.info(f"[{sample}] Uniquely mapped reads: {total_mapped:,}")
if total_mapped == 0:
    raise ValueError(f"[{sample}] 0 uniquely mapped reads - cannot normalise to RPM")

log.info(f"[{sample}] Parsing CIRCexplorer2: {cx2_path}")
cx2 = parse_circexplorer2(cx2_path)
log.info(f"[{sample}] CIRCexplorer2 circRNAs: {len(cx2)}")

log.info(f"[{sample}] Parsing CIRI2: {ciri_full_path}")
ciri2 = parse_ciri_full(ciri_full_path)
log.info(f"[{sample}] CIRI2 circRNAs: {len(ciri2)}")

# 2-of-2 consensus: inner join on coordinate key (works fine if either side is empty)
consensus = pd.merge(cx2, ciri2, on="coord", how="inner")
log.info(f"[{sample}] Consensus (2-of-2) circRNAs: {len(consensus)}")

# RPM normalisation
consensus["cx2_RPM"]   = (consensus["cx2_BSJ"]  / total_mapped) * 1e6
consensus["ciri2_RPM"] = (consensus["ciri2_BSJ"] / total_mapped) * 1e6
consensus["mean_RPM"]  = (consensus["cx2_RPM"] + consensus["ciri2_RPM"]) / 2

consensus["sample"]       = sample
consensus["total_mapped"] = total_mapped

consensus = consensus[OUT_COLS].sort_values(["chrom", "start", "end"])
consensus.to_csv(out_path, sep="\t", index=False)   # header-only file if 0 rows
log.info(f"[{sample}] Written to {out_path} ({len(consensus)} rows)")