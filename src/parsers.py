"""Parse MSP output files from RFMix and Gnomix into a common format."""

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class AncestryInterval:
    """A contiguous genomic interval assigned to one ancestry on one haplotype."""

    individual: str
    haplotype: int  # 0 or 1
    start: int
    end: int
    ancestry: str


@dataclass
class MSPResult:
    """Parsed MSP output with ancestry code mapping and per-window calls."""

    tool: str
    pop_order: dict[int, str]  # code -> label (e.g. {0: "AFR", 1: "EUR"})
    windows: pd.DataFrame  # columns: chm, spos, epos, + haplotype columns
    hap_columns: list[str]  # e.g. ["tsk_0.0", "tsk_0.1", ...]


def parse_msp(msp_path: Path, tool: str) -> MSPResult:
    """Parse an MSP file from either RFMix or Gnomix.

    Both tools produce the same format:
    - Line 1: #Subpopulation order/codes: AFR=0  EUR=1 ...
    - Line 2: header with chm, spos, epos, sgpos, egpos, n snps, then haplotype columns
    - Remaining lines: data rows
    """
    with open(msp_path) as f:
        header_line = f.readline().strip()
        col_line = f.readline().strip()

    # Parse population order from first line
    # Format: "#Subpopulation order/codes: AFR=0\tEUR=1" or similar
    pop_order = {}
    codes_part = header_line.split(":", 1)[1].strip()
    for pair in re.split(r"\s+", codes_part):
        pair = pair.strip()
        if "=" in pair:
            label, code = pair.split("=")
            pop_order[int(code)] = label

    # Parse column names from second line (strip leading #)
    columns = col_line.lstrip("#").split("\t")

    # Read data, skipping both header lines, using our cleaned column names
    df = pd.read_csv(msp_path, sep="\t", skiprows=[0, 1], header=None, names=columns)

    # Identify haplotype columns (everything after the 6 metadata columns)
    meta_cols = {"chm", "spos", "epos", "sgpos", "egpos", "n snps"}
    hap_cols = [c for c in columns if c not in meta_cols]

    return MSPResult(
        tool=tool,
        pop_order=pop_order,
        windows=df,
        hap_columns=hap_cols,
    )


def msp_to_intervals(result: MSPResult) -> list[AncestryInterval]:
    """Convert MSP window-level calls to merged ancestry intervals.

    Merges consecutive windows with the same ancestry for each haplotype.
    """
    intervals = []
    df = result.windows

    for hap_col in result.hap_columns:
        # Parse individual name and haplotype from column name (e.g. "tsk_0.0")
        parts = hap_col.rsplit(".", 1)
        indiv = parts[0]
        hap = int(parts[1])

        # Walk through windows, merging consecutive same-ancestry spans
        cur_anc = None
        cur_start = None
        cur_end = None

        for _, row in df.iterrows():
            anc_code = int(row[hap_col])
            anc_label = result.pop_order[anc_code]
            spos = int(row["spos"])
            epos = int(row["epos"])

            if anc_label == cur_anc:
                cur_end = epos
            else:
                if cur_anc is not None:
                    intervals.append(
                        AncestryInterval(indiv, hap, cur_start, cur_end, cur_anc)
                    )
                cur_anc = anc_label
                cur_start = spos
                cur_end = epos

        if cur_anc is not None:
            intervals.append(AncestryInterval(indiv, hap, cur_start, cur_end, cur_anc))

    return intervals


def intervals_to_df(intervals: list[AncestryInterval]) -> pd.DataFrame:
    """Convert a list of AncestryIntervals to a DataFrame."""
    return pd.DataFrame(
        [
            {
                "individual": iv.individual,
                "haplotype": iv.haplotype,
                "start": iv.start,
                "end": iv.end,
                "ancestry": iv.ancestry,
            }
            for iv in intervals
        ]
    )


def get_rfmix_msp_path(scenario: str, generation: int | None, chrom: int) -> Path:
    """Get the MSP output path for an RFMix run."""
    import config

    if generation is not None:
        return config.RFMIX_DIR / scenario / f"gen_{generation}" / f"chr{chrom}.msp.tsv"
    return config.RFMIX_DIR / scenario / "real" / f"chr{chrom}.msp.tsv"


def get_gnomix_msp_path(scenario: str, generation: int | None) -> Path:
    """Get the MSP output path for a Gnomix run."""
    import config

    if generation is not None:
        return config.GNOMIX_DIR / scenario / f"gen_{generation}" / "query_results.msp"
    return config.GNOMIX_DIR / scenario / "real" / "query_results.msp"
