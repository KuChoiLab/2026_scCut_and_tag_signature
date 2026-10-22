#!/usr/bin/env python3
"""Collect per-sample QC (results/<sample>/qc/) into one table and flag samples against ENCODE standards.

Usage:
  python3 qc_aggregate.py <samples.txt> [--results DIR] [--out FILE]

  samples.txt : the sample sheet given to chip_run.sh (defines which samples, in which order)
  --results   : PROJECT_RESULTS (default: the path in config.sh)
  --out       : default <results>/QC_summary_<samples name>.tsv

One row per sample (single replicate per sample: chip_rep1 vs ctl_rep1, peaks of rep1).
Samples whose finalize step has not produced qc/ yet are listed with status=not_finished.

Thresholds (ENCODE histone ChIP-seq standards / ENCODE audit categories); edit below if needed:
  usable fragments (ChIP) : >= 20M for narrow marks such as H3K4me3 (low < 20M, very_low < 10M)
  NRF                     : < 0.5 concerning, < 0.8 low (0.8-0.9 compliant, > 0.9 ideal)
  PBC1                    : < 0.5 severe, < 0.8 moderate, < 0.9 mild bottlenecking
  PBC2                    : < 1 severe, < 3 moderate, < 10 mild bottlenecking
  NSC / RSC               : NSC < 1.05 or RSC < 0.8 -> weak enrichment (cross-correlation)
  FRiP (rep1 peaks)       : < 0.01 -> low (no histone-specific ENCODE cut-off; informational)
  control                 : fewer usable fragments than the ChIP -> ctl_shallower (informational)
"""
import argparse, csv, os, re, sys

USABLE_LOW, USABLE_VERY_LOW = 20e6, 10e6
NRF_CONCERN, NRF_LOW = 0.5, 0.8
PBC1_CUTS = [(0.5, "severe"), (0.8, "moderate"), (0.9, "mild")]
PBC2_CUTS = [(1, "severe"), (3, "moderate"), (10, "mild")]
NSC_MIN, RSC_MIN = 1.05, 0.8
FRIP_LOW = 0.01

PIPE = os.path.dirname(os.path.abspath(__file__))


def default_results():
    m = re.search(r"^PROJECT_RESULTS=\$\{PROJECT_RESULTS:-([^}]+)\}", open(os.path.join(PIPE, "config.sh")).read(), re.M)
    return os.path.expandvars(m.group(1))


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def bottleneck(pbc1, pbc2):
    rank = {"none": 0, "mild": 1, "moderate": 2, "severe": 3}
    worst = "none"
    for v, cuts in ((pbc1, PBC1_CUTS), (pbc2, PBC2_CUTS)):
        if v is None:
            continue
        for cut, lab in cuts:
            if v < cut:
                worst = max(worst, lab, key=rank.get)
                break
    return worst


def fmt(v, nd=3):
    if v is None:
        return "NA"
    if isinstance(v, float) and v.is_integer() and abs(v) >= 1000:
        return str(int(v))
    return f"{v:.{nd}f}" if isinstance(v, float) else str(v)


def sample_row(sid, sdir):
    parts = sid.split("_")
    row = {"sample": sid, "cell_type": parts[1] if len(parts) > 2 else "NA",
           "donor": parts[-1] if len(parts) > 2 else "NA"}
    q = os.path.join(sdir, "qc")
    if not (os.path.exists(os.path.join(q, "align_qc.tsv")) and os.path.exists(os.path.join(q, "peak_qc.tsv"))):
        row["status"] = "not_finished"
        return row
    row["status"] = "done"

    units = {r[0] + "_" + r[1]: r[2] == "1" for r in csv.reader(open(os.path.join(sdir, "units.tsv")), delimiter="\t")}
    align = {r["unit"]: r for r in read_tsv(os.path.join(q, "align_qc.tsv"))}
    peaks = {r["peak_set"]: r for r in read_tsv(os.path.join(q, "peak_qc.tsv"))}
    flags = []

    for role in ("chip", "ctl"):
        a = align.get(f"{role}_rep1", {})
        pe = units.get(f"{role}_rep1", False)
        reads = num(a.get("filt_nodup_reads"))
        usable = reads / 2 if (reads is not None and pe) else reads   # PE: two reads per fragment
        row[f"{role}_paired"] = "PE" if pe else "SE"
        row[f"{role}_total_reads"] = num(a.get("total_reads"))
        row[f"{role}_mapped_pct"] = num(a.get("mapped_pct"))
        row[f"{role}_usable_frags"] = usable
        row[f"{role}_dup_pct"] = num(a.get("dup_pct"))
        for k in ("NRF", "PBC1", "PBC2"):
            row[f"{role}_{k}"] = num(a.get(k))
        row[f"{role}_bottleneck"] = bottleneck(row[f"{role}_PBC1"], row[f"{role}_PBC2"])
        nrf = row[f"{role}_NRF"]
        if nrf is not None and nrf < NRF_LOW:
            flags.append(f"{role}_NRF_{'concerning' if nrf < NRF_CONCERN else 'low'}")
        if row[f"{role}_bottleneck"] in ("moderate", "severe"):
            flags.append(f"{role}_{row[f'{role}_bottleneck']}_bottleneck")

    chip = align.get("chip_rep1", {})
    row["fraglen"], row["NSC"], row["RSC"] = num(chip.get("fraglen")), num(chip.get("NSC")), num(chip.get("RSC"))
    u = row["chip_usable_frags"]
    if u is not None and u < USABLE_LOW:
        flags.append("chip_very_low_depth" if u < USABLE_VERY_LOW else "chip_low_depth")
    if row["NSC"] is not None and row["NSC"] < NSC_MIN:
        flags.append("low_NSC")
    if row["RSC"] is not None and row["RSC"] < RSC_MIN:
        flags.append("low_RSC")
    c = row["ctl_usable_frags"]
    if u is not None and c is not None and c < u:
        flags.append("ctl_shallower")

    p = peaks.get("rep1", {})
    row["num_peaks"], row["FRiP"] = num(p.get("num_peaks")), num(p.get("FRiP"))
    row["peak_median_bp"], row["peak_mean_bp"] = num(p.get("size_median")), num(p.get("size_mean"))
    row["pr_overlap_peaks"] = num(peaks.get("overlap:rep1-pr1_vs_rep1-pr2", {}).get("num_peaks"))
    row["optimal_peaks"] = num(peaks.get("final:optimal", {}).get("num_peaks"))
    if row["FRiP"] is not None and row["FRiP"] < FRIP_LOW:
        flags.append("low_FRiP")

    row["flags"] = ";".join(flags) if flags else "PASS"
    return row


COLS = ["sample", "cell_type", "donor", "status", "flags",
        "chip_paired", "chip_total_reads", "chip_mapped_pct", "chip_usable_frags", "chip_dup_pct",
        "chip_NRF", "chip_PBC1", "chip_PBC2", "chip_bottleneck", "fraglen", "NSC", "RSC",
        "ctl_paired", "ctl_total_reads", "ctl_mapped_pct", "ctl_usable_frags", "ctl_dup_pct",
        "ctl_NRF", "ctl_PBC1", "ctl_PBC2", "ctl_bottleneck",
        "num_peaks", "FRiP", "peak_median_bp", "peak_mean_bp", "pr_overlap_peaks", "optimal_peaks"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples")
    ap.add_argument("--results", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    res = a.results or default_results()
    ids = []
    for line in open(a.samples):
        if line.startswith("#") or not line.strip():
            continue
        sid = line.split("\t")[0]
        if sid not in ids:
            ids.append(sid)
    rows = [sample_row(s, os.path.join(res, s)) for s in ids]
    name = os.path.splitext(os.path.basename(a.samples))[0].replace("samples_", "")
    out = a.out or os.path.join(res, f"QC_summary_{name}.tsv")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c)) if not isinstance(r.get(c), str) else r[c] for c in COLS) + "\n")

    # short view on screen
    show = [("sample", 22), ("status", 12), ("chip_usable_frags", 10), ("chip_NRF", 6), ("chip_PBC1", 6),
            ("chip_PBC2", 7), ("NSC", 6), ("RSC", 6), ("FRiP", 6), ("num_peaks", 9), ("flags", 0)]
    hdr = {"chip_usable_frags": "usable(M)", "chip_NRF": "NRF", "chip_PBC1": "PBC1", "chip_PBC2": "PBC2", "num_peaks": "peaks"}
    print("  ".join(f"{hdr.get(c, c):<{w}}" if w else hdr.get(c, c) for c, w in show))
    for r in rows:
        cells = []
        for c, w in show:
            v = r.get(c)
            if c == "chip_usable_frags" and isinstance(v, float):
                v = f"{v / 1e6:.1f}"
            elif isinstance(v, float):
                v = fmt(v, 0 if c == "num_peaks" else 3 if c not in ("chip_PBC2",) else 1)
            v = "" if v is None else str(v)
            cells.append(f"{v:<{w}}" if w else v)
        print("  ".join(cells))
    print(f"wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
