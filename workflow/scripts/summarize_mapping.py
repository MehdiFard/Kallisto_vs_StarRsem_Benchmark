"""Per-sample mapping QC (STAR, RSEM, kallisto) and compute cost (Snakemake benchmark files)."""
import json
import os
import sys

import pandas as pd

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

samples = snakemake.params.samples
runs = snakemake.params.runs

STAR_FIELDS = {
    "Number of input reads": "star_input_reads",
    "Uniquely mapped reads %": "star_unique_pct",
    "% of reads mapped to multiple loci": "star_multi_pct",
    "% of reads unmapped: too short": "star_unmapped_short_pct",
}


def parse_star(path):
    out = {}
    with open(path) as f:
        for line in f:
            if "|" not in line:
                continue
            key, value = (x.strip() for x in line.split("|", 1))
            if key in STAR_FIELDS:
                out[STAR_FIELDS[key]] = float(value.rstrip("%"))
    return out


def parse_rsem_cnt(path):
    # line 1: N0 (unalignable) N1 (alignable) N2 (filtered) N_tot; line 2: nUnique nMulti nUncertain
    with open(path) as f:
        n0, n1, n2, n_tot = map(int, f.readline().split())
        n_unique, n_multi, _ = map(int, f.readline().split())
    return {"rsem_total_reads": n_tot, "rsem_alignable": n1, "rsem_alignable_pct": 100 * n1 / n_tot,
            "rsem_unique": n_unique, "rsem_multi": n_multi}


def parse_kallisto(path):
    with open(path) as f:
        d = json.load(f)
    return {"kallisto_processed": d["n_processed"], "kallisto_pseudoaligned": d["n_pseudoaligned"],
            "kallisto_pseudoaligned_pct": d["p_pseudoaligned"], "kallisto_unique_pct": d["p_unique"],
            "kallisto_version": d["kallisto_version"]}


rows = []
for s, star, rsem, kal in zip(samples, snakemake.input.star, snakemake.input.rsem, snakemake.input.kallisto):
    rows.append({"sample_id": s, "n_runs": len(runs[s]), **parse_star(star), **parse_rsem_cnt(rsem),
                 **parse_kallisto(kal)})
pd.DataFrame(rows).to_csv(snakemake.output.qc, sep="\t", index=False, float_format="%.4g")


def bench(path):
    """Snakemake benchmark: s = wall seconds, cpu_time = CPU seconds, max_rss = MB."""
    if not os.path.exists(path):
        return {"s": float("nan"), "cpu_time": float("nan"), "max_rss": float("nan")}
    return pd.read_csv(path, sep="\t").iloc[-1][["s", "cpu_time", "max_rss"]].to_dict()


rows = []
for s in samples:
    row = {"sample_id": s}
    trim = [bench(f"benchmarks/trimmomatic/{r}.tsv") for r in runs[s]]
    row["trim_wall_s"] = sum(t["s"] for t in trim)
    row["trim_cpu_s"] = sum(t["cpu_time"] for t in trim)
    for step in ("star", "rsem", "kallisto"):
        b = bench(f"benchmarks/{step}/{s}.tsv")
        row[f"{step}_wall_s"], row[f"{step}_cpu_s"], row[f"{step}_max_rss_mb"] = b["s"], b["cpu_time"], b["max_rss"]
    row["star_rsem_wall_s"] = row["star_wall_s"] + row["rsem_wall_s"]
    row["star_rsem_cpu_s"] = row["star_cpu_s"] + row["rsem_cpu_s"]
    rows.append(row)
pd.DataFrame(rows).to_csv(snakemake.output.compute, sep="\t", index=False, float_format="%.4g")
print(f"{len(samples)} samples summarized")
