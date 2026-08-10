# Accuracy test (Tier 3) — chr20 COLO829

The **pre-release accuracy gate**: run RUFUS on a real matched tumor/normal pair over chr20, then score
recall/precision against a truth set. This is the check that runs on the `:stage` image before tagging a
release. (The fast, self-contained tests live in [`../functional/`](../functional); this one needs real
data.)

## Why the data isn't in the repo

The inputs are **large** (~100 GB CRAMs, a ~3 GB reference) and **controlled-access** (SMaHT data under a
data-use agreement). So only the **code** lives here; the **data** stays external and is referenced through
`config.sh`. Do **not** commit the CRAMs, reference, or any slice of them — even a chr20 subset is still
~GBs and still controlled-access.

## What you need (data manifest)

Lay these out under one directory (`RUFUS_TEST_DATA` in `config.sh`):

```
<RUFUS_TEST_DATA>/
  crams/
    SMHTCOLO829T-...-GRCh38.aligned.sorted.cram   (+ .crai)   # tumor  (SMaHT COLO829T)
    SMHTCOLO829BL-...-GRCh38.aligned.sorted.cram  (+ .crai)   # normal (SMaHT COLO829BL)
  ref/
    GCA_000001405.15_GRCh38_no_alt_analysis_set.fa (+ .fai)   # the reference the CRAMs were aligned to
  truth_set/
    SMaHT_COLO829_SNV_truth_set_v1.0.vcf.gz        (+ .csi)   # SNV truth set (SNV-only; carries RGN_T tiers + VAF)
```

| File | Source | Notes |
|---|---|---|
| COLO829 T / BL CRAMs | SMaHT data portal | Controlled-access — obtain under your DUA. Filenames encode the accessions; a matched pair (same center/aligner/build) is required. |
| GRCh38 reference | Public (`GCA_000001405.15_GRCh38_no_alt_analysis_set`) | The **no-ALT** analysis set. The CRAMs' `@SQ` M5s must match it, or decode fails. |
| COLO829 SNV truth set | SMaHT | Provides `INFO/RGN_T` confidence tiers (Tier1 = high) and per-variant VAF. |

The image (`RUFUS_SIF`) is a pulled `.sif` of the container under test — pull it first:
`apptainer pull rufus_stage.sif docker://stefinfection/rufus:stage`.

## Running it

```bash
cd tests/accuracy
cp config.sh.example config.sh        # then edit config.sh: data paths, RUFUS_SIF, threshold
#                                       also edit --account/--partition in run_chr20.slurm for your cluster
sbatch run_chr20.slurm                # ~34 min on the -vs config; outputs go to $RUFUS_TEST_DATA/runs/
bash   compare.sh                     # prints recall/precision + a PASS/FAIL verdict; exit 0 = pass
```

`compare.sh` output ends with a verdict line, e.g.:

```
  m5    recall 104/118 = 0.881    precision 887/1089 = 0.815    [PASS]
==================== VERDICT: PASS  (1/1 met recall >= 0.85) ====================
```

## Things to know

- **`-vs` is required.** Full assembly (no `-vs`) **deadlocks** on real chr20 data. `-vs` is also the
  shipped `setup_slurm` config, so it's what users actually run — this measures the real product.
- **The recall number is a moving baseline.** It reflects the current interpreter; strand-bias and
  genotyper work will shift it (especially precision). Re-baseline `MIN_RECALL` when the interpreter
  settles rather than treating today's threshold as final.
- **Precision is a lower bound.** The SNV truth set is a high-confidence subset, not exhaustive, so real
  calls outside it count against precision. A confident-region BED would turn it into a true figure.
- **Major releases:** this chr20 run is the routine proxy. For a whole-genome check, use the sharded
  `setup_slurm` path (not this direct-exec script).
