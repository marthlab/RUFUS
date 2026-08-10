# RUFUS.interpret replay harness

Re-runs `RUFUS.interpret` on its own against a preserved run directory, and diffs the
result against a frozen baseline. A full run takes **~3 minutes** and touches none of
the 200 GB CRAMs or the multi-hour assembly pipeline.

The point is not to assert that current output is *correct* — much of it isn't (see
[`docs/RUFUS.interpret.audit.md`](../../docs/RUFUS.interpret.audit.md)). The point is that
once behaviour is frozen, every subsequent change produces a diff you either explain or
investigate, and abandoning a direction is `git revert` plus a re-check rather than
archaeology.

## Why this works

`RUFUS.interpret` is the last pipeline stage, and every input it takes is a file earlier
stages leave behind in the work directory. `chr20_m5` has all of them: the contigs BAM,
the `-s`/`-sR` sample hashes, both `ctrlhash` control files, the `RepRefHash`, the
`MOB.sam`, and the `HashList`. So interpret can be replayed standalone.

## Usage

```bash
export PATH=/uufs/chpc.utah.edu/sys/installdir/samtools/1.16/bin:$PATH
cd <repo root>

# Check current code against the frozen baseline. This is the one you'll run most.
./tests/replay/replay.sh check ../resources/reg_test_files/runs/chr20_m5/rufus_chr20 \
                               tests/replay/baselines/chr20_m5

# Replay once and keep the output for inspection
./tests/replay/replay.sh run   <run_dir> [out_dir]

# Re-freeze after an intentional, reviewed behaviour change
./tests/replay/replay.sh freeze <run_dir> tests/replay/baselines/chr20_m5

# Prove determinism (run twice, diff). Do this before trusting a new baseline.
./tests/replay/replay.sh twice  <run_dir>
```

Exit status: `0` match, `1` usage/setup error, `2` baseline mismatch — so `check` drops
straight into CI or a pre-commit hook.

`<run_dir>` is the RUFUS *work directory* of a completed run, e.g.
`../resources/reg_test_files/runs/chr20_m5/rufus_chr20`.

Environment overrides: `RUFUS_INTERPRET` (binary), `RUFUS_REPLAY_REF` (reference FASTA),
`SAMTOOLS`.

## What gets compared

| Artifact | Why |
|---|---|
| `calls.vcf` | The VCF with `##fileDate` and `##RUFUSCommandLine` stripped — the only two lines that change between identical runs. Record order is preserved, because order changes are real changes. |
| `signals.txt` | Error/warning lines scraped from interpret's stdout, with long numbers masked, deduped with counts. **This is the important one.** Interpret prints megabytes of debug and returns 0 even on fatal errors, so neither the raw log nor the exit code is usable alone. This catches "the model file stopped loading" — a change that alters thousands of records but might not look like a crash. |
| `exit_code` | Cheap, and will become meaningful once the `return 0`-on-fatal-error bug is fixed. |
| `summary.txt` | Record count, FILTER/GT distributions, QUAL stats, column count. Not a pass/fail gate — it's what makes a diff readable at a glance. |

The shadow work directory is built from **symlinks**, so the reg-test data is never
written to. It is the oracle; keep it pristine.

## Baselines

`baselines/chr20_m5/` — COLO829 tumour/normal, chr20, `-min 5`. Frozen 2026-07-28 from
`bin/RUFUS.interpret` as built 2026-06-08, on branch `docker`.

Verified deterministic: two consecutive runs produced byte-identical output.

Baseline at freeze time:

```
records:   1680
FILTER:    1156 "."   201 PA;   163 SB;   48 LCH;SB;   31 LCH;   ...   3 PASS   1 LMQ
GT:        1675 "."   4 0/1   1 1/1
QUAL:      mean=78.56  max=108.0  >100=4  missing=0
columns:   11 on every record
```

Two things worth knowing about those numbers:

- **`>100=4`** reproduces the QUAL overflow (audit A8). Max is 108.
- **3 PASS and 5 real genotypes.** The audit's §2.1 said PASS was unreachable and GT was
  `"."` everywhere; that is true of the **SNV/indel path** (which routes through
  `BayseanGenotyper` and needs the absent model file) but *not* of the SV/BND paths, which
  use `ShittyGenotyper` and set `FILTER=PASS` independently. The post-processed VCF shows
  1137 records with no PASS at all, so downstream filtering removes those 3. Raw interpret
  output and final output are different things — this harness measures the former.

## Adding a baseline

Any completed run directory works. Good candidates once available: a trio run (exercises
multi-control pairing) and a case where the `.7.7.dist` model **is** present (exercises
`BayseanGenotyper`, which the chr20_m5 baseline leaves completely untested).

Run `twice` first. If a run is non-deterministic, a baseline from it is worse than none.
