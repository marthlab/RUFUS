# RUFUS.interpret — functional documentation and defect audit

Audit of `src/RUFUS.interpret.cpp` (7,388 lines) and its pipeline boundary, July 2026.
Produced by a nine-way parallel review covering every region of the file plus the upstream contract.

**Confidence marking.** Findings marked **[V]** were verified directly against source, compiled
test programs, or real run output. Findings marked **[R]** are reported by review and are
well-argued but were not independently re-checked. Treat [R] as a strong lead, not a fact.

---

# Part 1 — What RUFUS.interpret does

## 1.1 Position in the pipeline

`RUFUS.interpret` is the final stage. It is invoked from exactly one place,
`scripts/Overlap.shorter.sh:373`:

```
samtools view -h <contigs>.bam | perl AddSAtoReadSame.pl | grep -v chrUn \
  | RUFUS.interpret -mob … -mod … -mQ 10 -r <ref> -hf <HashList> -o <stub> \
      -m <MaxAlleleSize> <parentCRString> -sR … -s … -e … -rp … -ip … -plct …
```

Input is **assembled contigs aligned back to the reference**, not reads. By this point the
read→contig mapping has been discarded by the overlap assembler. Everything the caller knows
about allele support arrives as k-mer count tables.

## 1.2 The five inputs and what they populate

| Flag | File | In-memory table | Notes |
|---|---|---|---|
| `-hf` | `<proband>.k25_c<N>.HashList` | `Hash` (`map<string,int>`) | Subject-unique k-mers. **Also sets global `HashSize` from the first line's string length.** |
| `-s` / `-sR` | `…fastq.sample` / `…fastq.Ref.sample` | `MutantHashes` (merged) | Subject counts for contig k-mers and for reference k-mers |
| `-c` / `-cR` ×N | `ctrlhash.…<control>` | `ParentHashes[i]` (merged) | One map per control; **paired with `-c` by ordinal position only** |
| `-e` | `<stub>.ref.RepRefHash` | `ExcludeHashes` | Repeat/artifact veto |
| `-mod` | `<stub>.Jhash.histo.7.7.dist` | `DistGlobal`, `GenPrior` | Copy-number depth model — **see §2.1, this usually fails to load** |

Two extra channels are smuggled through the SAM rather than passed explicitly:

- **The QUAL string is not base quality.** `AnnotateOverlap` overwrites it with the per-base
  count of subject-unique k-mers covering that base, ASCII+33, **capped at 93**. `createPeakMap()`
  reads it as a signal and marks local maxima; those peaks are what designate a variant "DeNovo"
  and what gate nearly every SV breakpoint decision.
- **Strand counts live in the contig name.** `parse()` splits QNAME on `:` and reads fields 1 and 2
  as forward/reverse read counts. This is the *only* strand information in the program.

## 1.3 The core data model — `SamRead`

`getRefSeq()` (2821) expands each contig into a **column-aligned** representation. One column per
expanded-CIGAR position, and every vector shares that index space:

| CIGAR op | `RefSeq` | `seq` | `cigarString` out | `Positions` |
|---|---|---|---|---|
| `M` | ref base | read base | `M` if match, **`X` if mismatch** | `pos+i-InsOffset` |
| `I` | `-` | read base | `I` | previous ref position |
| `D` | ref base | `-` | `D` | ref position |
| `S` | `-` | read base | `S` | fabricated (extrapolated) |
| `H` | `H` | `H` | `H` | `-1` |
| `Y` | `-` | read base | `Y` | (tandem dup, injected by `BetterWay`) |

**The invariant** (nothing asserts it, and it is broken in several places — §3):

```
seq == RefSeq == cigarString == qual == strand == Positions == ChrPositions == PeakMap
    == AltKmers == RefKmers == MutAltCounts == MutRefCounts == MutContigCounts
    == MutHashListCounts == ParAltCounts[pi] == ParRefCounts[pi]     (all same length)
```

**The sentinel convention** in the count vectors:

| Value | Meaning |
|---|---|
| `> 0` | a real k-mer count |
| `0` | alt k-mer is identical to the ref k-mer (this column is not variant) |
| `-1` | a k-mer exists here but is absent from the lookup table |
| `-3` | no k-mer at this column (`getHash` returned `""`) |
| `< 0` in `MutContigCounts` only | a *negated real count* — collides with the `-3` sentinel |

Live consumers filter with `> 0`, which correctly excludes all three sentinels. `CheckPhase` is
the exception and is broken because of it (§3, D3).

## 1.4 Processing flow in `main()`

After setup, `main()` makes **six sequential passes** over one `vector<SamRead> reads`.
Contigs are never explicitly grouped into events; grouping is emergent via the `alignments` list
and, for multi-contig SVs, via **adjacency in the vector**.

| Pass | Lines | What it does |
|---|---|---|
| 0 | 5437–5461 | Pair alignments by QNAME string match (O(n²)); build `alignments` lists |
| 1 | 5462–5466 | `LookUpKmers`, `CheckPhase`, `clipPattern = ClipPattern()` |
| 2 | 5471–5660 | Multi-contig mobile elements → `<INS:ME:MOB>` |
| 3 | 5661–5907 | Multi-contig large DEL/DUP → `<DEL>`, `<DUP>` (**only >1 kb**) |
| 4 | 5909–5982 | `BetterWay` collapse, then `parseMutations` → all SNVs/indels |
| 5 | 5986–7377 | Six independent SV blocks: translocation/BND, inversion, triple-align insertion, large insertion, orphan MOB, `LastDitch` BND |

### `clipPattern` vocabulary

`ClipPattern()` collapses the expanded CIGAR into `{c, m}`, keeping only runs longer than ~10.
The SV code branches on this constantly:

- `"mc"` — aligned block then clip; breakend at the **right** end
- `"cm"` — clip then aligned block; breakend at the **left** end
- `"mc"`/`"cm"` as a *pair* is the canonical split-contig signature
- `"mcm"`, `"cmc"` — multi-junction; mostly unhandled
- `"mm"`, `"cc"` — **non-alternating**, produced when an intervening short run is dropped without
  resetting the run character. Accepted by the `length()==2` gates as if they were real pairs. [R]

### `parseMutations` — the small-variant caller (2402–2778)

Single loop `for (i = 25; i < cigarString.size() - 25; i++)`. A record is emitted iff:

1. contig passed `mapQual > MinMapQual` and `alignments.size() <= 2`;
2. `cigarString[i] ∈ {X, I, D, Y}` and `RefSeq[i] != 'N'`;
3. at least one column in the maximal run has encoded depth `> '!'`.

**No filter suppresses emission** — filters only decorate the ID and FILTER columns. The run is
extended forward while the op stays in `{X,I,D,Y}`, so `XIX` and `XDDX` become single complex
records. Indels get a left-anchor base by walking backward to the first column with a real
`ChrPositions` entry.

### `BetterWay` — split-contig collapse (3158–4128)

Takes exactly two `SamRead`s (primary + supplementary of the same contig) **by value** and returns
one synthetic read with the inter-alignment gap materialised as explicit `D` or `Y` columns — the
form `parseMutations` understands. So a split-read SV becomes an ordinary indel call. Six phases:

1. Build a 2-row column alignment, inserting gaps so both reads share an index space
2. Find the first anchored column
3. `bestQual` — restore the k-mer depth track destroyed by hard clipping
4. The merge walk: emit `D` runs for deletions, `Y` runs for tandem dups, classify oversize events
5. Rewrite internal clips to `I`; **discard the entire merge if ≥150 clip bases remain**
6. Opposite-strand branch: writes side-channel files only, produces no merged record

Because the input is by value, nothing written to `reads[1]` escapes; and the caller assigns the
result into a *local copy*, so the master `reads` vector is never updated either.

### Genotyping and AO/RO

Three independent paths that share almost no code:

| Path | Function | AO/RO estimator | Genotyper |
|---|---|---|---|
| SNV/indel (WGS) | `GetModes3` (1651) | filtered arithmetic **mean** (`PickDepthSomatic`) | `BayseanGenotyper` |
| SNV/indel (exome) | `GetModes` (1770) | off-center **median** `[(size-2)/2]` | `ShittyGenotyper` |
| SV/BND | `createStructGenotype` (879) | sorted `[0]` — the **minimum** | `ShittyGenotyper` |

`resources/vcf_header.txt` describes all of them as "Mode of … kmer counts", which is wrong in
every path. **AO/RO are k-mer depth statistics, not read counts.** QUAL is
`SupportingHashes / PossibleAltKmer * 100` — a percentage, not a Phred score.

The Bayesian model, as implemented: for each column of `DistGlobal` it **sums** per-k-mer
likelihoods rather than multiplying them, so evidence never accumulates (3 k-mers and 300 give the
same posterior shape) and a single repeat-inflated k-mer can outvote fifty consistent ones. The
normalizer `Pb` omits the prior, so `PaB` is not a probability — harmless for `argmax`, but the
posterior is then discarded anyway, which is why there is no GQ or PL anywhere in the output. [R]

---

# Part 2 — The findings that change the plan

These are the four things worth knowing before touching anything else.

## 2.1 The copy-number model is never built, in any run mode, and the genotyper is inert **[V]**

`ModelDist` writes **both** `.7.7.model` and `.7.7.dist` (`src/ModelDist.cpp:392, 403`), but it
never runs. The model phase is gated on `[ -z "$_arg_min" ] && [ $_arg_exome == "FALSE" ]`, and
`_arg_min` is initialised to `5` in the defaults block and **never cleared** — so the first test
is never true and the `else` branch always executes, hand-writing a 4-line `.7.7.model`
placeholder and never creating `.7.7.dist` at all.

This is **not** a property of `-min` or exome mode; passing `-min` changes nothing, because the
default already defeats the condition. The commented-out prior condition one line above the gate
is the smoking gun: originally every non-exome run built the model, and adding the
`-z "$_arg_min"` conjunct silently disabled it for everyone. The help text compounded it by
claiming `-m,--min` has "no default", which is false.

Confirmed across **every** run in `resources/reg_test_files/runs/`: zero `*.7.7.dist` files,
fifteen `*.7.7.model` placeholders, zero logs containing "Starting model phase", six containing
"min was provided". The actual behaviour is now documented inline at the gate, at the `_arg_min`
default, and in both help blocks in `runRufus.sh`.

`Overlap.shorter.sh:373` unconditionally passes `-mod <stub>.Jhash.histo.7.7.dist`.
`ProcessDist` (4755) treats a failed open as **non-fatal** — it prints and returns, leaving
`DistGlobal` and `GenPrior` empty.

Verified in the `chr20_m5` regression run:

```
chr20_17559175.out:173:  Error no model file given, not worring abou this now
on disk: …Jhash.histo.7.7.model     (exists)
         …Jhash.histo.7.7.dist      (absent)
```

Downstream consequence, verified on that run's final post-processed VCF (1,137 records):

```
GT     distribution:  1137  "."      ← zero genotypes
FILTER distribution:  1137  "."      ← zero PASS
```

**Scope correction (added after building the replay harness).** Replaying interpret
directly on the same inputs gives 1,680 raw records, of which **3 carry `FILTER=PASS` and
5 carry real genotypes** (4× `0/1`, 1× `1/1`). So the failure is confined to the
**SNV/indel path**, which routes through `BayseanGenotyper` and therefore depends on the
model file. The SV/BND paths use `ShittyGenotyper` — pure arithmetic on two counts, no
model — and set `FILTER=PASS` independently, so they still produce genotypes and passing
calls. Post-processing then removes those few, which is why the final VCF looks totally
empty of both. Read every claim in this section as "in the SNV/indel path". Reproduce with:

```
./tests/replay/replay.sh run ../resources/reg_test_files/runs/chr20_m5/rufus_chr20
```

The cascade: `DistGlobal` empty → `sums` empty → `maxI = -1` → no genotype branch matches →
`ParseGenotype("","")` returns `"."` → line 2609 sees no `1` in the genotype and overwrites
`Denovo = "Mosaic"` on **every** variant, erasing the PeakMap-derived DeNovo signal → line 2670's
`PASS` condition (`Denovo == "DeNovo"`) becomes unreachable.

So three of the most visible output defects — no genotypes, everything labelled Mosaic, no PASS —
are one plumbing failure. `Dist1XCutoff` also falls back to `100000`, which disables
`PickDepthSomatic`'s repeat filter; this is exactly the value the author's TODO at line 1451
recorded observing in real COLO829 output without following up.

**Important caveat before fixing:** simply making the model load will *not* produce correct
genotypes. It activates three currently-dormant defects — the shadowing bug at 1586 (which
specifically breaks homozygous calls, since `maxI ≥ 3` ⇔ ≥2 copies), out-of-bounds reads at
1488/1509/1539, and unconstrained GT ploidy in `ParseGenotype`. Those must be fixed in the same
change or the first real model file makes output worse.

## 2.2 There is no "is this allele present in the control?" check for SNVs **[R]**

For a tool whose premise is subtracting the control, this is the biggest functional gap.

The genotype-based comparison is **commented out** at 2693–2696. Per-parent genotypes *are*
computed (2512–2522, at the cost of four full vector copies per parent per variant) and *are*
emitted in the sample columns — and then used for nothing.

What remains is a k-mer heuristic at 2583 that requires `parentCounts[k][j] ≤ ParLowCovThreshold`
(default 7). It has an **upper** bound, so a control carrying the variant at normal depth is
invisible to it. The review reports 129/1675 emitted calls in a production VCF with control
`AO > 7`, 17 of them carrying no filter at all.

The SV path uses a different function (`SVCheckParentsForLowCov`) with the **opposite** convention:
SNVs sum low counts *across* controls, SVs *reject* anything low in more than one control.

## 2.3 A large class of coordinate and VCF-validity defects, partly masked by downstream band-aids

`runRufus.sh:1496–1585` contains an awk sanitizer that drops records with empty/invalid REF or ALT,
a tabix retry loop that deletes up to 50 more malformed records, and
`bcftools +fill-from-fasta -c REF` which **overwrites REF from the reference FASTA**. Those exist
because interpret emits exactly those defects — and the `fill-from-fasta` step converts a
detectable inconsistency into a plausible-looking, undetectable wrong call.

Concretely: REF is fetched one base to the right of POS at 5 of 13 SV emission sites (the other 8
are consistent, so this is a divergence, not a uniform convention error); every non-PASS FILTER
value carries a trailing `;` producing an empty filter ID; `FEX=` embeds raw semicolons from
`filterSV()`; `EN=` emits six comma-separated values plus a 50 bp raw sequence under a
`Number=1,Type=String` declaration; `<DEL>`, `<DUP>`, `<INV>`, `<INS>` and `FILTER=fail` are
undeclared in the header. **100% of BND records are unmated** — the ID column is a composite
string like `OrphanBND-LC=0bnd_1-DeNovo` while MATEID is the bare `bnd_2`, so they can never
match, and the reciprocal `LastDitch` call reassigns both BNDids anyway. [R]

**Verified on a production run (SMHT004 3A, 300× whole-genome, `-m 5`, DSA + tech control, Jul 2026) [V]:**
the undeclared `FILTER=fail` is real and reaches downstream tools — `bcftools` emitted
`[W::vcf_parse_filter] FILTER 'fail' is not defined in the header` on the RUFUS.Interpret output,
and the `runRufus.sh` sanitizer dropped **13 of 2662 records** as malformed (2649 kept). The
surviving genotyped VCF (`gx.wg.vcf.gz`, 1,354 calls) carries `FILTER=.` on **100%** of records —
i.e. §2.1's inert-genotyper failure also reproduces here. So both the §2.3 header-validity defects
and the §2.1 model-absent cascade reproduce on real whole-genome data, not just the `chr20_m5`
regression. (This run did pass `-min`, but per the corrected §2.1 that is immaterial: the model
phase is unreachable regardless, so `-min` versus whole-genome was never the distinguishing
factor.)

## 2.4 Two independent k-mer extraction passes with different correctness properties **[R]**

`LookUpKmers` (2948) and `BuildUpHashCountTable` (1310) both walk the contig building k-mer/count
arrays, with different masking rules and different sentinels. `GetModes3`/`GetModes` consume the
first; the parent low-coverage heuristics consume the second. Only the second guards k-mer length.

This matters because `HashToLong` (100) is **length-blind** (verified by compiling it): unused bits
are zero and `00` encodes `A`, so a truncated 20-mer hashes identically to the full 25-mer with
`AAAAA` appended. `getHash` returns short strings at contig ends, so the trailing `HashSize-1`
columns of every contig look up whatever count the A-padded k-mer happens to have. It is also
silently degenerate above k=32 with no validation on `-hs`, and non-ACGT characters (including
**lowercase** — i.e. any soft-masked reference) encode as `A`.

Consolidating these two passes should precede any AO/RO work.

---

# Part 3 — Defect inventory

Ranked within each tier. `file:line` refers to `src/RUFUS.interpret.cpp` unless stated.

## Tier 1 — silently wrong output

| # | Location | Defect |
|---|---|---|
| A1 | `runRufus.sh:1177` + `4759` | Model file absent in `-min`/exome runs; failure is non-fatal → genotyper inert **[V]** |
| A2 | `1586–1590` | Variable shadowing leaves `C` uninitialized when `maxI > 2` → garbage AO/RO at copy number >2 **[V]** |
| A3 | `2840–2853` | Leading `H` followed by `S`: only the hard clip is subtracted, so **every coordinate in the contig shifts**. Author's own comment at 2900 flags it **[R]** |
| A4 | `2693–2696` | Parent genotype comparison commented out; surviving check has an upper bound (§2.2) **[R]** |
| A5 | `2609` / `2670` | `Denovo` overwritten to `Mosaic` on every record; `PASS` unreachable **[V]** |
| A6 | `1678` vs `1686` | AO filtered against `ExcludeHashes`, RO not → AF biased wherever the artifact control touches a locus **[V]** |
| A7 | `2454–2462` | Indel anchor accepts `I`/`D`/`S` columns → REF or ALT becomes `"-"` or empty; record dropped downstream and **the real indel is lost** **[R]** |
| A8 | `1431` vs `1651` | QUAL numerator/denominator use different windows and filters → QUAL > 100 **[V]** |
| A9 | `2419` + `2467` | Interior non-ACGT bases silently dropped from REF/ALT, shortening the allele **[R]** |
| A10 | `1826` | Exome path filters *parent* alt counts using the *reference* k-mer (copy-paste) **[R]** |
| A11 | `3200`, `3221` | `BetterWay` indexes `reads[B].cigarString` with `Acount` → rows de-phase, heap over-read **[R]** |
| A12 | `3897` | `BetterWay` hard-codes `'M'` for B's cigar → **all SNVs/indels on the B side of a collapsed split read are never called** **[R]** |
| A13 | `3432` | Deletion REF fill indexes by column `i` instead of reference coordinate `j` **[R]** |
| A14 | `2137–2158` | `CheckPhase` misreads the `0` sentinel; all branches dead → `PH=none` on every record **[R]** |
| A15 | `2306` | `tempPeakMap[i] == tempPeakMap[i-1];` — comparison not assignment; the deletion fixup never runs **[R]** |

## Tier 2 — crashes, UB, out-of-bounds

| # | Location | Defect |
|---|---|---|
| B1 | `100–121` | `HashToLong` length-blind, non-ACGT→`A`, silent aliasing above k=32, `-hs` unvalidated **[V]** |
| B2 | `4783–4791` | `ProcessDist` writes one past the end of `DistGlobal`; first data row stored one column off **[R]** |
| B3 | `1488`, `1509`, `1539` | Depth clamp uses `>` where `>=` needed, then indexes at exactly `size()` **[R]** |
| B4 | `6371`, `6582`, `6820` | Mutating loop variable `i` inside the inner `j` loop → `alignments[1]` out of bounds **[R]** |
| B5 | `2414`, `1313`, `2289` | `size() - 25` / `size() - HashSize` unsigned underflow on short contigs **[R]** |
| B6 | `3277`, `3314` | Unhandled CIGAR op doesn't advance the cursor → infinite loop + unbounded growth **[R]** |
| B7 | `2910` | `=`/`X`/`N`/`P` unhandled: `N`/`P` walk `seq` index past the end; `=`/`X` silently drop every base **[R]** |
| B8 | `3443` + 12 sites | `PeakMap[Abreak-1]` evaluated before the `Abreak > 0` guard (`and` short-circuits left to right) **[R]** |
| B9 | `6918` | `substr` with `pos = -1` → `SIZE_MAX` → uncaught `out_of_range` **aborts mid-VCF** **[R]** |
| B10 | `3352–3364` | Unbounded scan for first anchored column **[R]** |
| B11 | `3076–3083` | Guard is `size() >= 2` but body reads `temp2[2]` **[R]** |
| B12 | `2833–2837` | `getRefSeq` early-returns on unknown contig but the read is still admitted → empty `Positions`, OOB later. Directly relevant to chromosome-sharded runs **[R]** |
| B13 | `4914–4982` | Every flag reads `argv[i+1]` with no bounds check **[R]** |
| B14 | `5107–5109` | `-cR` loop writes into `ParentHashes[i]` sized by the `-c` count **[R]** |

## Tier 3 — silent-failure / operability

- **Fatal errors return 0.** Unknown flag, missing `-r`, missing `-hf`, unopenable SAM all exit
  successfully (`4987, 4994, 4998, 5276`), and interpret runs at the end of a pipe. The driver
  reports "no variants found for this region". Note `-hS` in the help text vs `-hs` in the parser. [R]
- **No `is_open()` check on any of the six hash inputs** (`5036, 5076, 5101, 5126, 5138, 5156`).
  An empty `ParentHashes` makes every variant look de novo, with no warning. `CheckJellyHashList.sh`
  wraps `jellyfish query` in `timeout 1h` with no `set -e`, so a timeout writes a **truncated file
  and reports success**. [R]
- **Early `return 0` before the `#CHROM` line** (`5404`) emits a VCF with meta-lines and no header
  line for any empty shard — `bcftools` rejects it. [R]
- **SV records emit one sample column** while the header declares `1 + controls` (`4744` et al.),
  so any trio/tumour-normal run that emits an SV produces a ragged VCF. [R]
- **`operator[]` used for read-only lookups** on `ExcludeHashes`, `Hash`, `ParentHashes`
  (`768, 829, 1672, 1691, 1802, 1838`) — inserts a zero entry on every miss, growing the tables
  without bound. Worse, at `1684` the membership gate is `Hash.count(...)` while `1691` does
  `Hash[...]`, so **AO becomes order-dependent**: a k-mer inserted by one variant passes the gate
  for a later one. [R]
- **`##fileDate` is epoch seconds**, not `YYYYMMDD` (`5301`). [R]
- Non-unique ID column throughout; `SVTYPE=TRANS` and `COPY:PASTE` are not valid VCF types. [R]

## Tier 4 — dead, stubbed, or never wired up

- **`src/SamRead.h` / `src/SamRead.cpp` are a second, divergent copy of the `SamRead` class that
  nothing compiles** (absent from `CMakeLists.txt`). A fix applied there will appear to do nothing.
  Worth resolving as part of the CI/CD consolidation. **[R]**
- `processMultiAlignment()` (2030) — empty stub, `//check if this is a mis-joined contig`. Never
  called. Mis-joined contigs are exactly what produces the dominant `PA` filter.
- `FixTandemRef()` (2798) — **never called** (only call site commented out at 4116). So `RefSeq`
  over `Y` columns stays `-` and every tandem dup is emitted as a plain insertion, built on a
  knowingly-wrong reference. The whole DUP:TANDEM path operates on it.
- `StructCall` is computed by `compressVar` and **omitted from the live output statement** — all
  `SVTYPE=DUP`/`END`/`SVLEN` annotation for tandem dups is computed and thrown away. It is also
  never cleared between variants, so uncommenting the alternative output line at 2740 would
  immediately mis-annotate SNVs as duplications.
- The mmap/binary-search block (`checkPage`/`ProcessPage`/`search`, 123–310) is entirely dead and
  contains a use-after-`munmap` — it has never successfully run. Delete it. **[R]**
- `CheckTranslocation`, `CheckMob`, `CheckPolyATail`, `CheckLargeInsert` (4299–4313) — all
  `return false;`, none called. Note `CheckMob` shadows the live member `SamRead::checkMob` by
  capitalization only.
- `StartsWithAlign`, `EndsWithAlign`, `StartsWithAlignAtPeak`, `EndsWithAlignAtPeak` — defined,
  never called, and each carries out-of-bounds indexing.
- `GetModes2` (1719–1769) — fully commented-out earlier draft of `GetModes3`, plus a live
  declaration and a commented call. Delete.
- `PickDepth` never called; `PickDepthAverage` reachable only via the shadowed-`C` branch.
- `-w`/`isWindowed` parsed and never read; `-hs` always overwritten at 5218; `-plct` undocumented.
- `totalDeleted`, `totalAdded`, `ScGlobal`, `candidateHash`, `MaxBND`, `UsedForBigVar` — never read.
- `MutHashListCounts` and `MutContigCounts` — computed at cost per base per contig, read only by
  the never-called `writeVertical` and by the inert `CheckPhase`.
- Entropy (`w1`–`w5`) is emitted into `EN=` and never thresholded or used.
- **Six auxiliary output files** (`BEDOutFile`, `BEDBigStuff`, `BEDNotHandled`, `Invertions`,
  `Translocations`, `Unaligned`) — nothing in the repository reads any of them. Two are
  misleadingly named: `Invertions` and `Translocations` are written only by `BetterWay`, never by
  the SV classifier that actually emits `<INV>` and BND records. **[R]**
- All SV classification inside `BetterWay` (inversions, translocations, mobile elements, oversize
  events) writes to those side channels and **produces no VCF output at all**.

## Hard-coded values that should be parameters

`0.997` (dist coverage), `100000` (no-model `Dist1XCutoff`), `GenPrior = 1/i`, `< 400` AO/RO
ceilings (×4, with the author's own "should be based on cov" comment), `150` (BetterWay merge
discard), `MaxVarentSize + 1000` (×8), `25` in `parseMutations`' loop bounds (should be
`HashSize`), `10` (clipPattern run length), `NumLowCov > 25`, `"hs37d5"` hard-coded as the decoy
contig in five branches (silently never fires on GRCh38/T2T), strand-bias thresholds of
`0.99/0.01` in `filterSV` versus `0.99999/0.00001` in `parseMutations` — the two paths disagree
with each other, and neither matches the 95/5 + DP≥30 model from the strand-bias work.

---

# Part 4 — Upstream data needs

This is the question that determines whether fixes are local or structural.

## 4.1 Fixable entirely inside interpret

Everything in Tier 1 A2/A5/A6/A7/A9/A10/A13, all of Tier 2, all of Tier 3, all of Tier 4.
Also: reference-orient `SB` using the SAM `0x10` flag before emitting it; validate `HashSize <= 32`
and cross-check k across all hash files at load; replace `operator[]` with `find()` on every table.

## 4.2 Requires an upstream change, ranked by value per unit cost

**Rank 1 — make the model file reach interpret, or make its absence fatal.**
`runRufus.sh` must either run `ModelDist` in `-min`/exome mode or pass an explicit
`--no-model` so that "model missing" and "model deliberately not used" are distinguishable from
"model failed to build". Today all three are silent. Cost: a few lines. This is the single
highest-value change in the list, but see the §2.1 caveat — fix A2/B3 and GT ploidy in the same
change.

**Rank 2 — control k-mers counted at `-L 1`, or a "was-filtered" marker.**
`jellyfish count -L 2` on controls means a k-mer seen **exactly once** in a parent does not exist
in the control table. `jellyfish query` returns 0, interpret reads "not in parent", and it
contributes to a de novo call. The whole `ParLowCovThreshold` band (1..7) can never observe count 1
at current defaults. Cost: one default change, or a cheap second query pass restricted to the
contig k-mer set. **No new file format.** [R]

**Rank 3 — stop capping counts at `MaxCov=100000` silently.**
k-mers above the ceiling are *deleted from the file*, so interpret sees them as absent rather than
as saturating repeats — a high-copy repeat k-mer looks like clean de novo evidence. Emit a
saturation sentinel instead and teach interpret to treat it as "repeat, do not call". [R]

**Rank 4 — parse `SA:Z:`.** It is already in the input (`AddSAtoReadSame.pl` puts it there) and
interpret ignores it entirely — `parse()` scans optional fields for `AS` only. It states directly
what several code paths spend hundreds of lines re-deriving by heuristic: which alignments belong
together, and where the clipped sequence goes. It would replace the O(n²) QNAME matching in Pass 0,
the vector-adjacency window scans in Passes 2/3/5 (currently `j ∈ [-2,2]`, so contig pairing
depends on the aligner emitting the two halves of a junction within two lines of each other), and
the `H`-vs-`S` position asymmetry. **This is the highest-value item for SV quality.** [R]

**Rank 5 — align with `bwa mem -Y`.** Hard-clipped supplementary alignments physically lack their
bases, so `getRefSeq` writes literal `'H'` into `seq`/`RefSeq` and every k-mer overlapping the clip
becomes `-3`. A supplementary alignment therefore contributes **no k-mer evidence near its own
breakpoint** — the most informative region. `-Y` supplies the bases and also removes the entire
`bestQual` phase in `BetterWay`. Cost: one flag. [R]

**Rank 6 — pass real read-level support.** Two tiers:
- *Cheap:* `F` and `R` in the QNAME are already a read count (`F+R` = reads assembled into the
  contig), and interpret uses them only for the strand ratio. Reporting `F+R` in INFO needs no
  upstream change. [R]
- *Real:* a sidecar `contig → [read, offset, orientation]` from the overlap assembler. This is the
  **only** route to true AO/RO, and it simultaneously enables position-level strand bias. Largest
  cost, largest quality win.

**Rank 7 — preserve per-base assembly depth.** `AnnotateOverlap` overwrites the depth-encoded QUAL
with the k-mer count, and the ASCII channel caps at 93 — already below the depth of the samples
this runs on. `'!'` doubles as both "depth zero" and "no data". An integer array tag would remove
the ambiguity and the ceiling.

**Rank 8 — replace filename-as-API with a manifest.** Sample names are reverse-engineered by string
surgery: the subject from `-o` (truncate at `.generator`, strip trailing `.chr`), controls from the
literal marker `"overlap.asembly.hash.fastq."`. In whole-genome mode there is no `.chr`, so the
same subject gets a **different VCF sample name** in WG versus region runs. `Overlap.shorter.sh`
carries a comment explicitly documenting that it is named to accommodate this C++ parsing quirk.
A small manifest (subject name, control names in `-c` order, k, cutoffs actually applied, region)
retires that coupling and makes the other fixes safe to ship. This is the same class of problem as
the 255-char filename bug already fixed. [R]

**Also needed:** hash files carry no header, so interpret *infers* k from the first line's length,
*assumes* the separator, and *assumes* canonicalization. A one-line header recording k, canonical
yes/no, min/max coverage and producer version would turn a whole class of silent-wrong-answer
failures into startup errors. And `MaxBND`/`CurrentSVeventID` are per-process globals starting at
0, so BND and SV event IDs **collide across shards** when per-shard VCFs are concatenated.

---

# Part 5 — Suggested sequencing

0. **Replay harness — done.** `tests/replay/` re-runs interpret standalone against a
   preserved run directory in ~3 minutes and diffs VCF, scraped log signals and exit code
   against a frozen baseline. `tests/replay/baselines/chr20_m5` is frozen and verified
   deterministic. Run `./tests/replay/replay.sh check <run_dir> <baseline>` after every
   change; exit 2 means behaviour moved. Note the chr20_m5 baseline exercises the
   *model-absent* path only — a baseline from a run with a real `.7.7.dist` is needed
   before touching `BayseanGenotyper`.

1. **Instrument first.** Turn the silent failures into loud ones: `is_open()` checks on all six
   inputs, non-zero exit on fatal paths, validate `HashSize <= 32` and k-consistency across hash
   files, fail on an unknown contig. Nothing else can be trusted until a broken run looks broken.
2. **Fix the model plumbing plus its three dormant dependents together** (A1 + A2 + B3 + GT ploidy).
   Expect the output to change substantially — genotypes and PASS will appear for the first time.
3. **Consolidate the two k-mer extraction passes** and fix `HashToLong` length handling (§2.4).
   This is a prerequisite for any trustworthy AO/RO work.
4. **Then** the AO/RO and QUAL fixes that started this investigation (A6, A8) — they are cheap once
   3 is done, and premature without it.
5. **Coordinate correctness pass** (A3, A7, A9, A13, the REF/POS off-by-ones) with a round-trip
   test against the reference, so the downstream `fill-from-fasta` band-aid can be removed.
6. **Delete the dead code** (Tier 4) before restructuring anything — roughly 1,500 lines, including
   one whole uncompiled duplicate class.
7. **Then** decide on the upstream changes. `SA:Z:` parsing (Rank 4) and `bwa mem -Y` (Rank 5) are
   the two that unlock the most, and neither requires a new file format.
