# SGSC load profiles on LV12_248

Measured Smart Grid Smart City (SGSC) household profiles assigned to the 121 customers of
LV12_248 (Geth et al. 2025), ready for OpenDSS and for the OPF. Built 12 Oct 2026.

## Files

| File | What it is |
|---|---|
| `make_lv12_loads.py` | Builds everything in `out/` from the SGSC JSON and your cleaned `Loads.dss` |
| **`lv12_opf.jl`** | **The OPF script for LV12_248 (Julia/PMD). Staged: checks itself against OpenDSS, then HC₀, HC₁, HC₂** |
| `lv12_diag.jl` | Diagnostic for a failed power-flow stage. `include` it in the same Julia session after `lv12_opf.jl` stops: checks parsed impedances against OpenDSS, solves with PMD's own model, with zero load, and warm-started |
| `hc0_oracle.py` | HC₀ check value in OpenDSS (unity PF, no volt-var, PV on all 121 customers), by bisection |
| `out/sgsc_profiles_hourly.csv` | 363 usable profiles × 24 clock hours (kW), with window start, season, weekday/weekend |
| `out/lv12_assignment.csv` | Which profile each customer gets: draw 0 (cyclic) and draws 1–100 (random) |
| `out/lv12_load_kw_draw000.csv` | Draw 0 flattened: each customer's 24 hourly kW |
| `out/Loads_h13_draw000.dss` | Drop-in `Loads.dss` for a 13:00–14:00 snapshot, draw 0 |
| `out/Loadshapes_draw000.dss`, `out/Loads_daily_draw000.dss` | Draw 0 as 24-hour daily loadshapes for OpenDSS daily mode |
| `out/oracle_draw000.txt`, `out/oracle_draws_h13.txt` | HC₀ results below |

## Method

- **Assignment** follows CSIRO's method (`load_matrix_to_dict` in the LVFT repo), where load i
  gets profile i, cyclically.
  - Load order is the order of `Loads.dss`.
  - Profile order is the order of the JSON file. CSIRO's Julia code iterates a `Dict`, so its
    order is arbitrary; here it is fixed.
- **Draw 0** is cyclic. **Draws 1–100** each take 121 distinct profiles at random (seed 2026).
  Report results as a distribution over draws, not from one draw.
- **Clock alignment, per customer.** Each SGSC customer in this file covers a *different* 24-hour
  window (see Corrections below). Each profile is placed on clock time using its own timestamps;
  the day is stitched from the two calendar days its window spans.
- **Hour convention.** Timestamps are taken to mark the end of each half-hour, so hour H is the
  mean of the half-hours stamped H:30 and H+1:00. The opposite convention shifts everything by
  30 minutes, which matters little around noon.
- **Excluded profiles:**
  - 4 windows that are not 48 consecutive half-hours.
  - 41 profiles that are zero all day (vacant, or not recording). `--keep-zero` keeps them.
- **Power factor** 0.95 lagging on every load (`--pf`). SGSC records active power only.
- **Season.** The usable pool is 31 summer, 124 autumn, 184 winter and 24 spring windows.
  `--season winter` (or autumn) restricts the pool. Summer has too few profiles for 121
  customers without repeats.

Regenerate:

```
python3 make_lv12_loads.py smartgridsmartcities.json <cleaned>/Loads.dss out --draws 100 --seed 2026
python3 hc0_oracle.py <cleaned>/Master.dss out --draws 0 --hours 8-17 --irr irradiance.csv
```

## Corrections to what the project notes assumed

1. **The SGSC file is not one 13–14 Dec 2012 window starting at 09:30.** That is only the first
   customer. The 408 windows run from Oct 2011 to Jul 2013, start at 32 different times, and are
   mostly autumn and winter.
2. **The earlier alignment scrambled clock time.** The earlier hourly profile (and CSIRO's
   notebooks, which index column 0 as midnight) used column positions as if every profile started
   at the same time. So the earlier "0.408 kW/customer at 13:00" is a mix of different clock
   times. Aligned properly, the shape is a normal residential day: overnight low near 04:00,
   morning peak around 08:00, evening peak 18:00–21:00. The midday level is similar, so earlier
   HC conclusions are shifted rather than overturned.
3. **Zero load is not a strict lower bound on HC here.** With voltage limits measured to the
   neutral, load on the other phases shifts the neutral and can raise the exporting phase's
   voltage. 1 of 101 draws (draw 82) fell 1.6% below the zero-load value. Zero load is a close,
   not a rigorous, lower bound.

## Validation

**Rebuilt cleaned feeder against `LV12_248_cleaning_log.md`** (base case: 1 kW at pf 0.88 on
every load):

| | This rebuild | Logged |
|---|---|---|
| Max neutral voltage | 1.310 V | 1.312 V |
| Phasor sum of MEN bond currents | 2.01 A | 2.013 A |
| Source earth current | 3.21 A | 3.22 A |
| Source power | 122.19 kW | 122.43 kW |
| Customer voltages | 239.5–243.7 V | 239.2–243.4 V |

Under the old (misaligned) loads, HC₀ comes out 2–3% higher than the logged values: 2.78 vs
2.71 at 13:00, and 2.74 vs 2.67 at zero load. Two causes were tested and ruled out: measuring
the limit to ground (gives 3.38) and keeping the full transformer impedance (gives 2.39). The
difference is somewhere in how the source is modelled. **So compute the check value with
`hc0_oracle.py` on your own cleaned folder**: the oracle must come from the exact files the OPF
reads.

**Loadshape files:** OpenDSS daily mode at hour 13 gives exactly the 13:00 snapshot (50.545 kW).

## HC₀ check values (this rebuild, aligned profiles, pf 0.95, limit 253.0 V to neutral)

| Case | Load at 13:00 (kW/customer) | HC₀ (kW/customer) | Binding |
|---|---|---|---|
| Draw 0, 13:00 | 0.464 | **2.999** | B733 phase A; VUF 0.56% |
| Draw 0, worst hour (14:00) | 0.401 | 2.941 | B733 phase A |
| Zero load | 0 | 2.739 | B3278 phase A |
| Draw 0, load ×2 | 0.928 | 3.210 | B733 phase A |
| Draw 0, load ×4 | 1.856 | 3.622 | B733 phase A; VUF 1.13% |
| 101 draws, 13:00 | 0.26–0.55 | min 2.695 / p10 2.891 / **median 3.012** / p90 3.232 / max 3.391 | 4 different buses |

**HC₀(t), draw 0, clear-sky irradiance:**

| Hour | 08 | 09 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|---|---|---|
| HC₀ | 9.49 | 5.66 | 3.64 | 3.40 | 3.07 | 3.00 | **2.94** | 3.13 | 4.21 | 6.91 |

The worst hour is 12:00–14:00 depending on the draw. Within that window the differences are about
2%, so a 13:00 snapshot is sound; checking 12–14 per draw is cheap.

Across draws HC₀ spans ±12% around its median, and its correlation with total load is only 0.11.
What matters is where the load sits relative to the binding buses, not how much there is. That is
why hosting capacity has to be reported over draws.

## The OPF script: `lv12_opf.jl`

Same pipeline as `Small_Network_PV.jl`: `variables.jl` → `constraints_unified.jl` →
`add_pv_control!` → `objectives_FIXED.jl`, core files compiled once. Every generic function is
copied verbatim from that script. New code covers only what the 4-bus layout can't express: one PV
unit per customer, loads from the SGSC files, and the staged checks.

**Where to put it.** Put the whole `lv12_sgsc/` folder next to `Small_Network_PV.jl`. The script
then finds `core/` and your `Project.toml` exactly as that script does. Run:

```
julia --project=<thesis root> lv12_sgsc/lv12_opf.jl
```

**Stages.** Each prints PASS/FAIL and the run stops at the first failure:

| Stage | Does | Must match (OpenDSS, `reference_rebuild/Master_pmd_opendss_check.dss`) |
|---|---|---|
| 1 | Parse `Master_pmd.dss` | 248 buses, 247 four-wire lines, no 5th conductor, 121 loads `[p, 4]`, Σ 56.829 kW / 18.679 kvar, b908 neutral grounded |
| 2 | MEN bonds via `add_men!`; source impedance | 121 bonds = load buses = `men_buses.csv`; `set_loads!` self-test; source self-impedance in the math model = 0.0202 + j0.0735 mΩ |
| 3 | Power flow, no PV | Source 57.304 kW / 18.842 kvar; V_pn 238.85–245.07 V; NEV 2.211 V; VUF 0.318%; every bus within 0.1 V |
| 4 | 121 × 3.0 kW PV, unity PF | Source −300.41 kW / 21.40 kvar; V_pn 243.30–252.97 V; NEV 3.122 V; VUF 0.582% |
| 5 | HC₀ (uniform, unity PF) | **3.025 kVA/customer**. Same 253.05 V pass test as the script; 3.010 at a strict 253.0 V |
| 6 | HC₁ (volt-var + volt-watt) | New result; droop error checked |
| 7 | HC₂ per STATCOM mode, 30 kVA at b733 | New result |

Results go to `results/`: `lv12_stage_checks.csv`, `lv12_hc.csv` and `lv12_points.csv`.

**When a check fails, the pattern says where to look:**
- All voltages offset by the same amount → source voltage or voltage base.
- NEV wrong but phase voltages right → earthing (MEN shunts, or the source neutral).
- Errors growing with distance → linecode units (kft) or line lengths.
- A total-load error in stage 1 → load units in the ENGINEERING model.
- Solves at zero load but LOCALLY_INFEASIBLE under load → too much impedance in series with the
  loads; check the source first (stage 2 prints it). `lv12_diag.jl` separates the cases.

**Decisions built in:**
- **Thermal limits off.** The dataset has no conductor ratings; the 400 A a parser assigns is a
  default. At HC₀ the head cable (400 mm² Al triplex) carries about 422 A, so a thermal check would
  bind first on an invented number. Currents are still reported. Set `CHECK_THERMAL = true` once
  real ratings exist.
- **Stiff source** (OpenDSS default impedance, 0.08 mΩ). The neutral is solidly earthed with the
  `4_Bus_Small.dss` near-zero reactor, which PMD turns into a grounded flag. The impedance is
  written on the circuit line of both PMD masters **and** set in Julia (`SRC_Z1`, `SRC_Z0`,
  `set_source_impedance!`), because PMD's parser, given no impedance, computes OpenDSS's default
  fault level at 115 kV: 1.60 + j6.42 Ω. See "Source impedance" below.
- **MEN 10 Ω per customer**, one customer per bus, so no paralleling.
- **One PV unit per customer**, uniform size, p_avail = S. Count-based HC with placement draws is the
  next step, once stages 1–5 pass.

First run on Patrick's machine (13 Oct): stages 1–2 passed; stage 3 failed on the source impedance
above, now fixed.

## `reference_rebuild/`

The cleaned LV12_248, rebuilt from the raw dataset by `rebuild_clean_lv12.py` following the
cleaning log. Every number in this README was computed on it.

| File | Used by | What it is |
|---|---|---|
| `Master_pmd.dss` | `lv12_opf.jl` | The PMD master: draw-45 loads, source neutral earthed, **no customer MEN reactors** (`add_men!` adds those) |
| `Master_pmd_opendss_check.dss` | OpenDSS | The same plus `Groundings.dss`. The reference numbers come from this |
| `Master.dss` | OpenDSS, `hc0_oracle.py` | The cleaning-log version (ideal source, 1 kW placeholder loads) |
| `Linecodes.dss`, `Lines.dss`, `Switches.dss` | all | Cables and topology. The 7 ex-switches are 1 m lines; without them every customer is cut off |
| `Loads.dss` | `Master.dss` | Placeholder: 1 kW, pf 0.88 per customer (the cleaning-log case) |
| `Loads_h13_draw045.dss` | `Master_pmd*.dss` | Real loads: SGSC draw 45, 13:00, pf 0.95 |
| `Groundings.dss` | OpenDSS only | The 121 MEN bonds as reactors (PMD can't take this form) |
| `men_buses.csv` | stage 2 cross-check | The same 121 buses |
| `opendss_reference_h13_draw045.csv` | stages 3–4 | Per-bus V_pn, NEV and VUF, cases `no_pv` and `pv_3kW_unity` |

## Source impedance: the 2–3% gap, explained

Patrick's earlier cleaned folder used the 4-bus file's `ISC3=100000 ISC1=100000` source (about
2.3 mΩ). Repeating that reproduces his logged OpenDSS runs exactly: HC₀ 2.670 at zero load
(logged 2.67) and max NEV 1.312 V (logged 1.312). The source model matters for HC₀ (zero load):

| Source | HC₀ |
|---|---|
| Stiff (this folder) | 2.74 |
| ISC 100 kA | 2.67 |
| The dataset's transformer resistance kept (3.0 mΩ) | 2.39 |

**PMD parser trap (found 13 Oct).** With no impedance on the `new circuit` line, PMD builds the
source from OpenDSS's default fault level (MVAsc3 2000, MVAsc1 2100) evaluated at OpenDSS's default
115 kV, not the 0.4 kV given: 1.60 + j6.42 Ω instead of 0.019 + j0.078 mΩ. At zero load this only
lifts every bus by 0.14 V (cable charging current through 6.4 Ω), which is how it was spotted:
PMD 244.03–244.07 V against OpenDSS 243.90 V. Under 57 kW the feeder has no solution (OpenDSS
with the same impedance collapses to 71–91 V), so stage 3 came back LOCALLY_INFEASIBLE from PMD's
own model as well as from the core files. The 4-bus never hit this because it states `ISC3`/`ISC1`.
The fix states the impedance in the DSS and sets it in Julia, and stage 2 checks the value that
reaches the math model. Switching to the ISC 100 kA or transformer source means changing `SRC_Z1`
/`SRC_Z0` and the circuit line together, then regenerating the OpenDSS targets.

Reverse power flow through any source impedance raises the LV busbar. Which source is right is a
modelling decision (ask Fred). Note that the dataset's transformer reactance, 0.039%, looks far
too small for a 500 kVA unit; typical values are about 4%.
