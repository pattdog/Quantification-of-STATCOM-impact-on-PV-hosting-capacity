# Quantification of STATCOM impact on PV hosting capacity

This repository holds the code and results behind my thesis. It computes the rooftop-PV hosting capacity of four-wire low-voltage feeders with MEN earthing, with and without a STATCOM, as an unbalanced optimal power flow in PowerModelsDistribution solved with Ipopt. Every PV inverter follows the AS/NZS 4777.2 volt-var and volt-watt response, and the STATCOM can exchange reactive power, real power between phases, or both. Two networks are studied: a four-bus test feeder (Chapter 5) and the LV12_248 feeder from the Geth et al. (2025) dataset with measured Smart Grid Smart City household loads (Chapter 6).

## Where things are

All the code written for the thesis is under `rosetta_distribution_opf.jl/script_inverter_PSCC2024/`. That folder sits inside the `rosetta_distribution_opf.jl` package because the scripts find the shared model files, the network data and the Julia project by searching upward from their own location, so the layout should not be changed. `core/` holds the network, inverter, STATCOM and limit model that every final script loads (`variables.jl`, `constraints_unified.jl` and `objectives_FIXED.jl`); the other files in `core/` are earlier versions kept for the record. `ThesisContribution/` holds the study scripts. The rest of `rosetta_distribution_opf.jl` is the package of that name by Rahmat Heidari, which every script loads with `using rosetta_distribution_opf`, and its `data/` folder holds the network files the scripts read.

The table maps each part of the thesis to the script that produces it. Script paths are relative to `ThesisContribution/`. Results paths are relative to the repository root, except for LV12_248, whose results sit next to its script. The plotting scripts read the CSV files and never solve.

| Thesis part | Script | Results | Plots from the CSVs |
|---|---|---|---|
| Chapter 5, load-only study on the four-bus feeder | `SmallNetworkTest/No_PV/Small_Network_Sweep2.jl` | `results/small_network_2/` | `SmallNetworkTest/No_PV/Figures2.jl` |
| Chapter 5, PV hosting capacity on the four-bus feeder | `SmallNetworkTest/PV/SmallNetworkPV_RUNTIME_SHRINK.jl` | `results/small_network_pv/` | `SmallNetworkTest/PV/Figures.jl` |
| Chapter 6, LV12_248 case study | `LargeCaseStudy/lv12_sgsc/lv12_opf.jl` | `LargeCaseStudy/lv12_sgsc/results/_cableEQ/` | — |
| Three-bus preliminary study | `SmallNetworkTest/No_PV/Small_Network_Sweep.jl` | `results/toy_unbalance/` | `SmallNetworkTest/No_PV/figures.jl` |
| Early exploratory scripts on the ENWL feeder (September 2026) | the `.jl` files directly in `ThesisContribution/` | `plots/` | — |

The Chapter 6 results are in `results/_cableEQ/`, all on the Energex STNW3369 cable ratings: `n5x10/` holds the main set at 5 existing fleets × 10 applicant orders, and the `adopt_*` folders directly under `_cableEQ/` hold the circuit-A siting runs at b733 at 5 × 5. Sensitivity runs write to folders named after their settings (for example `pen80/` for 80% existing penetration). The `adopt_*` folders, `nothermal/` and the CSV files directly under `lv12_sgsc/results/` are runs from before 26 October 2026 at a uniform 400 A rating (or with thermal limits off), which the thesis does not use.

## Setting up

The scripts are written for Julia 1.11. After cloning, register the bundled package with the project and install the dependencies once:

```
git clone https://github.com/pattdog/Quantification-of-STATCOM-impact-on-PV-hosting-capacity.git
cd Quantification-of-STATCOM-impact-on-PV-hosting-capacity
julia --project=. -e 'using Pkg; Pkg.develop(path="rosetta_distribution_opf.jl"); Pkg.instantiate()'
```

Ipopt is installed by Ipopt.jl with its default linear solver. The Python helpers for LV12_248 need Python 3, and `hc0_oracle.py` also needs OpenDSSDirect.py (`pip install opendssdirect.py`).

## Running the studies

Run every command from the repository root. Each script activates the project itself, so the `--project` flag is optional. The four-bus studies are deterministic: they use no random numbers, and all settings are constants at the top of each script.

```
S=rosetta_distribution_opf.jl/script_inverter_PSCC2024/ThesisContribution
julia --project=. $S/SmallNetworkTest/No_PV/Small_Network_Sweep2.jl
julia --project=. $S/SmallNetworkTest/PV/SmallNetworkPV_RUNTIME_SHRINK.jl
```

`lv12_opf.jl` runs in stages. Stages 1 to 5 check the PowerModelsDistribution model against OpenDSS (network size, earthing, power flow with and without PV, and HC₀) and stop at the first failure, so nothing new is computed until the port is confirmed. Stage 11 is the household count used in Chapter 6: 55% of customers start with a 5 kVA inverter, the rest apply one at a time in random order, and hosting capacity is the number of applicants who connect without breaching a limit. Fleets and orders are drawn with fixed seeds (base seed 2026, one per fleet and per order), so a run repeats exactly, and the first five orders of each fleet in the 5 × 10 set reproduce the 5 × 5 runs. Each configuration is run separately and writes to its own folder:

```
L=$S/LargeCaseStudy/lv12_sgsc
LV12_NINIT=5 LV12_NNEW=10 LV12_CONFIGS=HC1 julia --project=. $L/lv12_opf.jl
```

The main set uses `HC1`, `HC2_D_one_coord`, `HC2_D_one_vbal`, `HC2_D_pair_coord`, `HC2_D_pair_vbal`, `HC2_C_one_coord`, `HC2_C_one_vbal` and `HC2_B_one_coord`, which write to `results/_cableEQ/n5x10/adopt_<configuration>/`. The b733 siting runs use `HC2_D_pair733_coord`, `HC2_D_pair733_vbal`, `HC2_D_oneA733_coord` and `HC2_D_oneA733_vbal` with `LV12_NINIT=5 LV12_NNEW=5`, and write to `results/_cableEQ/adopt_<configuration>/`. In the names, D, C and B are the STATCOM modes (reactive and real power, real-power exchange only, reactive only), `one` is a single 30 kVA unit (at b3278, or at b733 in `oneA733`) and `pair` is two (b3278 with b3162, or with b733 in `pair733`), and `coord` and `vbal` are the coordinated and local control tiers. The sensitivities (cables unrated, other penetrations, transformer rating and source voltage) are set with the environment variables documented at the top of `lv12_opf.jl`. The overnight runs were launched with `run_queue.sh <queue file> <N>`, which runs one configuration per line of the queue file, N at a time, and writes the logs to `logs/`.

The prepared LV12_248 inputs are committed, so `lv12_opf.jl` runs without rebuilding anything: `reference_rebuild/` holds the cleaned network and the OpenDSS reference values, and `out/` holds the load profiles and their assignment to customers. The scripts that made them are included. Run from `lv12_sgsc/`, `reference_rebuild/rebuild_clean_lv12.py <raw LV12_248 folder> reference_rebuild` rebuilds the cleaned network from the dataset's raw files; `make_lv12_loads.py smartgridsmartcities.json reference_rebuild/Loads.dss out` assigns the measured profiles (draw 0 cyclic, draws 1 to 100 random with seed 2026); and `hc0_oracle.py reference_rebuild/Master.dss out --draws 0-100 --hours 13` finds HC₀ in OpenDSS for each draw. `lv12_sgsc/README.md` describes the load method in full.

## Data and credits

The `rosetta_distribution_opf` package (`rosetta_distribution_opf.jl/src/`) is by Rahmat Heidari and is included as it was imported. The models are built on PowerModelsDistribution and JuMP and solved with Ipopt. The LV12_248 network comes from the Geth et al. (2025) CSIRO dataset of Australian low-voltage feeders (doi:10.25919/ghnz-bk28). The household loads are the Smart Grid Smart City profiles in `smartgridsmartcities.json` from CSIRO's RepresentativeLVNetworks repository, which is not redistributed here. Cable ratings are from Energex STNW3369, Release 5 (12 June 2026). The four-bus network's parameters and their sources are documented in the header of `rosetta_distribution_opf.jl/data/4_Bus_Small.dss`.
