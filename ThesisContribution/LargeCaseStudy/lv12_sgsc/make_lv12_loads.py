"""
Assign measured SGSC household load profiles to the 121 customers of LV12_248.

Source: smartgridsmartcities.json from CSIRO's RepresentativeLVNetworks repo (408 Smart Grid
Smart City customers, half-hourly, local time). Each customer is a DIFFERENT 24-hour window
(dates Oct 2011 - Jul 2013, 32 different start times, mostly autumn/winter). Each profile is
therefore aligned to clock time using its own timestamps. Windows that are not 48 consecutive
half-hours (4 customers) are dropped, and so are profiles that are zero for all 48 half-hours
(41 customers: vacant or not recording; keep them with --keep-zero).

Assignment follows CSIRO (src/core/load.jl :: load_matrix_to_dict): load i -> profile i, cyclic.
Load order = Loads.dss order; profile order = JSON file order (CSIRO's Julia code iterates a Dict,
whose order is arbitrary, so we fix it here).
  draw 0  : cyclic over the usable pool, in file order
  draw 1+ : 121 distinct profiles drawn at random without replacement (seeded)

Time convention: timestamps mark the END of each half-hour, so clock hour H (H:00-H+1:00) is the
mean of the half-hours stamped H:30 and H+1:00. Each customer's day is stitched from the two
calendar days its window spans.

Usage:
  python3 make_lv12_loads.py SGSC_JSON LOADS_DSS OUTDIR [--draws 100] [--seed 2026] [--pf 0.95]
                             [--season all|summer|autumn|winter|spring] [--hour 13]
"""
import argparse, csv, json, math, os, random, re
from datetime import datetime, timedelta

ap = argparse.ArgumentParser()
ap.add_argument('sgsc_json'); ap.add_argument('loads_dss'); ap.add_argument('outdir')
ap.add_argument('--draws', type=int, default=100, help='random draws in addition to draw 0')
ap.add_argument('--seed', type=int, default=2026)
ap.add_argument('--pf', type=float, default=0.95, help='lagging power factor applied to every load')
ap.add_argument('--season', default='all', choices=['all', 'summer', 'autumn', 'winter', 'spring'])
ap.add_argument('--hour', type=int, default=13, help='clock hour written to the snapshot Loads file')
ap.add_argument('--keep-zero', action='store_true', help='keep profiles that are zero all day')
ap.add_argument('--snapshot-draws', default='0,45', help='draws written as snapshot Loads files (45 = median HC0)')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
H = [f'h{h:02d}' for h in range(24)]
SEASON = {12: 'summer', 1: 'summer', 2: 'summer', 3: 'autumn', 4: 'autumn', 5: 'autumn',
          6: 'winter', 7: 'winter', 8: 'winter', 9: 'spring', 10: 'spring', 11: 'spring'}

# ---- 1. SGSC profiles, each aligned to clock time by its own timestamps -----------------
raw = json.load(open(a.sgsc_json))                 # Python dicts keep file order
prof, meta, dropped, zero = {}, {}, [], []
for k, v in raw.items():
    t = [datetime.fromisoformat(s) for s in v['time']]
    p = v['p']
    if len(t) != 48 or any(t[i + 1] - t[i] != timedelta(minutes=30) for i in range(47)):
        dropped.append(k); continue
    assert min(p) >= 0, f'negative power in {k}: PV netted in?'
    if max(p) == 0 and not a.keep_zero:
        zero.append(k); continue
    at = {(x.hour, x.minute): pi for x, pi in zip(t, p)}   # each clock half-hour exactly once
    assert len(at) == 48
    prof[k] = [(at[(h, 30)] + at[((h + 1) % 24, 0)]) / 2 for h in range(24)]
    mid = t[0] + (t[-1] - t[0]) / 2
    meta[k] = dict(window_start=t[0].isoformat(timespec='minutes'), season=SEASON[mid.month],
                   day='weekend' if t[0].weekday() >= 5 else 'weekday')
with open(os.path.join(a.outdir, 'sgsc_profiles_hourly.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['profile_id', 'window_start', 'season', 'day'] + H)
    for k in prof:
        w.writerow([k, meta[k]['window_start'], meta[k]['season'], meta[k]['day']]
                   + [f'{x:.4f}' for x in prof[k]])
pool = [k for k in prof if a.season == 'all' or meta[k]['season'] == a.season]

# ---- 2. LV12_248 customers, in Loads.dss order ------------------------------------------
loads = []
for line in open(a.loads_dss):
    m = re.search(r'(?i)^\s*new\s+load\.(\S+).*?bus1=(\w+)\.(\d)\.4', line)
    if m: loads.append((m.group(1), m.group(2), 'ABC'[int(m.group(3)) - 1]))
assert len(loads) == 121, f'expected 121 loads, found {len(loads)}'
assert len(pool) >= len(loads), f'only {len(pool)} {a.season} profiles for {len(loads)} loads'

# ---- 3. Assignment draws ---------------------------------------------------------------
rng = random.Random(a.seed)
draws = {0: [pool[i % len(pool)] for i in range(len(loads))]}
for d in range(1, a.draws + 1):
    draws[d] = rng.sample(pool, len(loads))
with open(os.path.join(a.outdir, 'lv12_assignment.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['draw', 'load', 'bus', 'phase', 'profile_id'])
    for d, pick in draws.items():
        for (n, b, ph), k in zip(loads, pick): w.writerow([d, n, b, ph, k])

# ---- 4. Convenience files for draw 0 ---------------------------------------------------
tanphi = math.tan(math.acos(a.pf))
with open(os.path.join(a.outdir, 'lv12_load_kw_draw000.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['load', 'bus', 'phase', 'profile_id'] + H)
    for (n, b, ph), k in zip(loads, draws[0]): w.writerow([n, b, ph, k] + [f'{x:.4f}' for x in prof[k]])

def load_line(n, b, ph, kw, extra=''):
    return (f'new load.{n} phases=1 bus1={b}.{"ABC".index(ph) + 1}.4 kv=0.230940 '
            f'kw={kw:.4f} kvar={kw * tanphi:.4f} model=1 vminpu=0.1 vmaxpu=2{extra}')

for d in [int(x) for x in a.snapshot_draws.split(',')]:
    with open(os.path.join(a.outdir, f'Loads_h{a.hour:02d}_draw{d:03d}.dss'), 'w') as f:
        f.write(f'! LV12_248 loads, clock hour {a.hour:02d}:00-{a.hour + 1:02d}:00, SGSC draw {d}, pf {a.pf}\n')
        for (n, b, ph), k in zip(loads, draws[d]): f.write(load_line(n, b, ph, prof[k][a.hour]) + '\n')
with open(os.path.join(a.outdir, 'Loadshapes_draw000.dss'), 'w') as f:
    f.write('! 24 hourly shapes in kW and kvar (useactual=yes); hour 0 = 00:00-01:00\n')
    for (n, b, ph), k in zip(loads, draws[0]):
        f.write(f'new loadshape.ls_{n} npts=24 interval=1 useactual=yes '
                f'mult=[{" ".join(f"{x:.4f}" for x in prof[k])}] '
                f'qmult=[{" ".join(f"{x * tanphi:.4f}" for x in prof[k])}]\n')
with open(os.path.join(a.outdir, 'Loads_daily_draw000.dss'), 'w') as f:
    f.write(f'! Redirect Loadshapes_draw000.dss first. pf {a.pf}\n')
    for (n, b, ph), k in zip(loads, draws[0]):
        f.write(load_line(n, b, ph, 1.0, f' daily=ls_{n}') + '\n')

print(f'{len(raw)} profiles read, {len(dropped)} dropped (not 48 consecutive half-hours): {dropped}')
print(f'{len(zero)} all-zero profiles dropped' + ('' if not zero else ' (use --keep-zero to keep)'))
print(f'pool ({a.season}): {len(pool)}; {len(loads)} loads; {len(draws)} draws -> {a.outdir}')
