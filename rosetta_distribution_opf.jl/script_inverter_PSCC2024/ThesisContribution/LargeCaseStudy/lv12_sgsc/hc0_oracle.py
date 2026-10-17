"""
HC0 oracle for LV12_248 in OpenDSS: the largest uniform PV size per customer (unity PF, no
volt-var/volt-watt, PV on all 121 customers on their own phase) before any bus exceeds the
upper phase-to-neutral limit. Found by bisection. This is the number the OPF's HC0 must reproduce
under the same loads, source and earthing.

Usage:
  python3 hc0_oracle.py MASTER.dss OUTDIR_OF_make_lv12_loads [--draws 0] [--hours 13]
                        [--pf 0.95] [--vmax 253.0] [--load-scale 1.0] [--irr IRRADIANCE_CSV]
  --draws accepts e.g. 0  or  0-20 ;  --hours accepts e.g. 13  or  9-17
"""
import argparse, csv, math, cmath
import opendssdirect as dss

ap = argparse.ArgumentParser()
ap.add_argument('master'); ap.add_argument('loaddir')
ap.add_argument('--draws', default='0'); ap.add_argument('--hours', default='13')
ap.add_argument('--pf', type=float, default=0.95); ap.add_argument('--vmax', type=float, default=253.0)
ap.add_argument('--load-scale', type=float, default=1.0)
ap.add_argument('--ref', default='neutral', choices=['neutral', 'ground'],
                help='voltage limit referenced to the neutral (customer sees this; default) or to ground')
ap.add_argument('--irr', default=None, help='24-row irradiance csv (midnight-indexed); default 1.0')
a = ap.parse_args()
rng = lambda s: list(range(int(s.split('-')[0]), int(s.split('-')[-1]) + 1))
A = cmath.exp(2j * math.pi / 3)

prof = {r['profile_id']: [float(r[f'h{h:02d}']) for h in range(24)]
        for r in csv.DictReader(open(f'{a.loaddir}/sgsc_profiles_hourly.csv'))}
assign = {}
for r in csv.DictReader(open(f'{a.loaddir}/lv12_assignment.csv')):
    assign.setdefault(int(r['draw']), []).append((r['load'], r['bus'], r['phase'], r['profile_id']))
irr = [float(x) for x in open(a.irr).read().split()] if a.irr else [1.0] * 24
tanphi = math.tan(math.acos(a.pf))

dss.Text.Command(f'compile "{a.master}"')
loads = assign[0]
for n, b, ph, _ in loads:   # one PV generator per customer, same bus and phase
    dss.Text.Command(f'new generator.pv_{n} phases=1 bus1={b}.{"ABC".index(ph)+1}.4 kv=0.230940 '
                     f'kw=0 pf=1 model=1 vminpu=0.1 vmaxpu=2')

def set_case(draw, hour, size):
    for n, b, ph, k in assign[draw]:
        kw = prof[k][hour] * a.load_scale
        dss.Loads.Name(n); dss.Loads.kW(kw); dss.Loads.kvar(kw * tanphi)
        dss.Generators.Name(f'pv_{n}'); dss.Generators.kW(size * irr[hour])
    dss.Solution.Solve()
    return dss.Solution.Converged()

def worst():
    vmax, at, vuf = 0.0, None, 0.0
    for bus in dss.Circuit.AllBusNames():
        dss.Circuit.SetActiveBus(bus); nodes = dss.Bus.Nodes(); v = dss.Bus.Voltages()
        ph = {n: complex(v[2*i], v[2*i+1]) for i, n in enumerate(nodes)}
        if 4 not in ph: continue
        pv = {p: ph[p] - ph[4] for p in (1, 2, 3) if p in ph}
        for p, x in pv.items():
            m = abs(x) if a.ref == 'neutral' else abs(ph[p])
            if m > vmax: vmax, at = m, f'{bus}.{p}'
        if len(pv) == 3:
            v1 = (pv[1] + A*pv[2] + A*A*pv[3]) / 3; v2 = (pv[1] + A*A*pv[2] + A*pv[3]) / 3
            vuf = max(vuf, 100 * abs(v2) / abs(v1))
    return vmax, at, vuf

print('draw hour irr  load_kW/cust  HC0_kW/cust  Vmax@HC   at        VUF@HC%')
for d in rng(a.draws):
    for h in rng(a.hours):
        if irr[h] <= 0: print(f'{d:4d} {h:4d} {irr[h]:4.2f}  no PV output at this hour'); continue
        ok = lambda x: set_case(d, h, x) and worst()[0] <= a.vmax   # non-convergence counts as a fail
        lo, hi = 0.0, 15.0
        assert ok(lo), 'limit already violated with no PV'
        if ok(hi): print(f'{d:4d} {h:4d}  HC0 > {hi}'); continue
        while hi - lo > 0.002:
            mid = (lo + hi) / 2
            (lo, hi) = (mid, hi) if ok(mid) else (lo, mid)
        set_case(d, h, lo); vm, at, vuf = worst()
        lk = sum(prof[k][h] for *_, k in assign[d]) / len(assign[d]) * a.load_scale
        print(f'{d:4d} {h:4d} {irr[h]:4.2f}  {lk:12.3f}  {lo:11.3f}  {vm:7.2f}  {at:10s} {vuf:7.3f}')
