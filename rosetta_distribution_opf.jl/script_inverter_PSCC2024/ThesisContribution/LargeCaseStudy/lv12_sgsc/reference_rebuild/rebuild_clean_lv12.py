"""Rebuild cleaned LV12_248 from the raw Geth 2025 files, following LV12_248_cleaning_log.md."""
import re, os, sys
RAW, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
rd = lambda f: open(os.path.join(RAW, f)).read()
# Linecodes: drop duplicate definitions (identical data), keep first
blocks, seen = [], set()
for blk in re.split(r'(?im)^(?=new\s+linecode\.)', rd('Linecodes.dss')):
    m = re.match(r'(?i)new\s+linecode\.(\S+)', blk.strip())
    if not m: continue
    if m.group(1).lower() in seen: continue
    seen.add(m.group(1).lower()); blocks.append(blk.rstrip() + '\n')
open(os.path.join(OUT, 'Linecodes.dss'), 'w').write(''.join(blocks))
open(os.path.join(OUT, 'Lines.dss'), 'w').write(rd('lines.dss'))
# Switches: closed ones become ordinary 1 m lines; commented-out (open points) dropped
sw = [l.replace(' switch=y', '') for l in rd('Switches.dss').splitlines()
      if l.strip().lower().startswith('new line')]
open(os.path.join(OUT, 'Switches.dss'), 'w').write('\n'.join(sw) + '\n')
# Groundings: keep customer MEN bonds (two-terminal .4 -> .0, 10 ohm); drop the 0.3 ohm tx reactor
g = [l for l in rd('Groundings.dss').splitlines() if l.strip().lower().startswith('new') and 'tx_busgrounding' not in l]
open(os.path.join(OUT, 'Groundings.dss'), 'w').write('\n'.join(g) + '\n')
# Loads: names/buses as raw, explicit power, LN kV. Base case = 1 kW, pf 0.88 (cleaning log)
L = []
for l in rd('Loads.dss').splitlines():
    if not l.strip().lower().startswith('new load'): continue
    name = re.search(r'(?i)load\.(\S+)', l).group(1); bus = re.search(r'(?i)bus1=(\S+)', l).group(1)
    L.append(f'new load.{name} phases=1 bus1={bus} kv=0.230940 kw=1 kvar=0.5397 model=1 vminpu=0.1 vmaxpu=2')
open(os.path.join(OUT, 'Loads.dss'), 'w').write('\n'.join(L) + '\n')
open(os.path.join(OUT, 'Master.dss'), 'w').write("""clear
set defaultbasefrequency=50
new circuit.lv12_248 basekv=0.4 pu=1.05610 phases=3 bus1=B908 model=ideal
redirect Linecodes.dss
redirect Lines.dss
redirect Switches.dss
redirect Loads.dss
redirect Groundings.dss
! solid neutral-earth at the source (cleaning log: replaces the 0.3 ohm reactor)
new reactor.source_earth phases=1 bus1=B908.4 bus2=B908.0 r=0.000001 x=0
set voltagebases=[0.4]
calcvoltagebases
""")
print('rebuilt', OUT, len(L), 'loads')
