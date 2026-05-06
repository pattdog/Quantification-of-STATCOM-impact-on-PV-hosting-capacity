using OpenDSSDirect
const ODSS = OpenDSSDirect


dss_script = """
clear
new circuit.SolarTest basekv=0.4 phases=3 bus1=sourcebus
new linecode.abc nphases=3 rmatrix=[0.25 | 0.02 0.25 | 0.02 0.02 0.25] xmatrix=[0.05 | 0.01 0.05 | 0.01 0.01 0.05] units=km
new line.line1 bus1=sourcebus bus2=loadbus linecode=abc length=0.1 units=km

new generator.solar_pv phases=1 bus1=loadbus.1 kv=0.23 kw=120 kvar=0 model=7
solve
"""

ODSS.dss(dss_script)
println("--- 1. High Solar PV (The Violation) ---")
ODSS.Circuit.SetActiveBus("loadbus")
v_high = round.(abs.(ODSS.Bus.Voltages()[1:3]), digits=2)
println("Phase Magnitudes (V): ", v_high)
println("Status: Phase A is at ", v_high[1], "V (Legal Limit is 253V!)")


statcom_script = """
new generator.statcom.1 phases=1 bus1=loadbus.1 kv=0.23 kw=0 kvar=-380 model=7
new generator.statcom.2 phases=1 bus1=loadbus.2 kv=0.23 kw=0 kvar=10  model=7
new generator.statcom.3 phases=1 bus1=loadbus.3 kv=0.23 kw=0 kvar=-150  model=7
solve
"""

ODSS.dss(statcom_script)
println("\n--- 2. With STATCOM Mitigation (Absorption) ---")
ODSS.Circuit.SetActiveBus("loadbus")
v_fixed = round.(abs.(ODSS.Bus.Voltages()[1:3]), digits=2)
println("Phase Magnitudes (V): ", v_fixed)