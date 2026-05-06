using OpenDSSDirect
const ODSS = OpenDSSDirect

# 1. Define a tiny 2-bus network directly in the script
dss_script = """
clear
new circuit.SimpleTest basekv=0.4 phases=3 bus1=sourcebus

! Define a short line (100 meters)
new linecode.abc nphases=3 rmatrix=[0.25 | 0.02 0.25 | 0.02 0.02 0.25] xmatrix=[0.05 | 0.01 0.05 | 0.01 0.01 0.05] units=km
new line.line1 bus1=sourcebus bus2=loadbus linecode=abc length=0.1 units=km

! Add a heavy unbalanced load (Phase A is much heavier)
new load.load1 phases=1 bus1=loadbus.1 kv=0.23 kw=50 kvar=20
new load.load2 phases=1 bus1=loadbus.2 kv=0.23 kw=10 kvar=5
new load.load3 phases=1 bus1=loadbus.3 kv=0.23 kw=10 kvar=5

solve
"""

ODSS.dss(dss_script)

# 2. Check Baseline Voltages
println("--- Baseline (No STATCOM) ---")
ODSS.Circuit.SetActiveBus("loadbus")
magnitudes = abs.(ODSS.Bus.Voltages())
println("Phase Magnitudes (V): ", round.(magnitudes[1:3], digits=2))
# 3. Inject STATCOM support on Phase A (Reactive power only)
# Model=7 in OpenDSS acts like a constant kVAR source
statcom_script = """
new generator.my_statcom.1 phases=1 bus1=loadbus.1 kv=0.23 kw=0 kvar=400 model=7
new generator.my_statcom.2 phases=1 bus1=loadbus.2 kv=0.23 kw=0 kvar=100 model=7
new generator.my_statcom.3 phases=1 bus1=loadbus.3 kv=0.23 kw=0 kvar=400 model=7


solve
"""

ODSS.dss(statcom_script)

# 4. Check Improved Voltages
println("\n--- With STATCOM Support ---")
ODSS.Circuit.SetActiveBus("loadbus")
magnitudes = abs.(ODSS.Bus.Voltages())
println("Phase Magnitudes (V): ", round.(magnitudes[1:3], digits=2))