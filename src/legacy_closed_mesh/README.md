# Superseded closed-convex-hull driver

`bulge_dent_driver.py` is the original driver that built a **closed** convex-hull
interface mesh (via the external `st3d_interface.compute_interface`). It produced
the results in `results/legacy_closed_mesh/` and is kept only for provenance.
The current analysis uses `src/open_interface.py` (genuine open interface sheet).
