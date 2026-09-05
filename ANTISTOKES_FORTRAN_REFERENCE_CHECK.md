# Anti-Stokes feature check against the full-production Fortran reference

**What this is**: a data-analysis check (not a code change) run in a separate chat session,
independent of the Mode B debugging/fix work. It directly inspects the paper-1 Fortran reference
code's own **full production-scale** output (not the shrunk toy validation used elsewhere in this
project) for the anti-Stokes feature the whole project is ultimately targeting. Flagging it here in
case it's useful once Mode B is validated at production scale — it identifies a specific, large,
previously-unnoted feature of the reference data itself that's worth checking the corrected code
against.

## Where the data lives

Real, full-production Fortran runs (**not** `validation_modeB/fortran_toy_reference*/` — those are
deliberately shrunk `nb=nt=1` toy copies for fast local iteration). These are genuine
`nb=4, nt=8, timesteps=320` production runs, i.e. the same scale validated for occupation dynamics
in `validation_modeB/occupation_d8_prxfig2/`:

- **XY polarization**: `~/Documents/GradSchool/Devereaux Group/Herr_Jost/tr_spectroscopy/Ni-L-RIXS/new_shape_function/XY/trSPEC_XY_{0010..1810}.dat` (19 files, already unzipped, ready to use).
- **XX polarization**: same directory's `XX/16_points_new_shape_function_XX.zip` — unzip it to get
  the same 19-file structure (`trSPEC_XX_{0010..1810}.dat`).
- **YY polarization**: `XX/new_shape_function_3_pts_XX.tar.gz`-style archives exist under `YY/` too,
  not yet checked.

Both XX and XY confirmed to share the identical `TD_Cluster_Xray.f90` config: `timesteps=320,
nb=4, nt=8`, `allowed_tm = [-320,-52,0,32,140,160,200,284,362,446,536,626,710,788,878,968,1052,
1136,1214]` (19 delays). Delay (eV⁻¹) = `tm * 2*pi/320`; file tlabel = `100*count + 10` for
`count = 0..18` in the same order as `allowed_tm` (i.e. `trSPEC_XY_0010.dat` = first/most-negative
delay, `_1810.dat` = last/most-positive) — same convention as everywhere else in this project.
Each file is a flat `(incident_energy, energy_loss, intensity)` table, 256 incident-energy points
× 128 loss points, loss ranging from **−1.9375 to +6.0 eV** (so the anti-Stokes/energy-gain side is
present in this data, unlike some of the toy validation runs).

Analysis scripts/plots for this check are saved in `reference_papers/anti_stokes_check/` (not
committed as of this writing — copy them into the repo if useful).

## Finding 1: the anti-Stokes feature is real and pump-induced, confirmed directly in the legacy code's own production output

At the main dd-excitation resonance (incident energy `ab ≈ -4.99 eV`, the same resonance tracked
throughout this project's Mode B debugging), the energy-loss profile:

- **Pre-pump** (`Δt = -6.3 eV⁻¹`, before the pump acts): smooth, featureless, monotonically
  decaying tail on the anti-Stokes (loss<0) side. No structure anywhere.
- **Long post-pump delay** (`Δt ≈ 19-21 eV⁻¹`): a **distinct peak** appears at
  **loss ≈ -0.8 to -0.9 eV**, sitting measurably above where the pre-pump curve is at the same loss
  value, with the two curves crossing right around `-0.85` to `-0.9 eV` — i.e. the pump specifically
  redistributes spectral weight past that point. This is present in **both XX and XY** polarization,
  slightly cleaner/more peak-like in XX, consistent with paper 1's own statement that this feature
  is more distinct in XX. The crossover location matches paper 1's stated anti-Stokes maximum at
  `-0.9 eV` closely.

This is an independent, from-the-reference-code confirmation that the target physics is real and
present in the Fortran implementation the whole project has been validating C++ Mode B against —
useful as a sanity check that the target actually exists in the data, separate from the question of
whether the corrected C++ implementation reproduces it.

## Finding 2 (the one worth specifically flagging): the anti-Stokes signal is far more intense *off* the main resonance

This is the new, more actionable result. Scanning the anti-Stokes peak intensity (max intensity for
loss in `[-1.0,-0.75]` eV) across the **entire incident-energy axis** at long post-pump delay reveals:

- The dominant anti-Stokes resonance sits at **`ab ≈ -7.3 eV`** — about **2.3 eV below** the main
  dd-excitation resonance at `-4.99 eV` — with peak intensity **~10-30× larger** than the
  anti-Stokes signal measured right at the main resonance (both XX and XY show this; XX:
  `8.2e-6` at `ab=-7.3` vs. `2.8e-7` at `ab=-4.99`; XY: `2.0e-6` at `ab=-7.3` vs. `1.1e-7` at
  `ab=-4.99`).
- This `-7.3 eV` peak is a **genuine local maximum**, not an edge-of-grid artifact — intensity rises
  from the grid boundary (`ab=-8.70`), peaks sharply at `-7.3`, then falls back down through
  `-6.5`, `-5.8` before reaching the much smaller main-resonance signal at `-4.99`. It is present at
  both long delays checked (`Δt=19.0` and `20.7`) and essentially absent pre-pump (pre-pump
  intensity at the same `ab` is roughly 10× smaller) — i.e. it's pump-induced, not a pre-existing
  equilibrium feature that happens to be large.
- Two smaller but still clearly pump-enhanced secondary features also appear around
  `ab ≈ +10 eV` and `ab ≈ +13.5 eV` in both polarizations.

**Why this matters for the Mode B work**: every numerical comparison against the Fortran reference
run in this project so far (the toy-validation debugging chased through `SESSION_HANDOFF.md`,
`THEORY_CROSSCHECK_PRB99.md`, etc.) has focused on the **main resonance** (`ab≈-4.99 eV`,
`loss≈2 eV` for the dd-excitation feature specifically). This analysis shows that, at least in the
real production-scale reference data, the anti-Stokes signal is **dominated by a completely
different, unexamined incident-energy region**. Once Mode B is confirmed correct at the main
resonance (per the fix in `THEORY_CROSSCHECK_PRB99.md`), it would be worth specifically checking
whether the corrected C++ implementation also reproduces this **off-resonance** enhancement at
`ab≈-7.3 eV` — this has not been checked at all yet, in either the buggy or the fixed version of
Mode B, and could plausibly behave differently from the main-resonance feature (different
intermediate core-hole states are presumably involved at a resonance that far detuned from the main
peak).

## Suggested next step (data-analysis, not code)

Once Mode B produces production-scale output, rerun the same "anti-Stokes intensity vs. incident
energy" scan (the method above: for each incident-energy grid point, take the max intensity in a
loss window around `[-1.0,-0.75] eV`, plot vs. incident energy, compare pre-pump vs. long-delay
post-pump) on the C++ output and check for the same `ab≈-7.3 eV` dominant peak. If it's absent or
much weaker in the corrected C++ output than in this Fortran reference, that would be a new,
specific, well-localized discrepancy to chase — distinct from (and possibly unrelated to) the
main-resonance issue this project has spent the most effort on.
