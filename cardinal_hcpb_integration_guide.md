# Using your HCPB OpenMC model in Cardinal (OpenMC + nekRS + MOOSE)

This guide shows how to move from your current standalone OpenMC workflow to a Cardinal workflow.

---

## 1) What to keep from your current script

From your `hcpb_temperature_sweep.py`, the parts you should preserve are:

- geometry builders:
  - `make_pebble_unit_cell(...)`
  - `fill_box_with_lattice(...)`
  - `make_first_wall_cooling_pipes_universe(...)`
  - `build_hcpb_module_gap(...)`
- material builder:
  - `build_materials(...)`
- source and tally logic:
  - `run_case_and_get_tbr(...)`

In Cardinal, these should become a reusable OpenMC model factory (a Python module) that can export XML or return an OpenMC model object.

---

## 2) Refactor your model into a reusable OpenMC factory

Create a module, for example `hcpb_openmc_model.py`, with a function like:

```python
def make_hcpb_openmc_model(temp_c: float, geom, gap_front: float,
                           breeder_csv: str, be_csv: str,
                           batches: int = 50, particles: int = 200000):
    materials, mats, temp_k, density_meta = build_materials(temp_c)
    breeder_uc, breeder_pitch = make_pebble_unit_cell(breeder_csv, mats["breeder"], mats["he_purge_mix"])
    be_uc, be_pitch = make_pebble_unit_cell(be_csv, mats["be"], mats["he_purge_mix"])
    geometry, surfs = build_hcpb_module_gap(geom, breeder_uc, breeder_pitch, be_uc, be_pitch, mats, gap_front)

    model = openmc.Model(materials=materials, geometry=geometry, settings=openmc.Settings())
    model.settings.run_mode = "fixed source"
    model.settings.batches = batches
    model.settings.inactive = 0
    model.settings.particles = particles
    model.settings.temperature = {"method": "interpolation", "default": temp_k}

    # add source and tallies exactly as in your current script
    return model
```

Why: Cardinal expects your neutronics setup to be stable and callable each coupled step.

---

## 3) Cardinal-side architecture (conceptual)

Cardinal couples:
- OpenMC for neutronics
- NekRS for CFD (fluid temperature/velocity)
- MOOSE for heat conduction/solid mechanics and coupling orchestration

Typical multiphysics loop:
1. MOOSE/NekRS provide temperature field(s)
2. OpenMC updates material temperatures/densities
3. OpenMC computes heating / reaction rates / TBR
4. Power/heating fields are transferred back to MOOSE/NekRS
5. Iterate until converged (or march transiently)

---

## 4) Mapping your HCPB data into Cardinal coupling quantities

At minimum define consistent fields and region IDs:

- OpenMC side:
  - cell/material IDs for breeder, Be beds, first-wall channels, armour, etc.
  - tallies for:
    - `heating` / `heating-local` (power deposition)
    - `H3-production` (TBR)

- Thermal/CFD side:
  - temperature field(s) on fluid and solid domains
  - optional density updates from your CTE model for solids

Recommendation:
- Keep your current per-material CTE / density law in one Python module.
- Drive it from the temperature fields Cardinal passes into OpenMC each iteration.

---

## 5) Practical migration path (lowest risk)

### Phase A — OpenMC-only in Cardinal

- Use Cardinal to run OpenMC with your HCPB model unchanged.
- Verify:
  - source definition
  - geometry cell naming/IDs
  - TBR tally equals your standalone result within statistics.

### Phase B — one-way coupling (thermal -> neutronics)

- Pass a prescribed temperature field from MOOSE/NekRS to OpenMC.
- Update material temperatures (and densities via your CTE model) each iteration.
- Recompute TBR and heating.

### Phase C — two-way coupling

- Return OpenMC heating to MOOSE/NekRS.
- Solve temperature, feed updated temperature back to OpenMC.
- Iterate until residual criteria are met.

---

## 6) Where your current script needs adaptation for Cardinal

1. **Remove hard-coded single-run loop assumptions**
   - Cardinal controls solve iterations/time steps.

2. **Avoid deleting statepoints globally each run**
   - Current script removes `statepoint.*.h5`; in coupled workflows, let Cardinal handle outputs.

3. **Separate model construction from execution**
   - `make_model(...)` should not call `model.run()` directly.

4. **Keep tallies minimal for coupled runs**
   - heavy mesh tallies/plots can be too expensive each coupling step.

5. **Keep postprocessing external**
   - your standalone plotting scripts should read outputs after coupled execution.

---

## 7) Suggested directory layout

```text
.
├─ hcpb_openmc_model.py          # model builder (geometry/materials/source/tallies)
├─ hcpb_temperature_sweep.py     # standalone sensitivity runs (already present)
├─ plot_tbr_vs_temperature.py    # standalone plotting (already present)
├─ cardinal/
│  ├─ hcpb_cardinal.i            # main Cardinal input
│  ├─ materials/                 # optional material property tables
│  └─ coupling/                  # transfer configs / mappings
└─ cardinal_hcpb_integration_guide.md
```

---

## 8) Verification checklist before full coupling

- Neutronics-only Cardinal result reproduces standalone OpenMC:
  - TBR within 1–2 sigma
  - heating totals consistent
- Temperature feedback test:
  - monotonic density change with temperature from your CTE law
  - expected TBR trend vs temperature
- Energy consistency:
  - deposited heating accounted for in thermal solve

---

## 9) Common pitfalls

- Mismatch between OpenMC cell IDs and Cardinal transfer blocks.
- Overly fine mesh tallies each iteration (slow coupling).
- Using CTE constants outside valid temperature range.
- Inconsistent reference temperature for density law vs material data.

---

## 10) Next concrete step

Implement `hcpb_openmc_model.py` first (pure refactor from your existing file):
- no physics changes,
- only reorganize into importable builder functions.

Then hook that module into a first Cardinal input for OpenMC-only execution and compare TBR against the standalone script.

---

## 11) Phase B code in this repository

Two helper files are now provided:

- `phase_b_temperature_feedback.py`
  - `apply_material_temperature_feedback(...)`
  - `run_phase_b_iteration(...)`
- `run_phase_b_example.py`
  - minimal one-way coupling example with prescribed region temperatures

### How to run the example

```bash
python run_phase_b_example.py
```

This performs one OpenMC solve with temperature feedback applied to solid
materials and reports TBR and total heating.
