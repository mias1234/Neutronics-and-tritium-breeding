# Convert `model.xml` to a CAD-compatible file

For OpenMC CSG models, direct STEP/IGES export is generally not available.
A practical workflow is to export a **voxel mesh** (`.vti` or `.vtk`) and optionally a
surface mesh (`.stl`) for CAD/visualization tools.

This repository includes `model_xml_to_cad.py` for that.

## 1) Export your OpenMC XML first

From Python:

```python
model.export_to_model_xml()
```

This writes `model.xml`.

## 2) Convert to VTI (default)

```bash
python model_xml_to_cad.py --model-xml model.xml --out-stem hcpb_mesh
```

Outputs:
- `hcpb_mesh.h5` (OpenMC voxel plot output)
- `hcpb_mesh.vti` (VTK image data)

## 3) Convert to legacy VTK (`.vtk`)

```bash
python model_xml_to_cad.py --model-xml model.xml --out-stem hcpb_mesh --vtk-format vtk
```

Outputs:
- `hcpb_mesh.vtk`

> Note: `openmc.voxel_to_vtk` produces `.vti` natively. The script now converts `.vti` to legacy `.vtk` when `--vtk-format vtk` is requested (requires `pyvista`).

## 4) Optional STL output

If `pyvista` is available (also needed for `--vtk-format vtk`):

```bash
python model_xml_to_cad.py --model-xml model.xml --out-stem hcpb_mesh --stl
```

Output:
- `hcpb_mesh.stl`

## 5) Resolution control

Increase mesh resolution if geometry appears coarse:

```bash
python model_xml_to_cad.py --model-xml model.xml --nx 400 --ny 240 --nz 240
```

Higher resolution increases runtime and file size.

## Notes

- This is a **mesh approximation** of CSG geometry (not exact boundary-rep CAD).
- If your geometry bounding box is unbounded, conversion will fail until you add
  finite boundaries in the model.


## 6) Export named volumes/surfaces for ParaView identification

To generate identifiable volume and boundary datasets:

```bash
python model_xml_to_cad.py --model-xml model.xml --out-stem hcpb_mesh --vtk-format vtk --export-named-sets
```

This writes `hcpb_mesh_named/` with:
- volume files (when present):
  - `volume_breeder.vtu`
  - `volume_be.vtu`
  - `volume_steel.vtu`
  - `volume_armour.vtu`
  - `volume_coolant.vtu`
- boundary files:
  - `boundary_cooling_surface.vtp`
  - `boundary_vacuum_surface.vtp`
  - symmetry planes (if detected):
    - `boundary_symmetry_ymin.vtp`
    - `boundary_symmetry_ymax.vtp`
    - `boundary_symmetry_zmin.vtp`
    - `boundary_symmetry_zmax.vtp`
- `named_sets_manifest.json`

In ParaView, open files in `hcpb_mesh_named/` to inspect each volume/surface by name.


## 7) Embed region tags directly in the `.vtk/.vti`

If you want region IDs available **inside** the main mesh file:

```bash
python model_xml_to_cad.py --model-xml model.xml --out-stem hcpb_mesh --vtk-format vtk --embed-volume-tags
```

This adds cell-data arrays to `hcpb_mesh.vtk`:
- `openmc_material_id`
- `volume_tag_id`

Default `volume_tag_id` mapping:
- 1: breeder
- 2: be
- 3: steel
- 4: armour
- 5: coolant

Mapping is written to `hcpb_mesh_named/volume_tag_manifest.json`.

In ParaView:
1. Open `hcpb_mesh.vtk`
2. In the coloring dropdown, choose `volume_tag_id`
3. Use `Threshold` filter to isolate each volume tag.
