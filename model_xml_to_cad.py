"""Convert OpenMC XML model geometry into CAD-compatible mesh artifacts.

What this does:
1) Loads OpenMC model geometry from `model.xml` (or geometry/materials XML)
2) Generates a voxelized geometry file with OpenMC plotting
3) Converts voxel H5 -> VTI (`openmc.voxel_to_vtk`)
4) Optionally converts VTI -> STL (if pyvista is installed)

This is a practical mesh-based CAD export path for CSG OpenMC models.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import openmc


def _load_model(model_xml: str | None, geometry_xml: str | None, materials_xml: str | None) -> openmc.Model:
    if model_xml is not None:
        model_path = Path(model_xml)
        if not model_path.exists():
            raise FileNotFoundError(f"model.xml not found: {model_xml}")
        return openmc.Model.from_model_xml(path=model_xml)

    if geometry_xml is None:
        raise ValueError("Provide --model-xml, or provide --geometry-xml (and optionally --materials-xml).")

    geometry = openmc.Geometry.from_xml(path=geometry_xml)
    materials = openmc.Materials.from_xml(path=materials_xml) if materials_xml else None
    return openmc.Model(geometry=geometry, materials=materials)


def _safe_bbox(geometry: openmc.Geometry):
    ll, ur = geometry.bounding_box
    if any(x == float("inf") or x == float("-inf") for x in (*ll, *ur)):
        raise ValueError(
            "Geometry bounding box is not finite. "
            "Please use finite vacuum/transmission boundaries or provide a bounded geometry."
        )
    return ll, ur


def convert_model_to_vti(
    model: openmc.Model,
    out_stem: str,
    pixels: tuple[int, int, int],
    color_by: str = "material",
    vtk_format: str = "vti",
) -> Path:
    ll, ur = _safe_bbox(model.geometry)

    plot = openmc.Plot()
    plot.type = "voxel"
    plot.filename = out_stem
    plot.origin = tuple((a + b) / 2.0 for a, b in zip(ll, ur))
    plot.width = tuple(b - a for a, b in zip(ll, ur))
    plot.pixels = pixels
    plot.color_by = color_by

    model.plots = openmc.Plots([plot])
    model.export_to_xml()
    openmc.plot_geometry(cwd=".")

    h5_path = Path(f"{out_stem}.h5")
    if not h5_path.exists():
        raise FileNotFoundError(f"Expected voxel file not produced: {h5_path}")

    if vtk_format not in {"vti", "vtk"}:
        raise ValueError("vtk_format must be 'vti' or 'vtk'")

    # OpenMC writes VTI from voxel H5. If the user requests legacy .vtk,
    # we convert the generated VTI to VTK using pyvista.
    vti_path = Path(f"{out_stem}.vti")
    openmc.voxel_to_vtk(str(h5_path), str(vti_path))

    if vtk_format == "vti":
        return vti_path

    vtk_path = Path(f"{out_stem}.vtk")
    ok = convert_vti_to_vtk(vti_path, vtk_path)
    if not ok:
        raise RuntimeError(
            "Requested --vtk-format vtk, but conversion from .vti to .vtk requires pyvista. "
            "Install pyvista or use --vtk-format vti."
        )
    return vtk_path


def convert_vti_to_vtk(vti_path: Path, vtk_path: Path) -> bool:
    try:
        import pyvista as pv  # type: ignore
    except Exception:
        return False

    grid = pv.read(str(vti_path))
    grid.save(str(vtk_path))
    return True


def maybe_convert_vti_to_stl(vtk_or_vti_path: Path, stl_path: Path) -> bool:
    try:
        import pyvista as pv  # type: ignore
    except Exception:
        return False

    grid = pv.read(str(vtk_or_vti_path))

    # Find first scalar array available
    arr_name = None
    for name in grid.array_names:
        arr_name = name
        break
    if arr_name is None:
        return False

    surf = grid.contour(isosurfaces=[0.5], scalars=arr_name)
    surf.save(str(stl_path))
    return True


def _canonical_volume_name(material_name: str) -> str:
    n = (material_name or "").lower()
    if "li4sio4" in n or "breeder" in n:
        return "breeder"
    if n == "be" or " beryll" in n or n.startswith("be"):
        return "be"
    if "eurofer" in n or "steel" in n:
        return "steel"
    if "tungsten" in n or "armour" in n or "armor" in n:
        return "armour"
    if "he" in n or "cool" in n:
        return "coolant"
    return material_name or "unknown"




def annotate_mesh_with_volume_tags(vtk_or_vti_path: Path, model: openmc.Model, out_stem: str) -> bool:
    """Embed numeric volume tags into the VTK/VTI cell data for MOOSE usage."""
    try:
        import numpy as np
        import pyvista as pv  # type: ignore
    except Exception:
        return False

    grid = pv.read(str(vtk_or_vti_path))
    if not grid.array_names:
        return False

    scalar_name = grid.array_names[0]
    arr = np.asarray(grid[scalar_name]).astype(int).ravel()

    id_to_name = {int(m.id): (m.name or f"material_{m.id}") for m in (model.materials or [])}

    canonical_ids = {"breeder": 1, "be": 2, "steel": 3, "armour": 4, "coolant": 5}
    volume_tag = np.zeros_like(arr)
    material_name = np.empty(arr.shape, dtype=object)

    for i, mid in enumerate(arr):
        mname = id_to_name.get(int(mid), "unknown")
        material_name[i] = mname
        cname = _canonical_volume_name(mname)
        volume_tag[i] = canonical_ids.get(cname, 0)

    grid.cell_data["openmc_material_id"] = arr
    grid.cell_data["volume_tag_id"] = volume_tag

    grid.save(str(vtk_or_vti_path))

    manifest = {
        "scalar_source": scalar_name,
        "volume_tag_id_map": canonical_ids,
        "material_id_to_name": {str(k): v for k, v in id_to_name.items()},
        "note": "Use volume_tag_id (or openmc_material_id) in downstream mesh-block generation.",
    }
    out_dir = Path(f"{out_stem}_named")
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "volume_tag_manifest.json").open("w", encoding="utf-8") as fp:
        json.dump(manifest, fp, indent=2)
    return True
def export_paraview_named_sets(vtk_or_vti_path: Path, model: openmc.Model, out_stem: str) -> bool:
    """Export named volume/boundary files for ParaView identification.

    Requires pyvista.
    """
    try:
        import numpy as np
        import pyvista as pv  # type: ignore
    except Exception:
        return False

    grid = pv.read(str(vtk_or_vti_path))
    if not grid.array_names:
        return False

    # Use first scalar array (material IDs for color_by=material).
    scalar_name = grid.array_names[0]
    arr = np.asarray(grid[scalar_name]).ravel()

    id_to_name = {int(m.id): (m.name or f"material_{m.id}") for m in (model.materials or [])}

    out_dir = Path(f"{out_stem}_named")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Volume blocks by canonical names.
    canonical_to_ids = {}
    for mid, mname in id_to_name.items():
        canonical = _canonical_volume_name(mname)
        canonical_to_ids.setdefault(canonical, []).append(mid)

    written = {"volumes": {}, "boundaries": {}}
    for canonical, ids in canonical_to_ids.items():
        mask = np.isin(arr.astype(int), np.array(ids, dtype=int))
        cid = np.where(mask)[0]
        if cid.size < 1:
            continue
        block = grid.extract_cells(cid)
        out = out_dir / f"volume_{canonical}.vtu"
        block.save(str(out))
        written["volumes"][canonical] = str(out)

        # cooling_surface from coolant volume surface if present
        if canonical == "coolant":
            surf = block.extract_surface()
            sfile = out_dir / "boundary_cooling_surface.vtp"
            surf.save(str(sfile))
            written["boundaries"]["cooling_surface"] = str(sfile)

    # vacuum_surface and symmetry surfaces from global outer surface.
    outer = grid.extract_surface()
    vac_file = out_dir / "boundary_vacuum_surface.vtp"
    outer.save(str(vac_file))
    written["boundaries"]["vacuum_surface"] = str(vac_file)

    b = outer.bounds
    x0, x1, y0, y1, z0, z1 = b
    tol = max(x1 - x0, y1 - y0, z1 - z0) * 1e-6
    pts = np.asarray(outer.points)

    def save_plane(mask, name):
        idx = np.where(mask)[0]
        if idx.size < 3:
            return
        sub = outer.extract_points(idx, adjacent_cells=True)
        f = out_dir / f"boundary_{name}.vtp"
        sub.save(str(f))
        written["boundaries"][name] = str(f)

    save_plane(np.abs(pts[:, 1] - y0) < tol, "symmetry_ymin")
    save_plane(np.abs(pts[:, 1] - y1) < tol, "symmetry_ymax")
    save_plane(np.abs(pts[:, 2] - z0) < tol, "symmetry_zmin")
    save_plane(np.abs(pts[:, 2] - z1) < tol, "symmetry_zmax")

    with (out_dir / "named_sets_manifest.json").open("w", encoding="utf-8") as fp:
        json.dump(written, fp, indent=2)

    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert OpenMC model XML to CAD-like mesh outputs.")
    parser.add_argument("--model-xml", default="model.xml", help="Path to OpenMC model.xml")
    parser.add_argument("--geometry-xml", default=None, help="Fallback geometry.xml path")
    parser.add_argument("--materials-xml", default=None, help="Fallback materials.xml path")
    parser.add_argument("--out-stem", default="openmc_geometry_mesh", help="Output stem for .h5/.vti/.stl")
    parser.add_argument("--nx", type=int, default=200, help="Voxel cells in x")
    parser.add_argument("--ny", type=int, default=120, help="Voxel cells in y")
    parser.add_argument("--nz", type=int, default=120, help="Voxel cells in z")
    parser.add_argument("--color-by", default="material", choices=["material", "cell"], help="Voxel coloring mode")
    parser.add_argument("--vtk-format", default="vti", choices=["vti", "vtk"], help="Output VTK file format")
    parser.add_argument("--stl", action="store_true", help="Also try to export STL (requires pyvista)")
    parser.add_argument(
        "--export-named-sets",
        action="store_true",
        help="Export named volume/boundary files for ParaView (requires pyvista)",
    )
    parser.add_argument(
        "--embed-volume-tags",
        action="store_true",
        help="Embed volume_tag_id/openmc_material_id arrays directly in .vtk/.vti (requires pyvista)",
    )
    args = parser.parse_args()

    model = _load_model(args.model_xml, args.geometry_xml, args.materials_xml)
    vtk_path = convert_model_to_vti(
        model,
        out_stem=args.out_stem,
        pixels=(args.nx, args.ny, args.nz),
        color_by=args.color_by,
        vtk_format=args.vtk_format,
    )
    print(f"Wrote {args.vtk_format.upper()}: {vtk_path}")

    if args.stl:
        stl_path = Path(f"{args.out_stem}.stl")
        ok = maybe_convert_vti_to_stl(vtk_path, stl_path)
        if ok:
            print(f"Wrote STL: {stl_path}")
        else:
            print("STL export skipped (pyvista not installed or no contourable scalar field).")

    if args.embed_volume_tags:
        ok = annotate_mesh_with_volume_tags(vtk_path, model, args.out_stem)
        if ok:
            print(f"Embedded volume_tag_id/openmc_material_id in: {vtk_path}")
        else:
            print("Volume-tag embedding skipped (pyvista unavailable or unreadable scalar field).")

    if args.export_named_sets:
        ok = export_paraview_named_sets(vtk_path, model, args.out_stem)
        if ok:
            print(f"Wrote named ParaView sets to: {args.out_stem}_named/")
        else:
            print("Named-set export skipped (pyvista unavailable or unreadable scalar field).")


if __name__ == "__main__":
    main()
