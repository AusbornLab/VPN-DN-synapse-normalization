import argparse
from pathlib import Path

import navis
from fafbseg import flywire


#in the terminal run the below command you can download the mesh for the specified root IDs and version.
#python Flywire_get_mesh.py --version 630 --ids 720575940622838154 

#List of DNs used from v630
#DNp01 720575940622838154
#DNp02: 720575940619654053
#DNp03: 720575940627645514
#DNp04: 720575940604954289
#DNp06: 720575940622673860

def main():
    parser = argparse.ArgumentParser(description="Download FlyWire FAFB meshes to OBJ")
    parser.add_argument("--version", type=int, choices=(630, 783), required=True)
    parser.add_argument(
        "--ids", type=int, nargs="+", default=[720575940622838154],
        help="One or more root IDs (default: 720575940622838154)",
    )
    parser.add_argument("--out", type=Path, default=Path("meshes"))
    parser.add_argument("--lod", type=int, choices=(0, 1, 2, 3), default=2)
    args = parser.parse_args()

    dataset = f"flat_{args.version}"
    args.out.mkdir(parents=True, exist_ok=True)

    for root_id in args.ids:
        path = args.out / f"{root_id}_v{args.version}.obj"
        print(f"Fetching {root_id} from {dataset} (LOD {args.lod})...")
        try:
            mesh = flywire.get_mesh_neuron(root_id, dataset=dataset, lod=args.lod)
            if not len(mesh.vertices) or not len(mesh.faces):
                raise ValueError("Mesh is empty; check the root ID and snapshot")
            navis.write_mesh(mesh, path)
        except Exception as exc:
            print(f"FAILED {root_id} ({dataset}): {exc}")
            continue
        print(f"Saved {path} ({len(mesh.vertices)} vertices, {len(mesh.faces)} faces)")


if __name__ == "__main__":
    main()