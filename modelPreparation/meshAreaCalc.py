import trimesh
from neuron import h

h.load_file("stdrun.hoc")

def Calculate_mesh_surface_area(mesh_file=None, dataset = None):
    
    # Load the mesh for the axons here, reach out to (am4946@drexel.edu for axonal mesh data, or access FANC data repository through CATMAID)
    mesh = trimesh.load(mesh_file)

    print("mesh loaded")
    print(type(mesh))
    surface_area = mesh.area

    print('Surface area:', surface_area, 'nm^2')
    if dataset == "FAFB":
        surface_area_um = surface_area / (4*4 * 40)/1000
    elif dataset == "Hemibrain":
        surface_area_um = surface_area / (8 *8 * 8)/1000
    elif dataset == "FANC":
        surface_area_um = surface_area / (4.3 *4.3 * 45)/1000
    print('Surface area:', surface_area_um, 'µm^2')


def skeleton_surface_area(swc_file):
    """
    Calculate the total surface area of a neuron from its SWC file.
    """
    h.load_file("import3d.hoc")

    # Import the morphology
    importer = h.Import3d_SWC_read()
    importer.input(swc_file)

    # Create a new cell from the morphology
    builder = h.Import3d_GUI(importer, 0)
    builder.instantiate(None)

    # Compute total surface area
    total_surface_area = 0.0

    for sec in h.allsec():
        for seg in sec:
            # NEURON automatically provides area in µm² for each segment
            total_surface_area += h.area(seg.x, sec=sec)

    print(f"Total skeleton neuron surface area: {total_surface_area:.2f} µm²")

swc_file = "datafiles/morphologyData/DNp03_morphData/Handtraced_DNp03_model.swc"  # Replace with your actual SWC file path
mesh_file = "meshes/DNp03_hemibrain_mesh.obj" # Replace with mesh path

print('Mesh data: DNp03')
Calculate_mesh_surface_area(mesh_file, dataset="Hemibrain") #this is an example for FANC dataset for the axon
skeleton_surface_area(swc_file)
