import skeletor as sk # Version 1.2.1
import trimesh as tm # Version 3.13.4
import fafbseg # Version 3.0.10
import cloudvolume as cv # Version 10.4.0
from fafbseg import flywire

print("skeletor version:", sk.__version__)
print("trimesh version:", tm.__version__)
print("fafbseg version:", fafbseg.__version__)
print("cloudvolume version:", cv.__version__)

### Note for the code to work you will need access to the flywire dataset, which is accessible after training: 
# The below commmand only needs to be done once, where you set your token, this creates a seperate cloudvolume file storing your token for future use. 
#fafbseg.flywire.set_chunkedgraph_secret("INSERT YOUR TOKEN HERE")

#Before getting the mesh of your neurons assure that it is the most up to date version within the dataset
x = flywire.update_ids(720575940622838154)
print(x)

#Here insert the ID of the neuron of interest, I am using the neuron ID of DNp01 as an example here
neuron_mesh = flywire.get_mesh_neuron(720575940622838154, threads=1)


# NOTE: You can use the script from here below to skeletonize any mesh as well from any other dataset if you have the mesh as an .obj or an object file. 
# These commands are for importing your own mesh if you obtain it from anywhere else:

# mesh_path ="C:/Path to your mesh/Mesh.obj"
# neuron_mesh =tm.load_mesh(mesh_path)

# Plotting and view of the original mesh
# navis.plot3d(neuron_mesh)

# Fixing any small disconnections 
fixed_mesh = sk.pre.fix_mesh(neuron_mesh, remove_disconnected=25, inplace=False)

#Keeping the step_size low usually creates more nodes, but maintains a high level of detail for the dendrites of the neuron, 
# thicker processes will end up having issues and need to be edited by hand
skel = sk.skeletonize.by_wavefront(fixed_mesh, waves=1, step_size=1) 

#Cleans the skeleton from excess branching, or recenters nodes outside the mesh
skel_clean= sk.post.clean_up(skel, mesh = fixed_mesh)

#View of the skeleton
skel_clean.show()

## Saving the skeleton in a form to then edit.
skel_clean.save_swc('720575940622838154_skeletonized.swc')
