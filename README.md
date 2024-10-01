# VPN-DN-synapse-normalization

This repository accompanys "Morphology and synapse topography optimize linear encoding of synapse numbers in Drosophila looming responsive descending neurons." (doi: 10.1101/2024.04.24.591016.)

It contains all code necessary for the modeling pipeline, skeletonizing DNs, modeling their synaptic activations, mapping receptive fields of VPNs, and investigating retinotopy. 

Scripts are written in R, matlab and python, and require changing between languages for  processing of data, the outline and order of files is outlined below.
The packages and appropriate versions for each package is located in the packages file. 
General flow of modeling pipeline:
1. Generate skeleton model (skeletonization.py) and get synapse data (Obtaining_syns.py)
2. Manually proofread skeleton models using Neutube (https://neutracing.com/) / Meshlab (https://www.meshlab.net/) 
3. Map synapses to skeleton model for simulations (syn2morphMapper.py) / Map synapses for dendrograms (Mapping syn to nodes.R) ---------> (DNp0X_dendrogram.m)
4. Generate Anatomical Receptive fields of VPN populations (LC_visual_field_startup.R) --------> (LC_Receptive_field.R)
5. Visualization of synapses by receptive fields (Visualizing VPN by synapses by receptive fields.R)
6. Obtain ephys data  --------> Process in matlab (ephysReaderScript.m) --------> optimize passive properties (passivePropVisualizer.py)
7. Activate synapses (synActivation.py)
8. Visualize/Analyze/Plot synapse activation data (synDepolPlotter.py) 


Note: Direct simulation data has been ommitted, however available upon request.
Note: In the synDepolPlotter.py, functions will need to be adjusted to use the appropriate simulation data generated from the individual simulations (synActivation.py)

For any questions about the pipeline or questions regarding simulations please direct any correspondance to Anthony Moreno-Sanchez (Am4946@drexel.edu)
For mesh data reach out to Anthony.
