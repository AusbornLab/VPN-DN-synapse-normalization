from neuron import h, gui
from neuron.units import ms, mV
import time as clock
import numpy
import pandas as pd
from tkinter import Tk
import tkinter.filedialog as fd

#0 = unlabeled | 1 = soma | 2 = axon | 3 = dendrite
h.load_file("stdrun.hoc")

def loadSynData():
    #select synapse data file using tkinter
    Tk().withdraw()
    fd_title = "Select synapse data file to use for synapse mapping"
    syn_file = fd.askopenfilename(filetypes=[("xlsx file", "*.xlsx"), ("csv file","*.csv")], initialdir=r"datafiles/morphologyData", title = fd_title)
    synData_df = pd.read_excel(syn_file, dtype={'pre': str})
    #column 1 of the generated synapse data xlsx files have indexing that exists from when the data
    #was pulled from FAFB and before low-confidence synapses were filtered out (which is why it appears out of order)
    #this column could just be ignored, but it is effectively meaningless and therefore is dropped for clarity/simplicity
    synData_df = synData_df.drop(synData_df.columns[0], axis=1)
    return synData_df

def mapSynapses(outFilename=None, neuron_name = None):
    #function for mapping synapses from coordinates in space to sections of the morphology

    #creates a dataframe to parse through while mapping with one row for each coordinate point
    #of the morphology - each row has 4 cols: x coord, y coord, z coord, ID of section that the point belongs to
    i = 0
    coordData_df = pd.DataFrame(columns=['x', 'y', 'z', 'sec']) 
    for sec in h.allsec():
        for j in range(0, sec.n3d()):
            coordData_df.loc[i] = [sec.x3d(j), sec.y3d(j), sec.z3d(j), sec]
            i += 1
    
    #load synapse data
    synData_df = loadSynData()

    #start timer for progress readout during mapping
    time_start = clock.time()

    #create holding lists that will contain section ID and range var value for each mapped synapse
    #for more information on range variable check: https://neuron.yale.edu/ftp/ted/book/revisions/chap5indexedref.pdf
    synSecCol = []
    synSegRangeVarCol = []

    ### begin main loop of mapping process ###

    #this outer for loop iterates over each synapse and assigns to it two properties: the ID of a section of the morphology and a range variable called synSegRangeVar
    #that indicates how far along the section the synapse should be located

    #the inner for loop then iterates over every coordinate point saved in the coordData_df dataframe and calculates the distance
    #from that coordinate to the synapse, updating a placeholder variable (synSec) whenever appropriate (currDist < minDist)
    #TODO: add if nan, skip
    for synIndex, synRow in synData_df.iterrows():
        synCoords = []
        #since we are mapping the synapse onto the morphology, we use the postsynaptic coordinates and extract them here
        #since we will be doing math to find the nearest section, we store the coords in a numpy array
        synCoords.append(synRow.loc["post_x"])
        synCoords.append(synRow.loc["post_y"])
        synCoords.append(synRow.loc["post_z"])
        synCoords = numpy.array(synCoords)
        #initialize a dummy minimum distance value 
        minDist = numpy.NaN
        #inner for loop iterates over every coordinate point within the morphology and calculates the distance from that
        #point to the location of the current synapse
        #whenever a new minimum distance is found between a coordinate point and the current synapse, the synSec variable is updated
        #to reflect the current section that the synapse would be assigned to
        #once every coordinate point has been iterated over for a given synapse, the synapse is assigned to whatever section is 
        #stored in the synSec variable at the end of the for loop, and the code iterates over the outer for loop to 
        #proceed to the next synapse 
        for coordIndex, coordRow in coordData_df.iterrows():
            currPt = []
            #current coordinate points are extracted from the dataframe and stored in a numpy array
            currPt.append(coordRow.loc["x"])
            currPt.append(coordRow.loc["y"])
            currPt.append(coordRow.loc["z"])
            currPt = numpy.array(currPt)
            #calculate euclidean distance from the current coordinate point to the coordiante point of the synapse
            #current being evaluated in the outer for loop
            currDist = numpy.linalg.norm(currPt-synCoords)
            #if the current minDist value is NaN, this means that the current coordinate point being evaluated is the
            #first coordinate point being evaluated for the given synapse
            #therefore, the minDist and synSec values are updated accordingly
            if numpy.isnan(minDist):
                minDist = currDist
                synSec = coordRow.loc["sec"]
            else:
                #if the minDist value is not NaN, then we check if the distance from the current coordinate point to the 
                #synapse is less than the minDist value
                #if it is, we update minDist and synSec, otherwise we pass and move onto the next coordinate point
                if currDist < minDist:
                    minDist = currDist
                    synSec = coordRow.loc["sec"]
        #after every coordinate point has been iterated over, we exit the inner for loop, append the final value of synSec
        #to the holding list we created earlier, and move onto the next parameter for the current synapse: synSegRangeVar
        synSecCol.append(synSec)
        print("mapping to", synSec)
        #synSegRangeVar represents the distance along the section for the given synapse
        #we calculate this value for each synapse via the following process:
        
        #first we identify the start and end coordinate points that represent the section to which this synapse belongs
        synSecStart = numpy.array([synSec.x3d(0), synSec.y3d(0), synSec.z3d(0)])
        synSecEnd = numpy.array([synSec.x3d(synSec.n3d()-1), synSec.y3d(synSec.n3d()-1), synSec.z3d(synSec.n3d()-1)])
        
        #we then calculate the direction vector synSecDir by subtracting synSecStart from synSecEnd
        #synSecDir represents the direction in which the line segment (the section to which this synapse belongs to) extends
        synSecDir = synSecEnd-synSecStart

        #synToSecVec is computed by subtracting synSecStart from synCoords (the coordinate point to which the synapse belongs) 
        #synToSecVec represents a vector from the start of the synapse section to the synapse itself
        synToSecVec = synCoords - synSecStart

        #we then calculate the scalar projection of synToSecVec onto synSecDir
        #we divide the dot product of synToSecVec and synSecDir by the dot product of synSecDir with itself, the resulting projection value indicates how far along the synSecDir vector the point synCoords lies.
        synToSecVecDist = numpy.dot(synToSecVec, synSecDir) / numpy.dot(synSecDir, synSecDir)
        
        #using the projection value calculated in the previous step, we then determine the closest point on the section to the current synapse
        #we do so by multiplying the projection value synToSecVec by synSecDir and then adding it to synSecStart
        synClosestSecPt = synSecStart + synToSecVecDist * synSecDir

        #we then calculate synSegRangeVar which defines the relative position along the section of the point that is closest to the synapse (synClosestSecPt) 
        #we do this by dividing the Euclidean distance between synClosestSecPt and synSecStart by the length of the section, which is the Euclidean distance between synSecEnd and synSecStart.
        synSegRangeVar = numpy.linalg.norm(synClosestSecPt - synSecStart) / numpy.linalg.norm(synSecEnd - synSecStart)

        #if the synSegRangeVar variable falls above the maximum (1.0) or minimum (0.0) value, it is instead set to the maximum/minimum possible value
        if synSegRangeVar > 1:
            synSegRangeVar = 1
        if synSegRangeVar < 0:
            synSegRangeVar = 0
        
        #the synSegRangeVar value is then appended to the holding list created earlier for this paramter
        synSegRangeVarCol.append(synSegRangeVar)
        

        #progress print statement that reports how many synapses have been mapped and calculates what % of synapses still need to be mapped to the morphology 
        print("{}/{} synapses mapped | mapping {}% done".format(str(synIndex+1), str(len(synData_df)+1), str(round((synIndex+1)/(len(synData_df)+1)*100, 2))))
        time_current = clock.time()
        steps_elapsed = synIndex+1
        steps_remaining = (len(synData_df)+1) - steps_elapsed
        time_elapsed = time_current - time_start
        time_remaining = (time_elapsed/steps_elapsed) * steps_remaining
        #progress print out of rough time estimate remaining until mapping is complete
        print("estimated time remaining: {} minutes".format(str(round(time_remaining/60, 2))))

    #once all synapses have been assigned to a section and have been assigned a range varaible along that section, all the data is stored in two new columns that are added to the synData dataframe
    synData_df = synData_df.assign(mappedSection=synSecCol)
    synData_df = synData_df.assign(mappedSegRangeVar=synSegRangeVarCol)
    

    #the synData dataframe is then written to a csv where it is either saved with the datetime if no file name is given, or under the given filename if one is provided
    if outFilename is not None:
        if neuron_name is not None:
            outFilename = "datafiles/morphologyData/" + neuron_name + "_morphData/" + outFilename
        else:
            outFilename = "datafiles/morphologyData/" + outFilename
        synData_df.to_csv(outFilename, index=False)
    else:
        datetime = clock.strftime("%Y%m%d-%H%M%S")
        if neuron_name is not None:
            outFilename = "datafiles/morphologyData/synMap" + datetime
        else:
            outFilename = "datafiles/morphologyData/" + neuron_name + "_morphData/" + "synMap" + datetime
        synData_df.to_csv(outFilename, index=False)

def instantiate_swc(filename):
    ''' 
    Load swc file and instantiate it as cell
    Code source: https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3257
    '''

    # load helper library, included with Neuron
    h.load_file('import3d.hoc')

    # load data
    cell = h.Import3d_SWC_read()
    cell.input(filename)

    # # instantiate
    i3d = h.Import3d_GUI(cell,0)
    i3d.instantiate(None)
    return cell
    
def main():
    #select morphology data file using tkinter
    Tk().withdraw()
    fd_title = "Select morphology file to use for synapse mapping"
    morph_file = fd.askopenfilename(filetypes=[("swc file", "*.swc"), ("hoc file","*.hoc")], initialdir=r"datafiles/morphologyData", title=fd_title)

    #cell is instantiated so h.allsec() can be called in the mapSynapses function
    #this allows sections to be accessed for mapping purposes
    cell = instantiate_swc(morph_file)

    #if a custom filename is desired for the output synapse map .csv, the outFilename parameter below should be changed
    #to change the directory to which the synapse map .csv is written, make edits to the outFilename parameter at the end of the mapSynapses function
    mapSynapses(outFilename="synMap_DNp01_021324.csv", neuron_name="DNp01")

main()

