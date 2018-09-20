from ovito.io import import_file, export_file
from ovito.modifiers import PythonScriptModifier, CalculateDisplacementsModifier
import numpy
import os
import glob
import re

# Load input data and create an ObjectNode with a data pipeline.
filenames = sorted(glob.glob('/home/leila/Leila_sndhard/codes/gb_mobility_test/gb_mob/mom_E0.1_T800/gb_dt .01/f0.01/zdump/zdump.*.out'))

# make a function that sort the input file names
def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
sort_nicely(filenames)
n_f=len(filenames)  #number of dump files which is running step/10 (it is dumped every 10 time steps)

for i in range(n_f-1): ### 40000 is
    node = import_file (filenames[i+1])

# Calculate per-particle displacements with respect to initial simulation frame:
    dmod = CalculateDisplacementsModifier()
    dmod.reference.load(filenames[i])
    node.modifiers.append(dmod)
    node.compute()
    coordinates = node.output.particle_properties.position.array

# Export calculated MSD value to a text file and let OVITO's data pipeline do the rest:
    export_file(node, "/home/leila/Leila_sndhard/codes/gb_mobility_test/gb_mob/mom_E0.1_T800/gb_dt .01/f0.01/msd/msd_data.%i.txt"%i,format = "lammps_dump",
            columns=["Particle Identifier",  "Position.X","Position.Y", "Position.Z",
                     "Displacement.X", "Displacement.Y", "Displacement.Z"], frame=i)
