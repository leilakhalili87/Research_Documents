from ovito.io import import_file, export_file
from ovito.data import NearestNeighborFinder
from ovito.modifiers import PythonScriptModifier
import numpy as np
import numpy
import glob
import re


filenames = sorted(glob.glob('/home/leila/Leila_sndhard/codes/gb_mobility_test/gb_mob/mom_E0.1_T800/dE0.07_centerofmass_correct/f0.01/zdump_modified_com/*.txt'))
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
n_f=len(filenames)
print(n_f)
# The lattice constant of the FCC crystal:
lattice_parameter = 3.570856

a = np.asarray([
    (0.000000,	2.524977,	-0.000000),
    (0.000000,	-2.524977,	-0.000000),
    (0.728898,	1.262488,	2.061635),
    (-0.728898,	-1.262488,	-2.061635),
    (2.186694,	1.262488,	0.000000),
    (-2.186694, -1.262488, 0.000000),
    (1.457796,	0.000000,	-2.061635),
    (-1.457796, 0.000000, 2.061635),
    (2.186694,	-1.262488,	0.000000),
    (-2.186694, 1.262488, 0.000000),
    (-0.728898,	1.262488,	-2.061635),
    (0.728898, -1.262488, 2.061635)
])
aa=len(a)


for i in range(aa):
 y=np.linalg.norm(a[i, :])
 a[i,:]=a[i,:]/y

reference_vectors=a
# Rescale ideal lattice vectors with lattice constant.
reference_vectors *= lattice_parameter

# The number of neighbors to take into account per atom:
num_neighbors = len(reference_vectors)

#############################################################################################################
def modify1(frame, input, output):
    # Show a text in the status bar:
    yield "Calculating order parameters"

    # Create output property.
    order_param = output.create_user_particle_property(
        "Order Parameter", "float").marray

    # Prepare neighbor lists.
    neigh_finder = NearestNeighborFinder(num_neighbors, input)

    # Loop over all input particles
    nparticles = input.number_of_particles
    for i in range(nparticles):

        # Update progress percentage indicator
        yield (i / nparticles)

        oparam = 0.0  # The order parameter of the current atom

        # Loop over neighbors of current atom.
        for neigh in neigh_finder.find(i):
            # Compute squared deviation of neighbor vector from every
            # reference vector.
            squared_deviations = np.linalg.norm(
                reference_vectors - neigh.delta, axis=1) ** 2

            # Sum up the contribution from the best-matching vector.
            oparam += np.min(squared_deviations)

        # Store result in output particle property.
        order_param[i] = oparam / num_neighbors
#################################################################################
for j in range(n_f): ### 40000 is
    node = import_file (filenames[j])
    node.modifiers.append(PythonScriptModifier(function=modify1))
    node.compute()
    export_file(node, "/home/leila/Leila_sndhard/codes/gb_mobility_test/gb_mob/mom_E0.1_T800/dE0.07_centerofmass_correct/f0.01/order_p/%i.out"%j, format = "lammps_dump",
            columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z", "Order Parameter"], frame=j)
