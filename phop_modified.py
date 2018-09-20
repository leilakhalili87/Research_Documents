import glob
import os
import re
#import matplotlib.pyplot as plt
import numpy
import numpy as np
import errno
import sys
from operator import add
import csv
numpy.set_printoptions(threshold=sys.maxsize)

#importing the path of all pos-displacement files
files_t = sorted(glob.glob('/home/leila/Leila_sndhard/codes/gb_mobility_test/gb_mob/mom_E0.1_T800/gb_dt .01/f0.01/msd/msd_data.*.txt'))
# sorting files based on the integer numbers
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
sort_nicely(files_t)
#######################################################################################################################
t_r=30  #t_r is t_r/2 in the paper
n_atoms=16315 # number of atoms
t_step=39999  # number of the time steps
t=0
j=0

# for i in range(1,(t_step-1)): # it should start from 1 (not zero) to n_step-1
for i in range(t_r,(t_step-t_r)): # it should start from 1 (not zero) to n_step-1
    print(i)
    ###########boundaris for first and last ten files
    if (i-t_r) < 0:
        t_r_min = 0
    else:
        t_r_min = int(i-t_r)
    if (i + t_r + 2)>t_step:
        t_r_max=t_step
    else:
        t_r_max=int((i + t_r + 2))
###########importing files inthe range of [t-t_r,t+t_r]
    file_1 = files_t[t_r_min:t_r_max]
    len_f = len(file_1)

    ff = np.zeros((len_f, n_atoms))
    ff_a = np.zeros((len_f, n_atoms))
    ff_b_1 = np.zeros((len_f, n_atoms))
    ff_a_1 = np.zeros((len_f, n_atoms))
    pos_x = np.zeros((len_f, n_atoms))
    pos_y = np.zeros((len_f, n_atoms))
    pos_z = np.zeros((len_f, n_atoms))
    dis_x = np.zeros((len_f, n_atoms))
    dis_y = np.zeros((len_f, n_atoms))
    dis_z = np.zeros((len_f, n_atoms))
    i_d = np.zeros((len_f, n_atoms))
##########reading x, y, z and displacements
    for idx, name in enumerate(file_1):
        try:
            with open(name) as f:
                result_x = []
                result_y = []
                result_z = []
                id=[]
                d_x = []
                d_y = []
                d_z = []
                lines = f.readlines()
                # l_line=len(lines)  #######333uncomment it
                l_line = n_atoms+9 #####################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Delete it
                line_1=lines[0]
                time_step=lines[1]
                line_2=lines[2]
                n_atomss=lines[3]
                line_3=lines[4]
                x_bound=lines[5]
                y_bound=lines[6]
                z_bound=lines[7]
                line_4=lines[8]

                for x in range(9, l_line):
                    y = lines[x]
                    id.append(int(y.split(' ')[0]))
                    result_x.append(y.split(' ')[1])
                    result_y.append(y.split(' ')[2])
                    result_z.append(y.split(' ')[3])
                    d_x.append(y.split(' ')[4])
                    d_y.append(y.split(' ')[5])
                    d_z.append(y.split(' ')[6])
        except IOError as exc:  # Not sure what error this is
            if exc.errno != errno.EISDIR:
                raise
        dis_x[idx, :] = d_x
        dis_y[idx, :] = d_y
        dis_z[idx, :] = d_z
        i_d[idx, :] = id
        if idx==0:
            pos_x[idx, :] = result_x
            pos_y[idx, :] = result_y
            pos_z[idx, :] = result_z

        else:
            pos_x[idx, :] = pos_x[(idx-1), :] + dis_x[idx, :]
            pos_y[idx, :] = pos_y[(idx-1), :] + dis_y[idx, :]
            pos_z[idx, :] = pos_z[(idx-1), :] + dis_z[idx, :]

############defining varialbles for A and B intervals

    if i<t_r:
        t_b=i+1   # t_b is maximum of A interval
        t_c=i    # t_c is the minimum of B intervval
        t_d=len_f-1  # is the max of B interval
    elif t_r <= i <= (t_step-t_r):
        t_b = t_r + 1
        t_c = t_r
        t_d = len_f - 1
    else:
        t_b=t_r+1
        t_c=t_r
        t_d=t_step
################calculating average in A and B interval
    a_x_i = (pos_x[0:t_b, :]).sum(axis=0) /t_b
    a_y_i = (pos_y[0:t_b, :]).sum(axis=0) / t_b
    a_z_i = (pos_z[0:t_b, :]).sum(axis=0) / t_b

    b_x_i = (pos_x[t_c:t_d, :]).sum(axis=0) /(t_d-t_c)
    b_y_i = (pos_y[t_c:t_d, :]).sum(axis=0) / (t_d-t_c)
    b_z_i = (pos_z[t_c:t_d, :]).sum(axis=0) /(t_d-t_c)
#########calculating average in r-r_a
    dif_a_x = (np.square(pos_x[0:t_b, :] - a_x_i)).sum(axis=0) / t_b
    dif_a_y = (np.square(pos_y[0:t_b, :] - a_y_i)).sum(axis=0) /t_b
    dif_a_z = (np.square(pos_z[0:t_b, :]- a_z_i)).sum(axis=0) / t_b
#calculating average in r-r_b
    dif_b_x = np.square(pos_x[t_c:t_d, :] - b_x_i).sum(axis=0) / (t_d-t_c)
    dif_b_y = np.square(pos_y[t_c:t_d, :] - b_y_i).sum(axis=0) /(t_d-t_c)
    dif_b_z = np.square(pos_z[t_c:t_d, :] - b_z_i).sum(axis=0) / (t_d-t_c)

    s_1 = dif_a_x + dif_a_y + dif_a_z
    s_2 = dif_b_x + dif_b_y + dif_b_z

    dd = np.multiply(s_1, s_2)
    qq= np.zeros(n_atoms)
    qq = np.sqrt(dd)

    ##########writing output file
    #text_file = open("/home/leila/Leila_sndhard/codes/gb_mobility_test/gb_mob/mom_E0.1_T800/gb_dt .01/f0.01/phop_30_modified/file%i.txt" % i, 'w')
    #text_file.write("%s\n" % qq   +'\n')
    text_file = open("/home/leila/Desktop/file%i.txt" % i,'w')
    text_file.write("ITEM: TIMESTEP\n")
    text_file.write("%s"  %time_step)
    text_file.write("ITEM: NUMBER OF ATOMS\n")
    text_file.write("%s" % n_atomss)
    text_file.write("ITEM: BOX BOUNDS ss pp pp\n")
    text_file.write("%s" % x_bound)
    text_file.write("%s" % y_bound)
    text_file.write("%s" % z_bound)
    text_file.write("ITEM: ATOMS  type x y z qq\n" )
    #text_file.write('time_step')
    np.savetxt(text_file, np.column_stack([i_d[i,:],pos_x[i,:], pos_y[i,:], pos_z[i,:],qq]))
    text_file.close()

    # with open('/home/leila/Leila_sndhard/codes/gb_mobility_test/gb_mob/mom_E0.1_T800/dE0.3/f0.01/phop_20/file.i%.csv'%i, 'w') as f:
    #     writer = csv.writer(f, delimiter='\t')
    #     writer.writerows(zip(qq))












