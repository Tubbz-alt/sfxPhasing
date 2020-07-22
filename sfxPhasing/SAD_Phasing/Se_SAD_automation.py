# This script is used to create the workflow for SHELXC/D and Crank2
# All paramters is parsed for SHELXC/D. Crank2 will take the result from SHELXD directly
# Partial result in the mid-steps is examined here but alse in SHELXC/D and Crank2 scripts.
# The result will be presented in the directory of the batch_submission.py
from __future__ import print_function
import sys
import subprocess
import os
import argparse
import re
import json
import ast
import numpy as np
import subprocess
import time
import shutil
import datetime


start_time = datetime.datetime.now()
current_path = os.getcwd()

if os.path.isfile("all_jobs.txt"):
    os.remove("all_jobs.txt")
###################################### Parse files and parameters ############################################    
parser= argparse.ArgumentParser()
parser.add_argument("-crank", "--crank-file-path", help="input the path for the crank file", type = str)
parser.add_argument("-shelx", "--shelx-file-path", help="input the path for the shelx file", type = str)
parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-hkl","--hkl-file", help="input shelx readeable file", type = str)
parser.add_argument("-seq","--sequence-file",help = "input the sequence file", type = str)
parser.add_argument("-FIND", "--number-of-atoms", default = 4, help='enter the number of atoms', type = str)
parser.add_argument("-resl","--resolution", default = '3.0', help = "input resolution", type = str)
parser.add_argument("-ESEL","--esel", default = '1.5', help = "input esel", type = str)
parser.add_argument("-thre","--threshold", default = '0.5', help = "input threshold of occupancy", type = str)
parser.add_argument("-DSUL","--disulfide", default = '4', help = "number of disulfide", type = str)
parser.add_argument("-SFAC", "--atom-type", help = "SAD atom type", type = str)
parser.add_argument("-MIND1","--mind-atom", default = '-3.5', help = "input minimum distance between atoms", type = str)
parser.add_argument("-MIND2","--mind-symm", default = '2.2',help = "input minimum distance between symmetry", type = str)
parser.add_argument("-lresl","--low-resolution-cut",default = '999', help = "input low resolution cutoff",type = str)
parser.add_argument("-P", "--path", help = "input the orginal path", type = str)
args = parser.parse_args()

if args.shelx_file_path:
    shelx = args.shelx_file_path
else:
    print("Please input the path for shelx python file")
    sys.exit()

if args.crank_file_path:
    crank = args.crank_file_path
else:
    print("Please input the path for the crank2 python file")
    sys.exit()

if args.reflection_mtz:
    reflectionFile = args.reflection_mtz
    print('Here is the reflection file')
    print(reflectionFile)
    job_file = reflectionFile.split('/')[-1]
    print('Here is the job file')
    print(job_file)
    #job_name = args.reflection_mtz.replace('.mtz','')
    job_name = job_file.replace('.mtz','')
else:
    print('Please input the mtz file for SHELXC/D')
    sys.exit()

if args.sequence_file:
    sequenceFile = args.sequence_file
    print('Here is the seq file')
    print(sequenceFile)
else:
    print('Please input sequence for crank2')
    sys.exit()

if args.resolution:
    resolution = args.resolution
    
if args.disulfide:
    dsul = args.disulfide
else:
    dsul = 0
    
if args.number_of_atoms:
    atom_find = args.number_of_atoms

if args.atom_type:
    atomType = args.atom_type
    
if args.esel:
    eSel = args.esel
print('Here is the reflection file again')    
print(reflectionFile)    
#shelx_CD = 'python SHELX_script.py  -rfl '+reflectionFile+' -MIND1 '+args.mind_atom+' -MIND2 '+args.mind_symm+' -resl '+str(resolution)+' -lresl '+args.low_resolution_cut+' -TEST 0 99 -ESEL '+eSel+' -SFAC '+atomType+' -NTRY 5000 -FIND '+str(atom_find)+' -DSUL '+str(dsul)+' -thre '+args.threshold
#shelx_CD = 'python '+shelx+ ' -rfl '+reflectionFile+' -MIND1 '+args.mind_atom+' -MIND2 '+args.mind_symm+' -resl '+str(resolution)+' -lresl '+args.low_resolution_cut+' -TEST 0 99 -ESEL '+eSel+' -SFAC '+atomType+' -NTRY 5000 -FIND '+str(atom_find)+' -DSUL '+str(dsul)+' -thre '+args.threshold
shelx_CD = 'python '+shelx+ ' -rfl '+reflectionFile+' -MIND1 '+args.mind_atom+' -MIND2 '+args.mind_symm+' -resl '+str(resolution)+' -lresl '+args.low_resolution_cut+' -TEST 0 99 -ESEL '+eSel+' -SFAC '+atomType+' -NTRY 5000 -FIND '+str(atom_find)+' -DSUL '+str(dsul)+' -thre '+args.threshold
print(shelx_CD)
f = open('shelxout.log','w')
os.system(shelx_CD)
#process = subprocess.Popen(shelx_CD, stdout = f, stderr=f, shell=True)

#f = open('shelxout.log','w')
CFOM_list = []
with open(job_name+'_fa.res','r') as f:
    for line in f:
        CFOM_list.append(line.rstrip('\n'))
for line in CFOM_list:
    if 'CFOM' in line:
        CFOM = line.split('CFOM')[-1].replace(' ','')

# To check whether SHELXC/D
if os.path.isfile(job_name+'_fa_cleaned.pdb'):
    new_pdb = job_name+'_fa_cleaned.pdb'
else:
    print('cleaned pdb has not been properly generated')
    os.chdir('../')
    sys.exit()

# run crank2
print('I m here')
#crank2_cl = 'python crank2_script.py -rfl '+reflectionFile+' -pdb '+new_pdb+' -seq '+sequenceFile+' -atype '+atomType+' -P '+args.path
#crank2_cl = 'python '+crank+' -rfl '+reflectionFile+' -pdb '+new_pdb+' -seq '+sequenceFile+' -atype '+atomType+' -P '+args.path
crank2_cl = 'python '+crank+' -rfl '+reflectionFile+' -pdb '+new_pdb+' -seq '+sequenceFile+' -atype '+atomType+' -P '+args.path
print(crank2_cl)
os.system(crank2_cl )
f = open('crank_output.log','w')
#process = subprocess.Popen(crank2_cl, stdout = f, stderr=f, shell=True)
#out,err = process.communicate()
#print(out)
#print(err)
# check whether crank2 works
end_time = datetime.datetime.now()
delta = end_time -  start_time
with open('times.log', 'w') as f1:
    f1.write('Total time for operation:' + str(float(delta.seconds) + float(delta.microseconds) / 1000000)+'\n')
os.chdir('crank2')
content = open('crank.loggraph','r').read()
for line  in content.split('\n'):
    if 'The final FOM is' in line:
        FOM = float(line.split('is ')[1])
        if FOM < 0.35:
            os.chdir(current_path)
            os.chdir('..')
            sys.exit()
        else:
            os.chdir(current_path)
            print(current_path+'/CFOM:'+str(CFOM)+'/FOM:'+str(FOM),file = open("FOM.txt","a"))
            if os.path.isfile('result.mtz'):
                new_mtz = 'result.mtz'
                print("Now it's time to start AutoBuild")
            else:
                print("Experimental reflection file has not been genereted")
                sys.exit()
print('I m exiting')

########################################## retrieve results ##############################################
# for i in os.listdir(os.getcwd()):
#     file_name = i.split('.')
#     if file_name[-1] == 'log':
#         if file_name[0] != 'shelxc' and file_name[0] != 'logfile':
#             log_file = i
            
# log_content = []
# with open(log_file) as f:
#     for line in f:
#         log_content.append(line.rstrip('\n'))
        
# R_list=[]
# R_free_list=[]
# if 'Majority of model was successfully built!' in log_content:
#     for i in log_content:
#         if 'R factor after refinement is ' in i:
#             R_list.append(i.split(' is ')[-1])
#         if 'R-free factor after refinement is ' in i:
#             R_free_list.append(i.split(' is ')[-1])


#     os.chdir(args.path)          
#     print(current_path.replace(args.path,'')+'/R:'+R_list[-1]+' ,R_free:'+R_free_list[-1],file=open('final_result.txt','a'))
    
            
############################################ To be modified !##############################################
# auto_build_cl = 'python autobuild.py -rfl '+new_mtz+' -orfl reflectionFile -seq '+sequenceFile+' -rfff 0.05 -nproc 12'

# os.system(auto_build_cl)

# os.chdir('AutoBuild_run_2_/')

# list=[]
# with open('AutoBuild_summary.dat','r') as f:
#     for line in f:
#         list.append(line.rstrip('\n'))

# os.chdir('..')
# for i in list:
#     if 'Best solution on cycle' in i:
#         print(os.getcwd()+'/'+i, file = open("best_solution_R_Rfree.txt",'a'))


