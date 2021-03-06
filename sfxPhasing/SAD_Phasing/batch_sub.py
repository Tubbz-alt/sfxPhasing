# This is the top file of the SAD pipeline. 
# It is used to carry out the batch submission
# This .py will do several automation and grid definition, and parse
# these information to SAD-automation file.
#
from __future__ import print_function
import sys
import subprocess
import os
import argparse
import re
import json
import ast
import numpy as np
import shutil
import random

############ This is the original directory which is the 'root' directory of your results ##########################
original_path = os.getcwd()
if os.path.isfile("final_result.txt"):
    os.remove("final_result.txt")
    
parser= argparse.ArgumentParser()

parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-seq","--sequence-file",help = "input the sequence file", type = str)
parser.add_argument("-SFAC","--atom-type", help = "input the name of atom of this SAD (case insensitive)" , type = str)
parser.add_argument("-q", "--queue", help = "input the computing queue you want to use", type = str)
parser.add_argument("-n", "--number-of-cores", help = "input the number of core you want to use", type = str)
parser.add_argument("-na", "--number-of-atoms", help = "input number of heavy anomalous scatterers", type = int)
parser.add_argument("-DSUL_R","--disulfide-range", nargs = '+', help='input the range of the disulfide range (optional)',type = str)
parser.add_argument("-RESOL_R","--resolution-range", nargs = '+', help='input the range of the resolution range (optional)',type = str)
parser.add_argument("-THRE_R","--threshold-range", nargs = '+', help='input the range of the threshold range (optional)',type = str)
parser.add_argument("-ATOM_R","--atom-range", nargs = '+', help='input the range of the atom number range (optional)',type = str)
parser.add_argument("-AutoBuild","--AutoBuild-polish", help = 'type N if you do not want it and type Y if you like it',type = str)
parser.add_argument("-Host","--host", help = 'must enter a host name, either type lcls or cori ',type = str)

args = parser.parse_args()

#parse reflection file
if args.reflection_mtz:
    reflectionFile = args.reflection_mtz
    job_name = args.reflection_mtz.replace('.mtz','')
else:
    print('Please input the mtz file')
    sys.exit()

if args.sequence_file:
    sequenceFile = os.path.abspath(args.sequence_file)
else:
    print('Please input sequence')
    sys.exit()
    
if args.atom_type:
    atomType = (args.atom_type).upper()
else:
    print('Please input the scattering atoms')
    sys.exit()

if args.queue:
    computeQueue = args.queue
else:
    print ('Please input the computing queue you want to use')
    sys.exit()
    
if args.number_of_cores:
    coreNumber = args.number_of_cores
else:
    print ('Please address the core number you want to use')
    sys.exit()
    
if args.host == 'lcls':
    job_submitter = 'bsub'
    python_run = 'python'
    core_selector = ' -n '
elif args.host == 'cori':
    job_submitter = 'srun -A m3506 -C haswell'
    python_run = 'python2.7'
    core_selector = ' --cores-per-socket '
else:
    print ('You have to enter the host: either lcls or cori')
    sys.exit()
###################### Setup the Grid Range #########################

print("This program utilize two softwares, CCP4i2 SHELXC/D and Crank2. Further refinement can be done by initiating Autobuild")

# Define the minimum distance between atoms and symmetry units. Users can change the definition in the parameter.json file in the 
# same directory.
with open('parameter.json') as json_file:
    data = json.load(json_file)
    for p in data['MIND1']:
        if p == atomType:
            mind_atom = str(data['MIND1'][atomType])
            
    for p in data['MIND2']:
        if p == atomType:
            mind_symm = str(data['MIND2'][atomType])
            
    for p in data["Low Resolution CutOff"]:
        if p == atomType:
            low_resolution_cut = str(data['Low Resolution CutOff'][atomType])
            
#Define the SE, S, and DSUL number from the sequence file
sequence=[]
with open(sequenceFile,'r') as f:
    for line in f:
        sequence.append(line.rstrip('\n'))
        
for line in sequence:
    if '>' in line:
        sequence.remove(line)

protein = ''.join(sequence)

single_S_or_SE_number=0
double_sulfur_number=0

for i in protein:
    if i=='M':
        single_S_or_SE_number+=1 #For S-SAD, M = S; For Se_SAD, M=SE
    elif i=='C':
        double_sulfur_number+=1
        
max_DSUL = 0
SE_num = 0
if single_S_or_SE_number == 0:
    if args.number_of_atoms:
        SE_num = args.number_of_atoms
        max_DSUL = double_sulfur_number/2
        max_S = single_S_or_SE_number+double_sulfur_number
    else:
        print("There is no methionine in this protein. Please enter the number of heavy atoms using -ATOM_R")
        #sys.exit()

else :
    max_DSUL = double_sulfur_number/2
    max_S = single_S_or_SE_number+double_sulfur_number
    max_SE = single_S_or_SE_number

#######################################################
process = subprocess.Popen('phenix.mtz.dump '+reflectionFile, 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()

split_out=out.splitlines()

for line in split_out:
    if 'Resolution range' in line:
        resolution = round(float(line.split(' ')[-1]),1)
        

DSUL_range = range(1,max_DSUL+1)
resolution_range = np.arange(resolution, resolution+1.0, 0.1) 

if atomType == 'S':
    atom_find = range(max_S/2, max_S+1)
else:
    if single_S_or_SE_number == 0:
        atom_find = range(SE_num-1, SE_num+2)
    else:
        atom_find = range(max_SE/2, max_SE+1)

    
    
thre_range = np.linspace(0.2,0.4,3) #(0.2,0.5,4)


################################# Update the grid range if asked by users ##############################

def get_range(x,y):
    if x == y:
        return np.array([x])
    else:
        return np.arange(x,y+0.1,0.1)

if args.disulfide_range:
    if atomType == 'S':
        DSUL_range = range(int(args.disulfide_range[0]),int(args.disulfide_range[1])+1)
    else:
        print('The atom type is not S. Input disulfide range will not be used')
else:
    print("Defualt DSUL search range will be used")
    
    
if args.resolution_range:
    resolution_range = get_range(round(float(args.resolution_range[0]),1),round(float(args.resolution_range[1]),1))
else:
    print("Defualt resolution search range will be used")
    
    
if args.atom_range:
    atom_find = range(int(args.atom_range[0]),int(args.atom_range[1])+1)
else:
    print("Defualt atom number search range will be used")
    
if args.threshold_range:
    thre_range = get_range(round(float(args.threshold_range[0]),1),round(float(args.threshold_range[1]),1))
else:
    print("Defualt threshold search range will be used")

################################# Randomize the job submission ########################################
print ("Disulfide search range: "+str(DSUL_range))
print ("Threshold cut-off range: "+str(thre_range))
print ("Resolution search range: "+str(resolution_range))
print ("Atom number search range:"+str(atom_find))

#sys.exit() #to be deleted



directory_list = []
command_list = []
#creat a list of directory and corresponding command
if max_DSUL > 0 and atomType == 'S':
    for dsul in DSUL_range:
        for thre in thre_range:
            for resolution in resolution_range:
                for number in atom_find:
                    directory = 'DSUL'+str(dsul)+'/threshold'+str(thre)+'/resolution'+str(resolution)+'/atom_number'+str(number)
                    directory_list.append(directory)

                    automation_cl = python_run+' Se_SAD_automation.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.3 -thre '+str(thre)+' -DSUL '+str(dsul)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm+' -lresl '+low_resolution_cut+' -P '+original_path+' -Host '+args.host
                    command_list.append(automation_cl)


#Do not consider DSUl parameter
elif max_DSUL == 0 or atomType != 'S':
    print ('0 disulfide bond found or the atom type you are looking for is not Sulfur. Do not consider the grid search for disulfied number')
    for thre in thre_range:
        for resolution in resolution_range:
            for number in atom_find:
                directory = 'threshold'+str(thre)+'/resolution'+str(resolution)+'/atom_number'+str(number)
                directory_list.append(directory)


                automation_cl = python_run+' Se_SAD_automation.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.5 -thre '+str(thre)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm+' -P '+original_path+' -Host '+args.host
                command_list.append(automation_cl)

# Shuffle the jobs
matching = list(zip(directory_list,command_list))
          
random.shuffle(matching) 

directory_list,command_list = zip(*matching)

# run jobs
for i in range(len(directory_list)):
    os.system('mkdir -p '+directory_list[i])
    os.system('cp Se_SAD_automation.py SHELX_script.py crank2_script.py autobuild.py '+sequenceFile+' '+reflectionFile+' '+directory_list[i])  
    os.chdir("./"+directory_list[i])
    os.system(job_submitter+' -q '+computeQueue+core_selector+coreNumber+' -o %J.log '+command_list[i])
    os.chdir(original_path)    
    
    
#################### Wait for to start Autobuild to polish the model ###########################
##helper method

##Jobs count method 
def job_count():
    finished_jobs = []
    try:
        with open('final_result.txt','r') as f:
            for line in f:
                finished_jobs.append(line.replace('\n',''))

        return len(finished_jobs)
    except:
        return 0

##Job selected to do autobuild
def case_select ():
    results = []
    with open('final_result.txt','r') as f:
        for line in f:
            results.append(line.replace('\n',''))

    R_score = []
    for result in results:
        R_score.append(round(float(result.split('R_free:')[-1].split('/')[0]),3))

    R_score = np.asarray(R_score)
    polish_cases_Rfree = np.where(R_score == R_score.min())[0]

    polish_cases_0 = []
    for i in polish_cases_Rfree:
        polish_cases_0.append(results[i])

    Residue_score = []
    for i in polish_cases_0:
        Residue_score.append(int(i.split('Residue:')[-1]))

    Residue_score = np.asarray(Residue_score)
    polish_cases_residue = np.where(Residue_score == Residue_score.max())[0]

    polish_cases_1 = []
    for i in polish_cases_residue:
        polish_cases_1.append(polish_cases_0[i])

    #return polish_cases_1
    resolution_filter = []
    if len(polish_cases_1) > 1:
        for i in polish_cases_1:
            resolution_filter.append(float(i.split('resolution')[-1].split('/atom')[0]))
        resolution_filter = np.asarray(resolution_filter)
        polish_cases_resolution = np.where(resolution_filter == resolution_filter.min())[0]
    
    
        polish_cases_2 = []
        for i in polish_cases_resolution:
            polish_cases_2.append(polish_cases_1[i])
            #print("have to pick randomly")
        if len(polish_cases_2) > 1:
            print("have to pick randomly")
            select_case = polish_cases_2[random.randint(0,len(polish_cases_2)-1)]

        else: select_case = polish_cases_2[0]
            
    else:
        select_case = polish_cases_1[0]

    print ("The best grid right now is "+select_case.split('R:')[0])
    print ("It has the score R:"+select_case.split('R:')[-1])
    return select_case.split('R:')[0]


Total_jobs = len(directory_list)
Half_total_jobs = Total_jobs // 2
percent97_total_jobs = Total_jobs * 97 //100
half_finished = False
percent97_finished = False

while half_finished == False:
    
    if job_count() < max(1,Half_total_jobs):
        pass
    else:

        selected_job_directory1 = case_select()
        if args.AutoBuild_polish != 'N':
            os.system("mkdir Autobuild1")
            os.system("cp autobuild.py "+sequenceFile+" "+reflectionFile+" "+selected_job_directory1+"result.pdb Autobuild1")
            os.chdir("Autobuild1")
            autobuild_cl = python_run+' autobuild.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -rfff 0.05 -nproc '+str(coreNumber)+' -pdb result.pdb'
            os.system(job_submitter+' -q '+computeQueue+core_selector+coreNumber+' -o %J.log '+autobuild_cl)
            os.chdir(original_path)
            half_finished = True


while percent97_finished == False:
    
    if job_count() < max(1,percent97_total_jobs):
        pass
    else:
        percent97_total_jobs = True
        selected_job_directory2 = case_select()
        if selected_job_directory2 == selected_job_directory1:
            break
        else:
            if args.AutoBuild_polish != 'N':
                os.system("mkdir Autobuild2")
                os.system("cp autobuild.py "+sequenceFile+" "+reflectionFile+" "+selected_job_directory2+"result.pdb Autobuild2")
                os.chdir("Autobuild2")
                autobuild_cl = python_run+' autobuild.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -rfff 0.05 -nproc '+coreNumber+' -pdb result.pdb'
                os.system(job_submitter+' -q '+computeQueue+core_selector+coreNumber+' -o %J.log '+autobuild_cl)
                os.chdir(original_path)
                percent97_finished = True








   
        
