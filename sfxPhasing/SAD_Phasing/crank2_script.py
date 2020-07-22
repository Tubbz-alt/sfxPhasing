from __future__ import print_function
import sys
import subprocess
import os
import argparse
import re
import json
import ast

current_path = os.getcwd()
#prevent repeatedly writing
if os.path.isfile("crank2.inp"):
    os.remove("crank2.inp")
#################################### Parse files and atom types ########################################
parser= argparse.ArgumentParser()


parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-pdb","--pdb-file",help="input the pdb file",type = str)
parser.add_argument("-seq","--sequence-file",help="input the sequence file",type = str)
parser.add_argument("-atype", "--atom-type", help="input the atom type", type = str)
parser.add_argument("-P", "--path", help = "input the orginal path", type = str)
args = parser.parse_args()

if args.reflection_mtz:
    reflection_file = args.reflection_mtz
    print('Here is the reflection file for crank2')
    print(reflection_file)
else:
    print('The reflection file is missing')
    sys.exit()

if args.pdb_file:
     pdb= args.pdb_file
else:
    print('The pdb file is missing')
    sys.exit
    
if args.sequence_file:
    seq_file = args.sequence_file
else:
    print('The sequence file is missing')
    sys.exit()

if args.atom_type:
    atomType = args.atom_type
################################################################################################################


################################### Modified the label of the input mtz file ####################################
os.system("cmtzsplit -mtzin "+reflection_file+" -mtzout reflection_out.mtz -colin 'F(+),SIGF(+),F(-),SIGF(-)' -colout Fplus,SIGFplus,Fminus,SIGFminus")
#################################################################################################################

################################### Extract protein chain information ###########################################
# By using xtriage, we retrieve
# monomer number per asymmetric unit, solvent content, and number of residues
process = subprocess.Popen('phenix.xtriage reflection_out.mtz '+seq_file, 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()

split_out=out.splitlines()

mylist = []
for i in range(len(split_out)):
    item = split_out[i].decode("utf-8")
    if item != '':
        mylist.append(split_out[i].decode("utf-8"))

for line in mylist:
    if 'Best guess' in line:
        best_guess = re.findall(r'\d+',line)[0]
        monomer_asu = best_guess
        
    if 'Crystallized molecule(s) defined as' and 'protein residues' in line:
        protein_residue_number = re.findall(r'\d+', line)[0]
        
    if 'Copies' and 'Solvent content' and 'Matthews coeff.' in line:
        statistical_title_index = mylist.index(line)
        Copies_slot = line.split('|')[1].replace(' ','')
        Solvent_content_slot = line.split('|')[2].replace(' ','')
        Matthews_coeff_slot = line.split('|')[3].replace(' ','')
        P_solvent_content_slot = line.split('|')[4].replace(' ','')


still_in_table = True
i = 2
while still_in_table:
    copy_check = float(mylist[statistical_title_index + i].split('|')[1].replace(' ',''))
    if copy_check == float(best_guess):
        solvent_content = mylist[statistical_title_index + i].split('|')[2].replace(' ','')
        still_in_table = False
    else:
        i+=1
    

num_atoms = str(4)
with open('Guessed_atom_number.txt','r') as f:
    for line in f:
        if 'exp_num_atoms' in line:
            num_atoms = line.split(" = ")[1].replace('\n','')
            print(type(num_atoms))
###########################################################################################################

###################################Draft Crank2.inp########################################################
#density modification run number is defined to be 50, dmcyc:50

line1 = 'fsigf plus dname=peak f=Fplus sigf=SIGFplus "file=reflection_out.mtz"'
line2 = 'fsigf minus dname=peak f=Fminus sigf=SIGFminus'
line3 = 'model substr atomtype='+atomType+' "file='+pdb+'" d_name=peak exp_num_atoms='+num_atoms+' monomers_asym='+monomer_asu+' residues_mon='+protein_residue_number+' solvent_content='+solvent_content
line4 = 'sequence "file='+seq_file+'"'
line5 = 'createfree no_output_to_next_step::True fraction::0.05'
line6 = 'handdet dmfull dm solomon phcomb multicomb'
line7 = 'dmfull dmcyc::50 threshold_stop::0.58 dm parrot phcomb refmac'
line8 = 'comb_phdmmb target::SAD maxbigcyc::100 exclude obj_from=0,typ=freeR dmfull dm parrot'

crank_inp_file = [line1, line2, line3, line4, line5, line6,line7,line8]

for i in crank_inp_file:
    print(i,file=open("crank2.inp", "a"))
#############################################################################################################

#################################### Run Crank2 ###############################################################
#command = 'python /reg/common/package/ccp4/ccp4-7.0/share/ccp4i2/pipelines/crank2/crank2/crank2.py --keyin crank2.inp --hklout result.mtz --xyzout result.pdb'
#command = 'python /project/projectdirs/m3506/ccp4-7.0/share/ccp4i2/pipelines/crank2/crank2/crank2.py --keyin crank2.inp --hklout result.mtz --xyzout result.pdb'
command = 'python /project/projectdirs/m3506/ccp4-7.0/share/ccp4i2/pipelines/crank2/crank2/crank2.py --keyin crank2.inp --hklout result.mtz --xyzout result.pdb'
f = open('crank2.log','w')
process = subprocess.Popen(command, 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()
f.write(out)
f.write(err)
#print(out)
#print(err)
split_out=out.splitlines()
#os.system(command)
# for line in split_out:
#     if 'Resolution range' in line:
#         resolution = round(float(line.split(' ')[-1]),1)


R_list=[]
R_free_list=[]
residue_report = []
Succeed = False
for line in split_out:
    if 'Majority of model was successfully built!' in line:
        Succeed = True
        for i in split_out:
            if 'R factor after refinement is ' in i:
                R_list.append(i.split(' is ')[-1])
            elif 'R-free factor after refinement is ' in i:
                R_free_list.append(i.split(' is ')[-1])
            elif ' fragments built.' in i:
                residue_report.append(i)
        break

#         os.chdir(args.path)          
#         print(current_path.replace(args.path,'')+'/R:'+R_list[-1]+' ,R_free:'+R_free_list[-1],file=open('final_result.txt','a'))
    else:
        Succeed = False
    
if Succeed == True:
    os.chdir(args.path)          
    print(current_path.replace(args.path+'/','')+'/R:'+R_list[-1]+'/R_free:'+R_free_list[-1]+'/Residue:'+residue_report[-1].split(' ')[0],file=open('final_result.txt','a'))
    

else:
  
    os.chdir(args.path)          
    print(current_path.replace(args.path+'/','')+'/R:0.55/R_free:0.55',file=open('final_result.txt','a'))





