#!/usr/bin/env python 
from ase import io, Atom, Atoms
import ase.db
from rdkit import Chem
import pybel
import os
from mordred import Calculator
import csv
import numpy as np
import pandas as pd
import progressbar
from time import sleep

from mordred import EState

#print out the process bar
bar = progressbar.ProgressBar(maxval= 12096, \
    widgets=[progressbar.Bar('=', 'Process bar[', ']'), ' ', progressbar.Percentage()])#12096

#Load in data                                                                                                                        
c = ase.db.connect('dssc.db')
data = pd.read_csv('data_visualization.csv').values
df = pd.DataFrame(data, columns = ['name', 'M', 'A', 'R1', 'R2', 'R3', 'E_gap'])

#define descriptor functions

Atom_types = ['sLi', 'ssBe', 'ssssBe', 'ssBH', 'sssB', 'ssssB', 'sCH3', 'dCH2', 'ssCH2', 'tCH', 'dsCH', 'aaCH', 'sssCH', 'ddC', 'tsC', 'dssC', 'aasC', 'aaaC', 'ssssC', 'sNH3', 'sNH2', 'ssNH2', 'dNH', 'ssNH', 'aaNH', 'tN', 'sssNH', 'dsN', 'aaN', 'sssN', 'ddsN', 'aasN', 'ssssN', 'sOH', 'dO', 'ssO', 'aaO', 'sF', 'sSiH3', 'ssSiH2', 'sssSiH', 'ssssSi', 'sPH2', 'ssPH', 'sssP', 'dsssP', 'sssssP', 'sSH', 'dS', 'ssS', 'aaS', 'dssS', 'ddssS', 'sCl', 'sGeH3', 'ssGeH2', 'sssGeH', 'ssssGe', 'sAsH2', 'ssAsH', 'sssAs', 'sssdAs', 'sssssAs', 'sSeH', 'dSe', 'ssSe', 'aaSe', 'dssSe', 'ddssSe', 'sBr', 'sSnH3', 'ssSnH2', 'sssSnH', 'ssssSn', 'sI', 'sPbH3', 'ssPbH2', 'sssPbH', 'ssssPb']

#Extract data 
M_list = set(('ZnP', 'TiOP', 'TiO2RP','FZnP', 'FH2P', 'FTiOP', 'FTiO2RP'))#, 'H2P'))
R_list = set(('Ph', 'DMP', 'TPA', 'MOTPA', 'TMP', 'DTA', 'DTBP', 'EthynPhM'))
A_list = set(('EthynPhA', '2CyanoPropenA', '2CarboxyPropenA', 'EthenThPCyanoAcryl', 'EthynBTDPhA', 'EthynDPhEPhA', 'EthynPhEPhA', 'EthynTPhEPhA', 'EthynPhDA', 'EthynFuA', 'EthynThPCyanoAcryl', 'EthynDThPCyanoAcryl', 'DThPCyanoAcryl', 'ThPCyanoAcryl', 'EthynThPA', 'EthynDThPA','rot-EthynPhA', 'rot-EthenThPCyanoAcryl', 'rot-EthynBTDPhA', 'rot-EthynDPhEPhA', 'rot-EthynPhEPhA', 'rot-EthynTPhEPhA', 'rot-EthynPhDA', 'rot-EthynFuA', 'rot-EthynThPCyanoAcryl', 'rot-EthynDThPCyanoAcryl', 'rot-EthynThPA', 'rot-EthynDThPA'))

Metal_list = ['ZnP', 'TiOP', 'TiO2RP','FZnP', 'FH2P', 'FTiOP', 'FTiO2RP','H2P']
Entry = []
bar.start()
n = 0
for Row in c.select(M= 'H2P'):
    Formula = Row.formula
    M= Row.M
    A = Row.A
    R1 = Row.R1
    R2 = Row.R2
    R3 = Row.R3
    E_gap = Row.E_gap
    
    #print 'Formula', Formula
    if A in A_list and R1 in R_list and R2 in R_list and R3 in R_list:        
        # Retrieve atom object form ase
        try:
            Atom =  c.get_atoms(Row.id)
            sample_name = str(Row.id+50161)
            
            io.write(sample_name +'.xyz', Atom, format='xyz')
            mol = pybel.readfile('xyz', sample_name+'.xyz').next()
            pybel.Outputfile('mol', sample_name +'.mol', overwrite = True).write(mol) # convert the input ase format to compatible 
            os.system('rm ' +  sample_name + '.xyz')       
        
            # Retrieve descriptors from mordred
            mols = Chem.MolFromMolFile(sample_name + '.mol')
            
            descriptor = []
            for k in range(len(Atom_types)):
                calc_ES = Calculator(EState.AtomTypeEState(type='sum', estate= Atom_types[k]))
                descriptor_atom = calc_ES(mols)[:1][0]
                descriptor.append(descriptor_atom)
                
            os.system('rm ' +  sample_name + '.mol')
        
            # Retrieve metal categorical index
            metals = np.zeros(len(Metal_list))
            index = Metal_list.index(M)
            metals[index] =1
            
            #rotation effect descriptor
            if A[:4]=='rot-':
                rot_A = [1]
            else:
                rot_A = [0]
        except:
            continue
        D = [Formula]+ [E_gap, M, A, R1, R2, R3] + descriptor + metals.tolist() + rot_A
        Entry.append(D)
        n+=1
        # Retrieve the E_gap for metal containning samples
        for i in M_list:
            Formula_ = df[(df['A']== A)&(df['R1']== R1)&(df['R2']== R2)&(df['R3']== R3)&(df['M']== i)]['name'].values
            E_gap_ = df[(df['A']== A)&(df['R1']== R1)&(df['R2']== R2)&(df['R3']== R3)&(df['M']== i)]['E_gap'].values.astype(float)
            if not len(Formula_)==0:
                #print 'Formula', Formula_[0]
                # Retrieve metal categorical index                                                                                            
                metals = np.zeros(len(Metal_list))
                index = Metal_list.index(i)
                metals[index] =1
        
                D_ = [Formula_[0]]+ [E_gap_[0], i, A, R1, R2, R3] + descriptor + metals.tolist() + rot_A
                Entry.append(D_)
                n+=1
                bar.update(n)
                sleep(0.1)
                #print 'n', n
bar.finish()
Entry_names = ['name', 'E_gap', 'M', 'A', 'R1', 'R2', 'R3'] + Atom_types + Metal_list +['Rot_A']
with open('E_State_Index.csv','w') as f:
    writer = csv.writer(f)
    writer.writerows([Entry_names])
    writer.writerows(Entry)

