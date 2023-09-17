#Initializes file paths and generates required directories
import os
import sys
os.environ['MKL_THREADING_LAYER'] = 'GNU'
GPU_NUMBER = [0]
os.environ["CUDA_VISIBLE_DEVICES"] = ",".join([str(s) for s in GPU_NUMBER])
os.environ["NCCL_DEBUG"] = "INFO"

filepath = os.getcwd()
homepath = os.path.expanduser('~')

if homepath not in filepath:
    
    if 'batchRun' not in filepath:
        homepath = filepath[:-15]
    else:
        homepath = filepath[:-25]
        filepath = filepath[:-10]
        
#Generates file for local job operation
try:
    os.mkdir(homepath + 'UnifiedPipeline/Pairs')
except:
    pass

#Package imports
from proteinTools import PT as p
import re
import csv
import ast
import random
from vina import Vina
import math
import traceback
import pandas as pd
import numpy as np
import pubchempy as pcp
import urllib
import urllib.request
import mygene 
from proteinTools import PT as p
import statistics
import argparse

os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'

#Regular expression for stripping numeric components from string
regex = re.compile('[0-9]')

#Regular expression for stripping non-numeric components from string
numconv = re.compile('^[0-9]')

# Obtains arguments
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', dest = 'file_name', type = str, help = 'Name of informational file that provides proteins/ligands for docking', default = '')
    parser.add_argument('-t', '--total_jobs', dest = "total_jobs", type = int, help = "Total number of jobs across sample space", default = 1)
    parser.add_argument('-i', '--instance', dest = 'instance', type = int, help = 'Current job instance of run', default = 1)
    parser.add_argument('-p', '--protein', dest = 'proteinID', type = str, help = 'Protein ID Type (PDB, Uniprot Accession, HGNC)', default = 'PDB')
    parser.add_argument('-l', '--ligand', dest = 'ligandID', type = str, help  = 'Ligand ID Type (SMILE, InChiKey', default = 'SMILE')
    args = parser.parse_args()
    
    return args

#Primary class for running pipeline
class run_pipeline:
    def __init__(self):
        
        #Generates output folder if one does not exist
        try:
            os.mkdir(homepath + 'UnifiedPipeline/Results')
        except:
            pass
            
        #Generates output folder for poses if one does not exist
        try:
            os.mkdir(homepath + 'UnifiedPipeline/OutPoses')
        except:
            pass
            
        # Reads user arguments
        args = get_arguments()
        
        #Determines job number
        self.totaljobs = args.total_jobs
        self.instance = args.instance
        
        #Reads configuration options
        testfile = args.file_name
        self.protID = args.proteinID
        self.ligID = args.ligandID

        #Generates list of acceptable trigrams and chains for future use
        res_dict = {'ALA':'A', 'ALX':'B', 'CYS': 'C', 'ASP' : 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS':'K', 'LEU':'L', 'MET': 'M', 'ASN':'N', 'PRO':'P', 'GLN': 'Q', 'ARG':'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR':'Y', 'GLX':'Z'}
        
        #Dictionary for atom masses
        self.atom_dict = {'H':1.01, 'C':12.01, 'O':16.00, 'N':14.01, 'A':12.01, 'P':30.97, 'F':18.998, 'FE':55.845, 'S':32.06, 'B':10.81, 'K':39.1, 'I':126.904, 'BR':79.904, 'CL':35.453, 'CA':40.08, 'NA':22.99, 'MG':24.305, 'AL':26.98, 'CR':51.996, 'NE':20.179, 'BE':9.01, 'FE':55.847, 'CO':58.933,'AG':107.868, 'CD':112.41, 'NI':58.693, 'ZN':65.39, 'BE':9.0122, 'IN':114.818, 'SI':28.085, 'SC':44.956, 'TI':47.867, 'V':50.941, 'MN':54.938, 'CU':63.546, 'GA':59.723, 'GE':72.64, 'SE':78.96, 'KR':83.8, 'ZR':91.224, 'NB':92.906, 'PD':106.42, 'SN':118.71, 'SB':121.76, 'XE':131.293, 'BA':137.327, 'LA':138.91, 'LI':6.941, 'HG':200.59, 'PB':207.2, 'BI':208.98, 'PO':209, 'TI':204.3833, 'AU':196.9665, 'IR':192.217, 'PT':195.078, 'RE':186.207, 'W':183.84, 'TA':180.948, 'YB':173.04, 'EU':151.964, 'ND':144.25, 'CE':140.116, 'TH':232.04, 'U':238.029, 'PU':244, 'FR':223, 'PA':231.04, 'HO':164.93, 'SM':150.36, 'PR':140.908, 'TE':127.6, 'TC':98, 'Y':88.906}
        
        #Creates protein and ligand lists and divides work up based on job number
        os.chdir(homepath + 'UnifiedPipeline')
        data = pd.read_csv(testfile)
        try:
            self.proteins = data['Protein'].tolist()
        except:
            try:
                self.proteins = data['Protein ID'].tolist()
            except:
                try:
                    self.proteins = data['Proteins'].tolist()
                except:
                    print('Invalid protein ID identifier column name! Must be Protein, Protein ID, or Proteins. Aborting program.')
                    sys.exit()
    
        if 'SMILE' in self.ligID.upper() or 'SMI' in self.ligID.upper():
            try:
                self.ligands = data['SMILE'].tolist()
            except:
                try:
                    self.ligands = data['smi'].tolist()
                except:
                    try:
                        self.ligands = data['Ligand SMILE'].tolist()
                    except:
                        print(data.columns)
                        print('Invalid ligand ID identifiers column name! Must be SMILE, smi, or Ligand SMILE. Aborting program.')
                        sys.exit()
        elif 'INCHIKEY' in self.ligID.upper() or 'ICK' in self.ligID.upper():
            try:
                self.ligands = data['InChiKey'].tolist()
            except:
                try:
                    self.ligands = data['ick'].tolist()
                except:
                    try:
                        self.ligands = data['Ligands'].tolist()
                    except:
                        print('Invalid ligand ID identifiers! Must be InChiKey, ick, or Ligands. Aborting program.')
                        sys.exit()
        else:
            print('Invalid ligand ID Type! Must be an InChiKey or a SMILE sequence')
            sys.exit()
            
        #Creates job queue based on job type parameters 
        datalen = len(self.ligands)        
        chunklen = datalen // self.totaljobs
        if self.instance == self.totaljobs:
            chunklen += datalen % self.totaljobs
            self.ligands = self.ligands[datalen - chunklen:]
            self.proteins = self.proteins[datalen - chunklen:]
        else:  
            self.ligands = self.ligands[(self.instance - 1) * chunklen:self.instance * chunklen]
            self.proteins = self.proteins[(self.instance - 1) * chunklen:self.instance * chunklen]
        
        #Rewrites PDB IDs incorrectly converted to exponential numbers from .csv format
        temp_proteins = []    
        for protein in self.proteins:
            try:
                if protein[4] == 'E' and len(protein) > 4:
                    try:
                        v = int(protein[0] + protein[6:8])
                        temp_proteins.append(protein[0] + protein[4] + protein[6:8])
                    except:
                        temp_proteins.append(protein)
                else:
                    temp_proteins.append(protein)
                    
            except:
                temp_proteins.append(protein)
        self.proteins = temp_proteins    
        self.joblen = len(self.proteins)
        print(str(self.joblen) + ' job(s) queued up!')
        
      #Downloads and prepares file suite for protein and ligand
    def generate_files(self, strip_ligands = False):     
        #Iterates through unique proteins in job list and generates .pdb and .pdbqt files
        self.protein_classes = []
        
        for protnum, protein in enumerate(set(self.proteins)):
            try:
                os.chdir(filepath + '/Pairs/' + protein)    
            except: 
                os.mkdir(filepath + '/Pairs/' + protein)
                os.chdir(filepath + '/Pairs/' + protein)
            
            print('Fetching data for ' + protein)
            if 'UNIPROT' in self.protID.upper():
                prot = p.Protein(protein, 'Uniprot')
            elif 'GENE' in self.protID.upper():
                prot = p.Protein(protein, 'Gene')
            else:
                prot = p.Protein(protein, 'PDB')
                
            prot.download(os.getcwd())
            self.protein_classes.append(prot)
        
            if strip_ligands == True:
                # Strips all ligands
                prot.strip_ligands()
            else:
                # Only strips non-cofactor ligands
                prot.strip_primary()
                
           #Generates partially charged 3D pdbqt structure
            try:
                # Clears existing file if needed
                try:
                    os.remove(protein + '.pdbqt')
                except:
                    pass
                
                python_executable = sys.executable

                # Check if the Python executable is within a Conda environment
                if "conda" in python_executable:
                    # Extract the Conda environment path from the executable path
                    conda_env_path = os.path.dirname(os.path.dirname(python_executable))
                
                    # Construct the path to the 'bin' directory within the Conda environment
                    bin_dir = os.path.join(conda_env_path, 'bin')
                    
                pdbqtconv = 'python2.7 ' + bin_dir + '/prepare_receptor4.py -r ' + protein + '.pdb -o ' + protein + '.pdbqt'  
              
                os.system(pdbqtconv)
                with open(protein + '.pdbqt', 'r') as f:
                    pass
            except:
               
                with open(protein + '.pdb') as f:
                    data = f.readlines()
                data = [i for i in data if i[0:4] == 'ATOM']
                with open(protein + '.pdb', 'w') as f:
                    for line in data:
                        f.write(line)
                pdbqtconv = 'obabel -ipdb ' + protein + '.pdb -opdbqt -O ' + protein + '1.pdbqt'
                os.system(pdbqtconv)

                try:
                    os.rename(protein + '1.pdbqt', protein + '.pdbqt')
                except:
                    pass
                      
        #Iterates through ligands in list and generates relevant directories in the associated protein directory, as well as generating .pdb and .pdbqt files
        self.ligand_smiles, self.ligand_names = [], []
        for ligandnum, ligand in enumerate(self.ligands):
            formatted_dir = re.sub('[()]', '', ligand)
            os.chdir(filepath + '/Pairs/' + self.proteins[ligandnum])
            try:
                os.mkdir(formatted_dir)
                os.chdir(formatted_dir)
            except:
                os.chdir(formatted_dir)
                
            if 'SMILE' in self.ligID.upper() or 'SMI' in self.ligID.upper():
                SMILE = ligand
                try:
                    with open(ligand + '.pdb', 'r') as f:
                        pass
                except:
                    with open('smile' + str(self.instance) + '.txt', 'w') as f:
                        f.write(SMILE)
                    os.system('obabel -ismi smile' + str(self.instance) + '.txt -O ' + 'lig' + str(self.instance) + '.pdb --gen3d')
                    os.remove('smile' + str(self.instance) + '.txt')
            else:
                pcp.download('SDF', 'lig' + str(self.instance) + '.sdf', ligand, 'inchikey', overwrite = True)
                pdbconv = 'obabel -isdf lig' + str(self.instance) + '.sdf -O lig' + str(self.instance) + '.pdb --gen3d'
                os.system(pdbconv)
                os.remove('lig' + str(self.instance) + '.sdf')
                SMILE = pcp.get_compounds(ligand, 'inchikey', as_dataframe = True)["isomeric_smiles"].to_list()[0]
            self.ligand_names.append(ligand)
            self.ligand_smiles.append(SMILE)
            self.ligands[ligandnum] = formatted_dir
            ligand = formatted_dir
        
            
            pdbqtconv = 'obabel -ipdb lig' + str(self.instance) + '.pdb -O ' + ligand + '.pdbqt --gen3d'
            os.system(pdbqtconv)
                
     #Calculates residue centers for each residue in every protein, as well as radius of gyration for each ligand
    def atomize_files(self):
        self.FASTA = []
        for protnum, protein in enumerate(set(self.proteins)):
    
            #Obtains FASTA sequence from protein
            self.FASTA.append(self.protein_classes[protnum].FASTA)

        #Obtains ligand radius of gyration for each ligand, then maps it to a dictionary
        self.ligand_radii = []
        for lignum, ligand in enumerate(self.ligand_names):
            ligand_mol = p.ligand(ligand)
            os.chdir(filepath + '/Pairs/' + self.proteins[lignum] + '/' + self.ligands[lignum])
            ligand_mol.download(directory = os.getcwd(), convert_pdb = True)
            self.ligand_radii.append(ligand_mol.radius)
            
    #Performs docking
    # Metric can be RMSD, Overlap (Overlap between docked positions), or Distance (deviation from true position) or None (no comparison performed)
    def run_docking(self, metric = None):
        self.metric = metric
        
        self.blind_outposes = []
        self.blind_distance = []
        self.blind_affinities = []
        
        #Instantiates outpose directory if not created
        try:
            os.mkdir(homepath + 'UnifiedPipeline/OutPoses')
            os.chdir(homepath + 'UnifiedPipeline/OutPoses')
        except:
            os.chdir(homepath + 'UnifiedPipeline/OutPoses')
        
        self.binding_affinities = []

        for jobnum in range(0, self.joblen):    
            v = Vina(sf_name = 'vina')
            protein = self.proteins[jobnum]
            ligand = self.ligands[jobnum]
           
            #Attempts to create folder for ligand poses
            try:
                os.mkdir(f'{homepath}UnifiedPipeline/OutPoses/{protein}{ligand}')
            except:
                pass
            v.set_receptor(f'{homepath}UnifiedPipeline/Pairs/{protein}/{protein}.pdbqt')
            v.set_ligand_from_file(f'{homepath}UnifiedPipeline/Pairs/{protein}/{ligand}/{ligand}.pdbqt')    

            # Performs Blind Docking
            protein_class = self.protein_classes[jobnum]
            center = protein_class.center
            x, y, z = [center[0], center[0]], [center[1], center[1]], [center[2], center[2]]
            for res in protein_class.residue_list:
                for atom in res.atoms:
                    if atom.x > x[1]:
                        x[1] = atom.x
                    elif atom.x < x[0]:
                        x[0] = atom.x
                    if atom.y > y[1]:
                        y[1] = atom.y
                    elif atom.y < y[0]:
                        y[0] = atom.y
                    if atom.z > z[1]:
                        z[1] = atom.z
                    elif atom.z < z[0]:
                        z[0] = atom.z
            size = [x[1] - x[0], y[1] - y[0], z[1] - z[0]]
            v.compute_vina_maps(center, size)
            v.dock(n_poses = 20)
            energy = v.energies()[0].tolist()[0]
            self.blind_affinities.append(energy)
            v.write_poses('Blind' + self.ligand_names[jobnum] + '.pdbqt', n_poses = 1, overwrite = True) 
            
            #Strips outpose information for blind docking
            ligand_obj = p.ligand(file_path = 'Blind' + self.ligand_names[jobnum] + '.pdbqt')    
            if metric == 'Distance':
                blind_pose = ligand_obj.center
                self.blind_outposes.append(blind_pose)
            elif metric == 'Overlap':
                blind_center = ligand_obj.center
                x, y, z = [blind_center[0], blind_center[0]], [blind_center[1], blind_center[1]], [blind_center[2], blind_center[2]]
                for atom in ligand_obj.atoms:
                    if atom.x > x[1]:
                        x[1] = atom.x
                    elif atom.x < x[0]:
                        x[0] = atom.x
                    if atom.y > y[1]:
                        y[1] = atom.y
                    elif atom.y < y[0]:
                        y[0] = atom.y
                    if atom.z > z[1]:
                        z[1] = atom.z
                    elif atom.z < z[0]:
                        z[0] = atom.z
                blind_volume = [x, y, z]
                self.blind_outposes.append(blind_volume)
            elif metric == 'RMSD':  
                atoms = atom_filter(ligand_obj)
                self.blind_outposes.append(atoms)
            
            
            try:
                if metric == 'Distance':
                    true_poses = [lig.center for lig in self.protein_classes[jobnum].ligands['Primary Ligand'].iloc[1]]
                    distance = pow(pow(true_pose[0] - blind_pose[0], 2) + pow(true_pose[1] - blind_pose[1], 2) + pow(true_pose[2] - blind_pose[2], 2), .5)
                    self.blind_distance.append(distance)
         
                elif metric == 'Overlap':
                
                    def calculate_overlap_and_proportion(box1, box2):
                        # Calculate the overlap between the boxes
                        overlap = []
                        for i in range(len(box1)):
                            min_overlap = max(box1[i][0], box2[i][0])
                            max_overlap = min(box1[i][1], box2[i][1])
                            if max_overlap < min_overlap:
                                # No overlap in this dimension, return 0 proportion
                                return 0.0
                            overlap.append(max_overlap - min_overlap)
                        
                        # Calculate the volumes of the boxes
                        volume = 1
                        for i in range(len(box1)):
                            volume *= overlap[i]
                        volume1, volume2 = 1, 1
                        for i in range(len(box1)):
                            volume1 *= (box1[i][1] - box1[i][0]) 
                            volume2 *= (box2[i][1] - box2[i][0]) 
                            
                        # Calculate the proportion of overlap to the smallest volume
                        #smallest_volume = min(volume1, volume2)
                        smallest_volume = volume1
                        overlap_proportion = volume / smallest_volume
                        
                        return overlap_proportion

                    blind_overlap = 0
                    for volume in true_volumes:
                        proportion = calculate_overlap_and_proportion(volume, blind_volume)
                        if proportion > blind_overlap:
                            blind_overlap = proportion
                    self.blind_distance.append(100 * blind_overlap)

                    
                elif metric == 'RMSD':    
                    smallest_RMSD = 10000
                 
                    lig_min = min([len(pose_atoms)] + [len(atom_filter(a)) for a in self.protein_classes[jobnum].ligands['Primary Ligand'].iloc[1]])
                  
                    smallest_RMSD = 10000
                    for lig in self.protein_classes[jobnum].ligands['Primary Ligand'].iloc[1]:
                        
                        lig_atoms = atom_filter(lig)
                        squared_diff = np.square(np.array(lig_atoms)[:lig_min] - np.array(atoms)[:lig_min])
                        mean_squared_diff = np.mean(squared_diff)
                        rmsd = np.sqrt(mean_squared_diff)
                        if rmsd < smallest_RMSD:
                            smallest_RMSD = rmsd
                    self.blind_distance.append(smallest_RMSD)
            except:
                self.blind_distance.append('N/A')
          
          
    def write_results(self):
    
        try:
            os.chdir(homepath + 'UnifiedPipeline/Unified_Results')
        except:
            os.mkdir(homepath + 'UnifiedPipeline/Unified_Results')
            os.chdir(homepath + 'UnifiedPipeline/Unified_Results')
         
        if self.metric == None:  
            results = pd.DataFrame.from_dict({'Ligands':self.ligands, 'Proteins':self.proteins, 'Blind Affinity':self.blind_affinities})
        elif self.metric == 'Distance':
            results = pd.DataFrame.from_dict({'Ligands':self.ligands, 'Proteins':self.proteins, 'Blind Distances': self.blind_distance, 'Blind Affinity':self.blind_affinities})
        elif self.metric.upper() == 'RMSD':
             results = pd.DataFrame.from_dict({'Ligands':self.ligands, 'Proteins':self.proteins, 'Blind RMSD': self.blind_distance, 'Blind Affinity':self.blind_affinities})
        elif self.metric == 'Overlap':
             results = pd.DataFrame.from_dict({'Ligands':self.ligands, 'Proteins':self.proteins, 'Blind Overlap': self.blind_distance, 'Blind Affinity':self.blind_affinities})
            
        for key in results.columns:
            print(f'{key}: {results[key]}')
        results.to_csv(str(self.instance) + 'Results.csv')
        
        
if __name__ == '__main__':
    pipe = run_pipeline()
    pipe.generate_files()
    pipe.atomize_files()
    pipe.run_docking()
    pipe.write_results()
