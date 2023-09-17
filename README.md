# AutoDocker
A powerful, lightweight wrapper for blind autodocking with Vina, which allows for efficient docking of large-scale datasets 

# Environment
Conda is required to adequately prepare the protein files. Once the environment is created, the requirements.txt file can be installed within the conda environment. Finally, install mgltools with

    conda install -c bioconda mgltools

And the environment is prepared!

# Requirements and Options
The primary function is blind.py, located in the UnifiedPipeline folder. The file can be run with the following flags:

    -t, --total_jobs (integer) = Total amount of jobs per input file, default = 1
   
    -i, --instance (integer) = Current instance of run (i.e. if there is a file with 100 protein/ligand pairs, -t of 50 and -i of 25 would mean that pair 49 and 50 woud be run from the file), default = 1
   
    -f, --file (string) = Input file with formatted columns. This flag is required.
   
    -p, --protein (string) = Input protein ID (PDB or Uniprot), default = PDB
   
    -l, --ligand (string) = Input ligand ID (SMILE/smi or InChiKey/ick), default = SMILE

    
The input data file must contain two columns, one with protein IDs (named Proteins), and one with ligand IDs (named Ligands). Make sure to place the input file in the same directory as blind.py
Upon running the pipeline, a folder named Unified_Results will be created that contains the results from the given run, with the format Results<job instance>.csv 

# Batch Runs with Slurm
A sample batch script is located in the batchRuns folder. This script can be ran, with the format:

    sbatch Unified.bash <total jobs> <starting job, usually 0> <file name>

With a slurm manager present, docking jobs will be run in parallel until completion.
  
