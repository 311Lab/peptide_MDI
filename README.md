# peptide_MDI

### Authors
- Tianfei YU  
- Tianshuo HU  
- Yan LIU  
- Qian LIN  

---

## Introduction
**peptide_MDI** is a workflow designed for the prediction and screening of bioactive peptides from proteins. This streamlined process integrates computational tools to analyze protein sequences and identify potential peptides with biological activity.

---

## Workflow
### Step 1: Configure the Environment
Before running the workflow, ensure you have **Anaconda** installed on your system.

1. Create a new Conda environment:
   ```bash
   conda create -n peptide_MDI

2. Activate the environment：
   ```bash
   conda activate peptide_MDI
   
3. Install Nextflow：
   ```bash
   conda install -c bioconda nextflow

### Step 2: Navigate to the Working Directory
Switch to the directory containing the workflow scripts and resources:
```bash
   cd /path/to/peptide_MDI
```

### Step 3: Upload the Required Files into the `in_put` Folder
Ensure the following files are placed in the `in_put` folder before running the workflow:

- **Protein sequence file**:  
  `Protein.fasta`  
  Contains the protein sequences to be analyzed.

- **Receptor structure file**:  
  `Receptor.pdbqt`  
  The 3D structure of the receptor used for docking studies.

- **Docking configuration file**:  
  `Docking_config.txt`  
  Provides information about the active pocket of the receptor.

---

### Step 4: Run the Workflow
Execute the workflow using Nextflow by running the following command:
```bash
nextflow run main.nf -with-conda 


   
