#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define output directory parameters
def baseDir = '/home/ytf/peptide_MDI/peptide'

// Defining Channels
def configChannel    = Channel.fromPath("${baseDir}/config.yaml")
def inputFileChannel = Channel.fromPath("${baseDir}/in_put/Wheat_protein.fasta")
def jsonFileChannel  = Channel.fromPath("${baseDir}/DeepDigest")
def Autodock_ProteinChannel  = Channel.fromPath("${baseDir}/in_put/protein.pdbqt")
def plot_Channel  = Channel.fromPath("${baseDir}/plot")
def gpfFileChannel = Channel.fromPath("${baseDir}/in_put/config.gpf")

workflow {
    def step1_ch = Step1(configChannel, inputFileChannel, jsonFileChannel)
    
    def step2_ch = Step2(configChannel, step1_ch)

    def step3_ch = Step3(configChannel, step2_ch)

    def step4_ch = Step4(configChannel, step3_ch)

    def step5_ch = Step5(configChannel, Autodock_ProteinChannel, gpfFileChannel, step4_ch)

    def step6_ch = Step6(configChannel,step5_ch,Autodock_ProteinChannel)

    def step7_ch = Step7(configChannel,step6_ch)

    def step8_ch = Step8(configChannel,step6_ch)

    def step9_ch = Step9(configChannel,step6_ch)

    def step10_ch = Step10(configChannel,step5_ch)

    Step11(configChannel,plot_Channel,step8_ch,step9_ch,step10_ch)
}

process Step1 {
    conda 'env1.yml'
    
    tag 'Predicting Enzymatic Peptides'

    input:
      path configFile, from: configChannel
      path inputFile,  from: inputFileChannel
      path jsonDir,    from: jsonFileChannel

    output:
      path "out_put/Step_1.txt", emit: step1_result

    publishDir "${baseDir}/", mode: 'copy'

    script:
    """
    echo "Config file: ${configFile}"
    echo "Input file: ${inputFile}"
    echo "JSON directory: ${jsonDir}"
    echo "Publishing results to ${baseDir}/out_put"
    
    python ${baseDir}/Step1_Predicting_Enzymatic_Peptides.py ${configFile} ${jsonDir}
    """
}

process Step2 {
    conda 'env1.yml'
    
    tag "Peptide Data Cleaning"

    input:
      path configFile
      path step1_file

    output:
      path "out_put/Step_2.txt", emit: step2_result

    publishDir "${baseDir}/", mode: 'copy'

    script:
    """
    python ${baseDir}/Step2_Peptide_Data_Cleaning.py ${configFile} ${step1_file}
    """
}

process Step3 {
    conda 'env1.yml'
    
    tag "Creating PDB files"

    input:
      path configFile
      path step2_file

    output:
      path "out_put/step3_results/", emit: step3_pdbs 

    publishDir "${baseDir}/", mode: 'copy'

    script:
    """
    python ${baseDir}/Step3_Creating_PDB_files.py ${configFile} ${step2_file}
    """
}


process Step4 {
    conda 'env1.yml'
    
    tag "Switch from PDB to PDBQT"

    input:
      path configFile
      path step3_pdbs

    output:
      path "out_put/step4_results/", emit: step4_results

    publishDir "${baseDir}/", mode: 'copy'

    errorStrategy 'ignore' 

    script:
    """
    python ${baseDir}/Step4_Switch_from_PDB_to_PDBQT.py ${configFile} ${step3_pdbs}
    """
}



process Step5 {
    conda 'env4.yml'
    
    tag "Molecular docking"

    input:
      path configFile, from: configChannel
      path Protein,  from: Autodock_ProteinChannel
      path gpf,    from: gpfFileChannel
      path step4_results

    output:
      path "out_put/step5_results/", emit: step5_results

    publishDir "${baseDir}", mode: 'copy'

    script:
    """
    python ${baseDir}/Step5_Moleccular_docking_GPU.py ${configFile} ${Protein} ${gpf} ${step4_results}
    """
}


process Step6 {
    conda 'env1.yml'
    
    tag "Merging molecular docking results"
    input:
      path configFile, from: configChannel
      path Protein,  from: Autodock_ProteinChannel
      path step5_results

    output:
      path "out_put/step6_results/", emit: step6_results

    publishDir "${baseDir}", mode: 'copy'

    script:
    """
    python ${baseDir}/Step6_Merging_molecular_docking_results.py ${configFile} ${step5_results} ${Protein} 
    """
}

process Step7 {
    conda 'env1.yml'
    
    tag "Visualization of molecular docking results"

    input:
      path configFile, from: configChannel
      path step6_results

    output:
      path "out_put/step7_results/", emit: step7_results

    publishDir "${baseDir}", mode: 'copy'

    script:
    """
    python ${baseDir}/Step7_Visualization_of_molecular_docking_results.py ${configFile} ${step6_results} 
    """
}

process Step8 {
    conda 'env1.yml'
    
    tag "Surface scoring"

    input:
      path configFile, from: configChannel
      path step6_results

    output:
      path "out_put/step8_results/", emit: step8_results

    publishDir "${baseDir}", mode: 'copy'

    script:
    """
    python ${baseDir}/Step8_Surface_scoring.py ${configFile} ${step6_results} 
    """
}


process Step9 {
    conda 'env2.yml'
    
    tag "Affinity prediction"
    input:
      path configFile, from: configChannel
      path step6_results

    output:
      path "out_put/step9_results/", emit: step9_results

    publishDir "${baseDir}", mode: 'copy'

    script:
    """
    python ${baseDir}/Step9_Affinity_prediction.py ${configFile} ${step6_results} 
    """
}


process Step10 {
    conda 'env1.yml'
    
    tag "Extraction of molecular docking binding energy"

    input:
      path configFile, from: configChannel
      path step5_results

    output:
      path "out_put/step10_results/", emit: step10_results

    publishDir "${baseDir}", mode: 'copy'

    script:
    """
    python ${baseDir}/Step10_Extraction_of_molecular_docking_binding_energy.py ${configFile} ${step5_results} 
    """
}

process Step11 {
    conda 'env3.yml'
    
    tag "plot"

    input:
      path configFile, from: configChannel
      path Rscript,from: plot_Channel
      path step8_results
      path step9_results
      path step10_results

    output:
      path "out_put/plot_results/", emit: step11_results

    publishDir "${baseDir}", mode: 'copy'

    script:
    """
    python ${baseDir}/Step11_plot.py ${configFile} ${Rscript} ${step8_results} ${step9_results}${step10_results}
    """
}


