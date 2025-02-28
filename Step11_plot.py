import os
import subprocess

# Defining the R script folder path
r_script_folder = "./plot"

# Get all R script file names
r_scripts = [
    "Plot_1_Docking.R",
    "Plot_2_Interface.R",
    "Plot_3_Prodigy.R"
]

# Execute each R script in turn
for r_script in r_scripts:
    r_script_path = os.path.join(r_script_folder, r_script)
    print(f"Executing: {r_script_path}")
    
    try:
        # Use subprocess to invoke the Rscript command to execute an R script.
        subprocess.run(
            ["Rscript", r_script_path],
            check=True
        )
        print(f"Successful execution: {r_script}")
    except subprocess.CalledProcessError as e:
        print(f"failure of execution: {r_script}")
        print(f"error message: {e}")
