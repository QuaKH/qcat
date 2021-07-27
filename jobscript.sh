#!/bin/tcsh
 
#SBATCH --job-name=qcat
#SBATCH -n 4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4096
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=sdontha2@illinois.edu
 
bash run_all__pd_codes.sh final_snappy_out_4_to_12.txt