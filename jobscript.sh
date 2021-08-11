#!/bin/tcsh
 
#SBATCH --job-name=qcat
#SBATCH -n 20
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=sdontha2@illinois.edu
 
bash run_all__pd_codes.sh snappy_out_12_to_20.txt