#!/bin/tcsh
 
#SBATCH --job-name=qcat_cuda
#SBATCH -n 100
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=anshula2@illinois.edu
#SBATCH --partition="gpu"

bash run_all__pd_codes.sh snappy_out_12.txt