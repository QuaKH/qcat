#!/bin/tcsh
 
#SBATCH --job-name=qcat
#SBATCH -n 20
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=jayr2@illinois.edu
 
bash run_all__pd_codes.sh pd_code_files/pretzel_pd_codes.txt db/eigs.db Pretzel