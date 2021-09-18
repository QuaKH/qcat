# Arguments:
# 1: pd code input file path

echo "cleaning up"
rm gp_pd_code_input
rm run_all__pd_codes_TIMES
rm -r KhoHo/differentials/*
rm -r eigs/*
rm -r laplacian_sparsity/*

# replace directories
mkdir KhoHo/differentials
mkdir eigs
mkdir laplacian_sparsity

echo "compiling notebooks"
./compile_nb.sh

# parse pd code input file and write output gp commands to file
python3 parse_pd_code.py $1 > gp_pd_code_input

# iterate through each command in the file; time each separately
index=0
while read p; do
    crossings=(`echo $p | grep -Po '], \K[^,]*'`)
    dir="KhoHo/differentials/knot_${crossings}_${index}"
    mkdir $dir

    echo "KNOT $crossings $index" >> run_all__pd_codes_TIMES
    echo "DIFFERENTIALS:" >> run_all__pd_codes_TIMES

    echo Computing differentials...
    cd KhoHo
    # { time echo "$p" | gp -f -q -s 1000000000 KH unpack_matrix.gp > ../garbage; } 2>> ../run_all__pd_codes_TIMES
    echo "$p" | /usr/bin/time -o ../run_all__pd_codes_TIMES -a --format='%Uuser %Ssystem %Eelapsed %PCPU %MmaxKB %tavgKB %Wswaps %ccontext_switch %wwaits' gp -f -q -s 1000000000 KH unpack_matrix.gp > /dev/null;
    cd ..

    echo "EIGENVALUES:" >> run_all__pd_codes_TIMES

    echo Getting eigenvalues...
    /usr/bin/time -o run_all__pd_codes_TIMES -a --format='%Uuser %Ssystem %Eelapsed %PCPU %MmaxKB %tavgKB %Wswaps %ccontext_switch %wwaits' python3 get_eigs.py ${dir} ${crossings} ${index}

    echo Deleting differentials...
    rm -r $dir

    printf "\n" >> run_all__pd_codes_TIMES
    index=$(($index + 1))
done < gp_pd_code_input

# clean up
rm gp_pd_code_input
