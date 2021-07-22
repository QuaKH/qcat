# clean previous data files
./clean_all.sh
# compile python notebook into script
./compile_nb.sh

arr=(1 1 2 3 7 21 49 165)

for i in 3 4 5 6 7 8 9 10;
do
    j=1
    while [ $j -le ${arr[$i-3]} ]
    do
        # create dir for knot differentials if it doesn't exist
        # TODO: fix dir creation for general pd code iteration
        dir="./KhoHo/differentials/knot_${i}_${j}"
        mkdir $dir

        # add compute differential call to gp_out
        echo "compute_knot_differential(${i}, ${j})" >> gp_out

        j=$(($j + 1))
    done
done

python3 parse_pd_code.py "./pd_code_input.txt"

# compute all knot differentials
cd KhoHo
echo "$(cat ../gp_out)" | gp -f -q -s 120000000 KH unpack_matrix.gp

# head -166 gp_out > gp_out_head
# head -83 gp_out_head > gp_out1
# tail -83 gp_out_head > gp_out2
# tail -83 gp_out > gp_out3

# echo "$(cat ../gp_out1)" | gp -f -q -s 120000000 KH unpack_matrix.gp
# echo "$(cat ../gp_out2)" | gp -f -q -s 120000000 KH unpack_matrix.gp
# echo "$(cat ../gp_out3)" | gp -f -q -s 120000000 KH unpack_matrix.gp

# parallel --jobs 6 gp -f -q -s 120000000 KH unpack_matrix.gp :::: gp_out

echo "finished parallel part"
cd ..

# compute all eigenvalues
for crossings in 3 4 5 6 7 8 9 10;
do
    index=1
    while [ $index -le ${arr[$crossings-3]} ]
    do
        dir="./KhoHo/differentials/knot_${crossings}_${index}"
        python3 get_eigs.py ${dir} ${crossings} ${index}
        echo ${crossings} ${index}

        # remove all the differentials
        rm -r "./KhoHo/differentials/knot_${crossings}_${index}"

        index=$(($index + 1))
    done
done
