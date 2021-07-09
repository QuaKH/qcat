arr=(1 1 2 3 7 21 49 165)

for i in 3 4 5 6 7 8 9 10;
do
    j=1
    while [ $j -le ${arr[$i-3]} ]
    do
        echo "compute_knot_differential(${i}, ${j})" >> gp_out
        
        # ./run_get_eigs.sh $i $j
        # echo $i $j

        dir="./KhoHo/differentials/knot_${i}_${j}"
        mkdir $dir

        j=$(($j + 1))
    done
done

cd KhoHo
cat gp_out | gp -s 120000000 KH unpack_matrix.gp

for i in 3 4 5 6 7 8 9 10;
do
    j=1
    while [ $j -le ${arr[$i-3]} ]
    do
        python3 -c "import get_eigs; get_eigs.get_knot_eigs(\"${dir}\", ${i}, ${j})"
        j=$(($j + 1))
    done
done
