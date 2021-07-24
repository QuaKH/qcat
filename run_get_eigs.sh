dir="./KhoHo/differentials/knot_${1}_${2}"
mkdir $dir
cd KhoHo
echo "compute_knot_differential(${1}, ${2})" | pari/gp -s 120000000 KH unpack_matrix.gp
cd ..
python3 -c "import get_eigs; get_eigs.get_knot_eigs(\"${dir}\", ${1}, ${2})"
