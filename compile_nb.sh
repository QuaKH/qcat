jupyter nbconvert --to python get_eigs.ipynb

echo "if __name__ == \"__main__\":
    import sys
    get_knot_eigs(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))" >> get_eigs.py