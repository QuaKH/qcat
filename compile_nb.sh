# compile get_eigs notebook
jupyter nbconvert --to python get_eigs.ipynb
echo "if __name__ == \"__main__\":
    import sys
    get_knot_eigs(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))" >> get_eigs.py

# compile parse_pd_code notebook
jupyter nbconvert --to python parse_pd_code.ipynb
echo "if __name__ == \"__main__\":
    import sys
    parse_pd_code(sys.argv[1])" >> parse_pd_code.py