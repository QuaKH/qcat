/* Given the index (datapos) of a diagram whose differentials have already been computed, and an i,j,   
 * returns the i,j-th differential matrix. i and j are the indices of the matrices in the big matrix of
 * differentials. NOT the i and j of the grading
 */
get_differential_m(datapos, i, j) =

{
	local(entry, mat);

	mat=matrix(chain_ranks[datapos][j, i + 1], chain_ranks[datapos][j, i]);

	for (pos = 1, allmatr_length[datapos][i, j],

		if (is_arch_64,

			entry = allmatr[datapos][i, j][pos];

			mat[abs(entry) \ arch64_mask, abs(entry) % arch64_mask] = sign(entry);

		,

			mat[allmatr[datapos][i, j][2*pos - 1], 

				    abs(allmatr[datapos][i, j][2 * pos])] = 

					sign(allmatr[datapos][i, j][2 * pos]);

		);

	);
	return(mat);

}

/* Given the index (datapos) of a diagram whose differentials have already been computed, and an i,j,   
 * returns the i,j-th differential matrix.
 */
get_differential(datapos, i, j) = 
{
	return(get_differential_m(datapos, i2m(datapos, i), j2m(datapos, j)))
}



/**
 * write differential matrix to binary file in sparse format
 */
get_sparse_differential(datapos,i,j) = {

	/*
		local(entry, row, col, val);
	*/

	i = i2m(datapos, i);
	j = j2m(datapos, j);

	/*
		TODO: fix convention; add negatives
	*/
	write(Str("sparse_d_",datapos,"_",i,"_",j), allmatr[datapos][i, j]);

	/*
		for (pos = 1, allmatr_length[datapos][i, j],

			entry = allmatr[datapos][i, j][pos];

		
			row = abs(entry) \ arch64_mask;
			col = abs(entry) % arch64_mask;
			val = sign(entry);
			str = Str(row, " ", col, " ", val)

			write("get_sparse_differential_output", str);
		);
	*/

	return 1;
}

