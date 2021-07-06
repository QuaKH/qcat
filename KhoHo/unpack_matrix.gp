/* Given the index (datapos) of a diagram whose differentials have already been computed, and an i,j,
 * returns the i,j-th differential matrix.
 */
matr_populate(datapos, i, j) =

{
	local(entry, mat);



	mat=

			matrix(chain_ranks[datapos][j, i + 1], chain_ranks[datapos][j, i]);

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

