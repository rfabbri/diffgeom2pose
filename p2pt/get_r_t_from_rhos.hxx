#ifndef get_r_t_from_rhos_hxx_
#define get_r_t_from_rhos_hxx_

namespace P2Pt {

template<typename T>
void
pose_poly<T>::
get_r_t_from_rhos(
	const int ts_len,
	const T (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN], const int (&sigmas1_len)[TS_MAX_LEN],
	const T (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN], const int (&sigmas2_len)[TS_MAX_LEN],
	const T (&rhos1)[ROOT_IDS_LEN], const T (&rhos2)[ROOT_IDS_LEN],
	const T (&gama1)[3], const T (&tgt1)[3],
	const T (&gama2)[3], const T (&tgt2)[3],
	const T (&Gama1)[3], const T (&Tgt1)[3],
	const T (&Gama2)[3], const T (&Tgt2)[3],
	T (*output)[RT_MAX_LEN][4][3], int *output_len
)
{
	T lambdas1[TS_MAX_LEN][TS_MAX_LEN];
	T lambdas2[TS_MAX_LEN][TS_MAX_LEN];
	T Gama_sub[3];
  T den1[3], den2[3], den3[3]; // Buffers

	for (int i = 0; i < ts_len; i++) {
		for (int j = 0; j < sigmas1_len[i]; j++) {
			common::vec1vec2_3el_sub(Gama1, Gama2, Gama_sub);            // Gama_sub = Gama1 - Gama2

			common::vec_3el_mult_by_scalar(rhos1[i], gama1, den1);       // den1 = rhos1(i) * gama1
			common::vec_3el_mult_by_scalar(rhos2[i], gama2, den2);       // den2 = rhos2(i) * gama2
			common::vec1vec2_3el_sub(den1, den2, den1);                  // den1 = den1 - den2

			common::vec_3el_mult_by_scalar(rhos1[i], tgt1, den2);        // den2 = rhos1(i) * tgt1
			common::vec_3el_mult_by_scalar(sigmas1[i][j], gama1, den3);  // den3 = sigmas1{i}(j) * gama1
			common::vec1vec2_3el_sum(den2, den3, den2);                  // den2 = den2 + den3

			// lambdas1{i}(j) = Gama_sub' * Tgt1 / den1' * den2
			lambdas1[i][j] = common::vec1vec2_3el_dot(Gama_sub, Tgt1) / common::vec1vec2_3el_dot(den1, den2);
		}
		for (int j = 0; j < sigmas1_len[i]; j++) { // XXX sigmas2_len ?? bug ??
			common::vec1vec2_3el_sub(Gama1, Gama2, Gama_sub);            // Gama_sub = Gama1 - Gama2

			common::vec_3el_mult_by_scalar(rhos1[i], gama1, den1);       // den1 = rhos1(i) * gama1
			common::vec_3el_mult_by_scalar(rhos2[i], gama2, den2);       // den2 = rhos2(i) * gama2
			common::vec1vec2_3el_sub(den1, den2, den1);                  // den1 = den1 - den2

			common::vec_3el_mult_by_scalar(rhos2[i], tgt2, den2);        // den2 = rhos2(i) * tgt2
			common::vec_3el_mult_by_scalar(sigmas2[i][j], gama2, den3);  // den3 = sigmas2{i}(j) * gama2
			common::vec1vec2_3el_sum(den2, den3, den2);                  // den2 = den2 + den3

			// lambdas2{i}(j) = Gama_sub' * Tgt2 / den1' * den2
			lambdas2[i][j] = common::vec1vec2_3el_dot(Gama_sub, Tgt2) / common::vec1vec2_3el_dot(den1, den2);
		}
	}

	//% Rotation:

	//% RA = B = > R = B / A
	//% TODO: use svd or some other way to be sure R is unique and orthogonal.
	//% right now just testing det(R) = 1

	const T A[3][3] = {
		Gama_sub[0], Tgt1[0], Tgt2[0],
		Gama_sub[1], Tgt1[1], Tgt2[1],
		Gama_sub[2], Tgt1[2], Tgt2[2],
	};

	T inv_A[3][3]; common::invm3x3(A, inv_A);

	// Matrix containing Rotations and Translations
	T (&RT)[RT_MAX_LEN][4][3] = *output;
	int &RT_len               = *output_len;

	RT_len = 0;
	for (int i = 0; i < ts_len; i++) {
		for (int j = 0; j < sigmas1_len[i]; j++, RT_len++) {
			T (&Rots)[4][3] = RT[RT_len];
			T (&Transls)[3] = RT[RT_len][3];

			#define B_row(r) \
				rhos1[i]*gama1[(r)] - rhos2[i]*gama2[(r)], \
				lambdas1[i][j]*(rhos1[i]*tgt1[(r)] + sigmas1[i][j]*gama1[(r)]), \
				lambdas2[i][j]*(rhos2[i]*tgt2[(r)] + sigmas2[i][j]*gama2[(r)])

			const T B[3][3] = {
				B_row(0),
				B_row(1),
				B_row(2)
			};

			common::multm3x3(B, inv_A, Rots);

			T buff1[3], buff2[3];

			// Transls{end+1} = rhos1(i)*gama1 - Rots{end}*Gama1;
			common::vec_3el_mult_by_scalar(rhos1[i], gama1, buff1);
			common::multm_3x3_3x1(Rots, Gama1, buff2);
			common::vec1vec2_3el_sub(buff1, buff2, Transls);
		}
	}
}
}

#endif // !get_r_t_from_rhos_hxx_

