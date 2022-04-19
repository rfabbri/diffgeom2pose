// % to be called from rf_pose_from_point_tangents_root_find_function_any.m
#include "common.hxx"

template<typename T>
void
rf_get_r_from_rhos(
	const int ts_len,
	T sigmas1[ts_len][ts_len], int end_sigmas1[ts_len],
	T sigmas2[ts_len][ts_len], int end_sigmas2[ts_len],
	T rhos1[t_vector_len], T rhos2[t_vector_len]
)
{
	//% to be called from rf_pose_from_point_tangents_root_find_function_any.m

	//% Lambdas:
	static T lambdas1[ts_len][ts_len] = {0};
	static T lambdas2[ts_len][ts_len] = {0};

	static T DGama[3];

	for (int i = 0; i < ts_len; i++) {
		if (end_sigmas1[i] == 0) {
			// Is this necessary?
			//lambdas1[i] = {};
			//lambdas2[i] = {};
		}
		//lambdas1[i] = zeros(size(sigmas1{i}));
		//lambdas2[i] = zeros(size(sigmas2{i}));

		// Buffers
		static T den1[3], den2[3], den3[3];

		// TODO: Check if separate loops are better for vectorization
		for (int j = 0; end_sigmas1[i]; j++) {
			common::vec1vec2_sub(Gama1, Gama2, DGama);               // DGama = Gama1 - Gama2

			// TODO: Modify to product of 3 elements
			common::vec_mult_by_scalar3(rhos1[i], gama1, den1);      // den1 = rhos1(i) * gama1
			common::vec_mult_by_scalar3(rhos2[i], gama2, den2);      // den2 = rhos2(i) * gama2
			common::vec1vec2_sub(den1, den2, den1);                  // den1 = den1 - den2

			common::vec_mult_by_scalar3(rhos1[i], tgt1, den2);       // den2 = rhos1(i) * tgt1
			common::vec_mult_by_scalar3(sigmas1[i][j], gama1, den3); // den3 = sigmas1{i}(j) * gama1
			common::vec1vec2_sum(den2, den3, den2);                  // den2 = den2 + den3

			// lambdas1{i}(j) = DGama' * Tgt1 / den2' * den3
			lambdas1[i][j] = common::vec1vec2_dot(DGama, Tgt1) / common::vec1vec2_dot(den2, den3);
		}
		for (int j = 0; end_sigmas1[i]; j++) {
			common::vec1vec2_sub(Gama1, Gama2, DGama);               // DGama = Gama1 - Gama2

			// TODO: Modify to product of 3 elements
			common::vec_mult_by_scalar3(rhos1[i], gama1, den1);      // den1 = rhos1(i) * gama1
			common::vec_mult_by_scalar3(rhos2[i], gama2, den2);      // den2 = rhos2(i) * gama2
			common::vec1vec2_sub(den1, den2, den1);                  // den1 = den1 - den2

			common::vec_mult_by_scalar3(rhos2[i], tgt2, den2);       // den2 = rhos2(i) * tgt2
			common::vec_mult_by_scalar3(sigmas2[i][j], gama2, den3); // den3 = sigmas2{i}(j) * gama2
			common::vec1vec2_sum(den2, den3, den2);                  // den2 = den2 + den3

			// lambdas2{i}(j) = DGama' * Tgt2 / den2' * den3
			lambdas2[i][j] = common::vec1vec2_dot(DGama, Tgt2) / common::vec1vec2_dot(den2, den3);
		}
	}

	//% Rotation:

	//% RA = B = > R = B / A
	//% TODO: use svd or some other way to be sure R is unique and orthogonal.
	//% right now just testing det(R) = 1

	T A[3][3] = {
		DGama[0], Tgt1[0], Tgt2[0],
		DGama[1], Tgt1[1], Tgt2[1],
		DGama[2], Tgt1[2], Tgt2[2],
	};
	T B[3][3];

	T inv_A[3][3]; common::invm3x3(A, inv_A);

	#define TEMP 10;

	T Rots[TEMP][3][3];
	T Transls[TEMP][3][3];
	int Rots_end    = 0;
	int Transls_end = 0;

	for (int i = 0; i < ts_len; i++) {
		for (int j = 0; j < end_sigmas1[0]; j++) {

			#define B_row(r) \
				rhos[i]*gama1[r] - rhos[i]*gama2[r], \
				lambdas1[i][j]*(rhos1[i]*tgt1[r] + sigmas1[i][j]*gama1[r]), \
				lambdas2[i][j]*(rhos2[i]*tgt2[r] + sigmas2[i][j]*gama2[r]) \

			T B[3][3] = {
				B_row(0),
				B_row(1),
				B_row(2)
			};

			// TODO: Check if `A` and `B` are always 3x3
			common::multm3x3(B, inv_A, Rots[Rots_end++]);

			//% be sure it's close to 1:
			//common::det3x3(Rots[Rots_end]);

			T buff1[3];
			T buff2[3];

			// Transls{end+1} = rhos1(i)*gama1 - Rots{end}*Gama1;
			common::vec_mult_by_scalar(rhos1[i], gama1, buff1);
			common::multm_3x3_3x1(Rots[Rots_end], Gama1, buff2);
			common::vec1vec2_sub(buff1, buff2, Transls[Transls_end++]);

			//% this should be the same
			// rhos2(i)*gama2 - Rots{end}*Gama2;
			//common::vec_mult_by_scalar(rhos2[i], gama2, buff1);
			//common::multm_3x3_3x1(Rots[Rots_end], Gama2, buff2);
			//common::vec1vec2_sub(buff1, buff2, buff1);

			//for (int i = 0; i < 3; i++) {
			//	assert(Transls[Transls_end - 1][i] - buff1[i] < 1.0e-4);
			//}
		}
	}


}

