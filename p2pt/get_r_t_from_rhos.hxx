// % to be called from pose_from_point_tangents_root_find_function_any.m
#include "common.hxx"

namespace P2Pt {

template<typename T>
void
pose_poly<T>::
get_r_t_from_rhos(
	const int ts_len,
	const T (&sigmas1)[TS_LEN][TS_LEN], const int (&sigmas1_end)[TS_LEN],
	const T (&sigmas2)[TS_LEN][TS_LEN], const int (&sigmas2_end)[TS_LEN],
	const T (&rhos1)[TS_LEN], const T (&rhos2)[TS_LEN],
	const T (&gama1)[3], const T (&tgt1)[3],
	const T (&gama2)[3], const T (&tgt2)[3],
	const T (&Gama1)[3], const T (&Tgt1)[3],
	const T (&Gama2)[3], const T (&Tgt2)[3],
	T (*output)[RT_LEN][4][3]
)
{
	//% to be called from pose_from_point_tangents_root_find_function_any.m

	//% Lambdas:
	static T lambdas1[TS_LEN][TS_LEN] = { {0} };
	static T lambdas2[TS_LEN][TS_LEN] = { {0} };

	static T Gama_sub[3];

	for (int i = 0; i < ts_len; i++) {
		//if (sigmas1_end[i] == 0) {
			// Is this necessary?
			//lambdas1[i] = {};
			//lambdas2[i] = {};
		//}
		//lambdas1[i] = zeros(size(sigmas1{i}));
		//lambdas2[i] = zeros(size(sigmas2{i}));

		// Buffers
		static T den1[3], den2[3], den3[3];

		// TODO: Check if separate loops are better for vectorization
		for (int j = 0; j < sigmas1_end[i]; j++) {
			common::vec1vec2_sub3(Gama1, Gama2, Gama_sub);            // Gama_sub = Gama1 - Gama2

			common::vec_mult_by_scalar3(rhos1[i], gama1, den1);       // den1 = rhos1(i) * gama1
			common::vec_mult_by_scalar3(rhos2[i], gama2, den2);       // den2 = rhos2(i) * gama2
			common::vec1vec2_sub3(den1, den2, den1);                  // den1 = den1 - den2

			common::vec_mult_by_scalar3(rhos1[i], tgt1, den2);        // den2 = rhos1(i) * tgt1
			common::vec_mult_by_scalar3(sigmas1[i][j], gama1, den3);  // den3 = sigmas1{i}(j) * gama1
			common::vec1vec2_sum3(den2, den3, den2);                  // den2 = den2 + den3

			// lambdas1{i}(j) = Gama_sub' * Tgt1 / den1' * den2
			lambdas1[i][j] = common::vec1vec2_dot3(Gama_sub, Tgt1) / common::vec1vec2_dot3(den1, den2);
		}
		for (int j = 0; j < sigmas1_end[i]; j++) {
			common::vec1vec2_sub3(Gama1, Gama2, Gama_sub);            // Gama_sub = Gama1 - Gama2

			common::vec_mult_by_scalar3(rhos1[i], gama1, den1);       // den1 = rhos1(i) * gama1
			common::vec_mult_by_scalar3(rhos2[i], gama2, den2);       // den2 = rhos2(i) * gama2
			common::vec1vec2_sub3(den1, den2, den1);                  // den1 = den1 - den2

			common::vec_mult_by_scalar3(rhos2[i], tgt2, den2);        // den2 = rhos2(i) * tgt2
			common::vec_mult_by_scalar3(sigmas2[i][j], gama2, den3);  // den3 = sigmas2{i}(j) * gama2
			common::vec1vec2_sum3(den2, den3, den2);                  // den2 = den2 + den3

			// lambdas2{i}(j) = Gama_sub' * Tgt2 / den1' * den2
			lambdas2[i][j] = common::vec1vec2_dot3(Gama_sub, Tgt2) / common::vec1vec2_dot3(den1, den2);
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
	T (&RT)[7][4][3] = *output;
	int RT_end = 0;

	for (int i = 0; i < ts_len; i++) {
		for (int j = 0; j < sigmas1_end[i]; j++, RT_end++) {
			T (&Rots)[4][3] = RT[RT_end];
			T (&Transls)[3] = RT[RT_end][3];

			#define B_row(r) \
				rhos1[i]*gama1[r] - rhos2[i]*gama2[r], \
				lambdas1[i][j]*(rhos1[i]*tgt1[r] + sigmas1[i][j]*gama1[r]), \
				lambdas2[i][j]*(rhos2[i]*tgt2[r] + sigmas2[i][j]*gama2[r]) \

			const T B[3][3] = {
				B_row(0),
				B_row(1),
				B_row(2)
			};

			common::multm3x3(B, inv_A, Rots);

			//% be sure it's close to 1:
			//common::det3x3(Rots[Rots_end]);

			T buff1[3], buff2[3];

			// Transls{end+1} = rhos1(i)*gama1 - Rots{end}*Gama1;
			common::vec_mult_by_scalar3(rhos1[i], gama1, buff1);
			common::multm_3x3_3x1(Rots, Gama1, buff2);
			common::vec1vec2_sub3(buff1, buff2, Transls);

			//% this should be the same
			// rhos2(i)*gama2 - Rots{end}*Gama2;
			//common::vec_mult_by_scalar(rhos2[i], gama2, buff1);
			//common::multm_3x3_3x1(Rots[Rots_end], Gama2, buff2);
			//common::vec1vec2_sub3(buff1, buff2, buff1);

			//for (int i = 0; i < 3; i++) {
			//	assert(Transls[Transls_end - 1][i] - buff1[i] < 1.0e-4);
			//}
		}
	}

	//std::copy(std::begin(Rots[0][0][0]), std::end(Rots[Rots_end][2][2]), std::begin(RT[0][0][0]));
	//std::copy(std::begin(Transls[0][0]), std::end(Transls[Transls_end][2]), std::begin(RT[3][0][0]));
}

}

