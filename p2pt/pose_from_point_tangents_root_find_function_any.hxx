#include "common.hxx"
#include "poly.h"

#include "pose_from_point_tangents_2.hxx"
#include "find_bounded_root_intervals.hxx"
#include "get_r_t_from_rhos.hxx"

namespace P2Pt {

template<typename T>
void pose_from_point_tangents_root_find_function_any(
	const T (&gama1)[3], const T (&tgt1)[3],
	const T (&gama2)[3], const T (&tgt2)[3],
	const T (&Gama1)[3], const T (&Tgt1)[3],
	const T (&Gama2)[3], const T (&Tgt2)[3],
	T (*output)[2][RT_LEN][4][3]
)
{
	// % This is the main routine to find roots. Can be used with any input.


	// % test for geometric degeneracy -------------------------------

	static T DGama[3];
	common::vec1vec2_3el_sub(Gama1, Gama2, DGama);
	common::vec_3el_div_by_scalar(common::norm(DGama, 3), DGama, DGama);

	// Matrix for degeneracy calculation
	T &degen = (*output)[1][0][0][0];
	const T degen_matrix[3][3] = {
		DGama[0], Tgt1[0], Tgt2[0],
		DGama[1], Tgt1[1], Tgt2[1],
		DGama[2], Tgt1[2], Tgt2[2]
	};
	degen = common::det3x3(degen_matrix);

	if (std::abs(degen) < 1.0e-3) {
		std::cout << "data point not reliable" << std::endl;
		//(*output)[0] = nullptr;
		//(*output)[1] = nullptr;
		//(*output)[2] = nullptr;
		return;
	}

	// % compute roots -------------------------------

	T t_vector[T_VECTOR_LEN]; common::colon(-1.0, 0.001, 1.0, t_vector);

	pose_poly<T> p;
	p.pose_from_point_tangents_2(gama1, tgt1, gama2, tgt2, Gama1, Tgt1, Gama2, Tgt2);

	static T root_ids[ROOT_IDS_LEN] = {0};
	p.find_bounded_root_intervals(t_vector, &root_ids);

	// % compute rhos, r, t --------------------------
	static T rhos[8][ROOT_IDS_LEN];
	p.rhos_from_root_ids(t_vector, root_ids, &rhos);
	T (&rhos1)[ROOT_IDS_LEN] = rhos[0];
	T (&rhos2)[ROOT_IDS_LEN] = rhos[3];
	T (&ts)[ROOT_IDS_LEN]    = rhos[6];
	T &ts_len                = rhos[7][0];


	// TODO: Check if this upper bound is correct (`for` loop and `sigmas1_end`)
	int sigma_len = ts_len * ts_len;

	T sigmas[4][sigma_len][sigma_len];

	// HACK: Turn multidimensional array into flat array
	T &sigmas_flat = sigmas[0][0][0];

	p.get_sigmas(ts_len, ts, &sigmas_flat);
	//T (&sigmas1)[sigma_len][sigma_len] = sigmas[0];
	//T (&sigmas2)[sigma_len][sigma_len] = sigmas[1];
	//T (&sigmas1_end)[sigma_len]        = sigmas[2][0];
	//T (&sigmas2_end)[sigma_len]        = sigmas[3][0];
	T *sigmas1     = &sigmas[0][0][0];
	T *sigmas2     = &sigmas[1][0][0];
	T *sigmas1_end = &sigmas[2][0][0];
	T *sigmas2_end = &sigmas[3][0][0];

	//T (&RT)[RT_LEN][4][3] = (*output)[0];
	T &RT = (*output)[0][0][0][0];

	p.get_r_t_from_rhos(
		ts_len,
		sigmas1, sigmas1_end,
		sigmas2, sigmas2_end,
		rhos1, rhos2,
		gama1, tgt1,
		gama2, tgt2,
		Gama1, Tgt1,
		Gama2, Tgt2,
		&RT
	);
}

}

