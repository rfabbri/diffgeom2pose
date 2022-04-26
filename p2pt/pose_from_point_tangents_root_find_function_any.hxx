#ifndef pose_from_point_tangents_root_find_function_any_hxx_
#define pose_from_point_tangents_root_find_function_any_hxx_

#include "common.hxx"
#include "poly.h"

#include "pose_from_point_tangents_2.hxx"
#include "find_bounded_root_intervals.hxx"
#include "rhos_from_root_ids.hxx"
#include "get_sigmas.hxx"
#include "get_r_t_from_rhos.hxx"

namespace P2Pt {

template<typename T>
void pose_from_point_tangents_root_find_function_any(
	const T (&gama1)[3], const T (&tgt1)[3],
	const T (&gama2)[3], const T (&tgt2)[3],
	const T (&Gama1)[3], const T (&Tgt1)[3],
	const T (&Gama2)[3], const T (&Tgt2)[3],
	T (*output_RT)[RT_MAX_LEN][4][3],
	int *output_RT_len,
	T *output_degen
)
{
	// % This is the main routine to find roots. Can be used with any input.


	// % test for geometric degeneracy -------------------------------
	static T DGama[3];
	common::vec1vec2_3el_sub(Gama1, Gama2, DGama);
	common::vec_3el_div_by_scalar(common::norm(DGama, 3), DGama, DGama);

	// Matrix for degeneracy calculation
	const T degen_matrix[3][3] = {
		DGama[0], Tgt1[0], Tgt2[0],
		DGama[1], Tgt1[1], Tgt2[1],
		DGama[2], Tgt1[2], Tgt2[2]
	};
	T &degen = *output_degen;
	degen = common::det3x3(degen_matrix);

	if (std::abs(degen) < 1.0e-3) {
		std::cout << "data point not reliable" << std::endl;
		output_RT     = nullptr;
		output_RT_len = nullptr;
		output_degen  = nullptr;
		return;
	}

	// % compute roots -------------------------------
	pose_poly<T> p;
	p.pose_from_point_tangents_2(
		gama1, tgt1,
		gama2, tgt2,
		Gama1, Tgt1,
		Gama2, Tgt2
	);

	static T t_vector[T_VECTOR_LEN]; common::colon(-1.0, 0.001, 1.0, t_vector);
	static T root_ids[ROOT_IDS_LEN] = {0};

	p.find_bounded_root_intervals(t_vector, &root_ids);


	// % compute rhos, r, t --------------------------
	static T rhos[3][ROOT_IDS_LEN];
	static int ts_len;

	p.rhos_from_root_ids(t_vector, root_ids, &rhos, &ts_len);

	T (&ts)[ROOT_IDS_LEN]    = rhos[0];
	T (&rhos1)[ROOT_IDS_LEN] = rhos[1];
	T (&rhos2)[ROOT_IDS_LEN] = rhos[2];


	static T sigmas[2][TS_MAX_LEN][TS_MAX_LEN];
	static int sigmas_len[2][TS_MAX_LEN];

	p.get_sigmas(ts_len, ts, &sigmas, &sigmas_len);

	T (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN] = sigmas[0];
	T (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN] = sigmas[1];
	int (&sigmas1_len)[TS_MAX_LEN]       = sigmas_len[0];
	int (&sigmas2_len)[TS_MAX_LEN]       = sigmas_len[1];


	T (&RT)[RT_MAX_LEN][4][3] = *output_RT;
	int &RT_len               = *output_RT_len;

	p.get_r_t_from_rhos(
		ts_len,
		sigmas1, sigmas1_len,
		sigmas2, sigmas2_len,
		rhos1, rhos2,
		gama1, tgt1,
		gama2, tgt2,
		Gama1, Tgt1,
		Gama2, Tgt2,
		&RT, &RT_len
	);
}

}

#endif // !pose_from_point_tangents_root_find_function_any_hxx_

