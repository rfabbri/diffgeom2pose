#include "common.hxx"
#include "poly.h"

#include "rf_pose_from_point_tangents_2.hxx"
#include "rf_find_bounded_root_intervals.hxx"

template<typename T>
void rf_pose_from_point_tangents_root_find_function_any(
	const T gama1[3], const T tgt1[3],
	const T gama2[3], const T tgt2[3],
	const T Gama1[3], const T Tgt1[3],
	const T Gama2[3], const T Tgt2[3],
	T output[3]
)
{
	// % This is the main routine to find roots. Can be used with any input.

	T *Rots    = &result[0];
	T *Transls = &result[1];
	T *degen   = &result[2];


	// % test for geometric degeneracy -------------------------------

	static T DGama[3];
	common::vec1vec2_sub(Gama1, Gama2, DGama);
	DGama = Dgama / common::norm(DGama, DGama_len);

	// Matrix for degeneracy calculation
	static T degen_matrix[3] = { DGama, Tgt1, Tgt2 };
	*degen = common::det3x3(degen_matrix);

	if (abs(*degen) < 1.0e-3) {
		std::cout << "data point not reliable" << std::endl;
		output[0] = nullptr;
		output[1] = nullptr;
		output[2] = nullptr;
		return;
	}

	// % compute roots -------------------------------

	T t_vector[2001];
	common::colon(-1, 0.001, 1, t_vector);

	// TODO: Improve later
	// TODO: Think if it needs to be static or not
	pose_poly poly;
	poly.rf_pose_from_point_tangents_2(gama1, tgt1, gama2, tgt2, Gama1, Tgt1, Gama2, Tgt2);

	static T root_ids[2001];
	poly.rf_find_bounded_root_intervals(t_vector, root_ids);

	// % compute rhos, r, t --------------------------
	//rf_rhos_from_root_ids(t_vector, root_ids); // TODO: implement `rf_rhos_from_root_ids()`

	//#include "rf_get_sigmas.hxx"
	//#include "rf_get_r_t_from_rhos.hxx"
}
