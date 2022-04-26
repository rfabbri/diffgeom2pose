#include "poly.h"
#include "pose_from_point_tangents_2_fn_t_for_root.hxx"
#include <boost/math/tools/roots.hpp>

namespace P2Pt {

template <typename T>
void
pose_poly<T>::
rhos_from_root_ids(
	const T (&t_vector)[T_VECTOR_LEN], const T (&root_ids)[ROOT_IDS_LEN],
	T (*output)[7][ROOT_IDS_LEN], /* = {rhos1, rhos1_minus, rhos1_plus, rhos2, rhos2_minus, rhos2_plus, ts} */
	int *output_ts_len
)
{
	T (&rhos1)[ROOT_IDS_LEN]       = (*output)[0];
	T (&rhos1_minus)[ROOT_IDS_LEN] = (*output)[1];
	T (&rhos1_plus)[ROOT_IDS_LEN]  = (*output)[2];
	T (&rhos2)[ROOT_IDS_LEN]       = (*output)[3];
	T (&rhos2_minus)[ROOT_IDS_LEN] = (*output)[4];
	T (&rhos2_plus)[ROOT_IDS_LEN]  = (*output)[5];
	T (&ts)[ROOT_IDS_LEN]          = (*output)[6];

	static int ts_end;
	ts_end = 0;

	for (int i = 0; i < ROOT_IDS_LEN; i++) {
		if (root_ids[i] == 1) {
			static T t_ref;

			unsigned long max_iter = 7; // TODO: 7 seems to be the ideal # of iters
			std::pair<T,T> t_ref_pair = boost::math::tools::toms748_solve(
				*this,
				t_vector[i],
				t_vector[i+1],
				boost::math::tools::eps_tolerance<T>(),
				max_iter
			);

			t_ref = (t_ref_pair.first + t_ref_pair.second)/2.0;

			// TODO: Possibly optimize the size of `ts[]`.
			// What is the max number of 1s than can appear in `root_ids[]`?
			ts[ts_end++] = t_ref;
		}
	}
	*output_ts_len = ts_end;

	T t_stddev = t_vector[1] - t_vector[0];

	#define ALPHA_TS_COS(x) (2 * alpha * (x) * cos(theta))
	#define ALPHA_TS_SIN(x) (-2 * alpha * (x) * sin(theta))
	#define BETA_TS_SIN(x) (beta * (1 - (x) * (x)) * sin(theta))
	#define BETA_TS_COS(x) (beta * (1 - (x) * (x)) * cos(theta))
	#define TS_DEN(x) (1 + (x) * (x))

	//% Each root is now ts(i), plus minus t_stddev.
	//% Now get rho1(t):

	for (int i = 0; i < ts_end; i++) {
		T ts_new = ts[i];
		rhos1[i] = ALPHA_TS_COS(ts_new) + BETA_TS_SIN(ts_new);
		rhos1[i] /= TS_DEN(ts_new);
		rhos2[i] = ALPHA_TS_SIN(ts_new) + BETA_TS_COS(ts_new);
		rhos2[i] /= TS_DEN(ts_new);
	}
	for (int i = 0; i < ts_end; i++) {
		T ts_new = ts[i] - t_stddev;
		rhos1_minus[i] = ALPHA_TS_COS(ts_new) + BETA_TS_SIN(ts_new);
		rhos1_minus[i] /= TS_DEN(ts_new);
		rhos2_minus[i] = ALPHA_TS_SIN(ts_new) + BETA_TS_COS(ts_new);
		rhos2_minus[i] /= TS_DEN(ts_new);
	}
	for (int i = 0; i < ts_end; i++) {
		T ts_new = ts[i] + 2 * t_stddev;
		rhos1_plus[i] = ALPHA_TS_COS(ts_new) + BETA_TS_SIN(ts_new);
		rhos1_plus[i] /= TS_DEN(ts_new);
		rhos2_plus[i] = ALPHA_TS_SIN(ts_new) + BETA_TS_COS(ts_new);
		rhos2_plus[i] /= TS_DEN(ts_new);
	}

}

}


