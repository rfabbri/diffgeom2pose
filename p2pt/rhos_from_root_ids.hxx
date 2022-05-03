#ifndef rhos_from_root_ids_hxx_
#define rhos_from_root_ids_hxx_

#include "pose_from_point_tangents_2_fn_t.hxx"
#include <boost/math/tools/roots.hpp>

namespace P2Pt {

template <typename T>
void
pose_poly<T>::
rhos_from_root_ids(
	const T (&t_vector)[T_VECTOR_LEN], const T (&root_ids)[ROOT_IDS_LEN],
	T (*output)[3][ROOT_IDS_LEN],
	int *output_ts_len
)
{
	T (&ts)[ROOT_IDS_LEN] = (*output)[0];

	int ts_end = 0;
	std::pair<T,T> t_ref_pair;

	for (int i = 0; i < ROOT_IDS_LEN; i++) {
		if (root_ids[i] == 1) {
			static unsigned long max_iter = 7; // INFO: 7 seems to be the ideal # of iters

			try {
				t_ref_pair = boost::math::tools::toms748_solve(
					*this,
					t_vector[i],
					t_vector[i+1],
					boost::math::tools::eps_tolerance<T>(),
					max_iter
				);
			} catch (const std::exception& err) {
        // TODO(OpenMVG) take care of this error without cerr
				std::cerr << err.what() << std::endl;
			}
      
			// TODO: Possibly optimize the size of `ts[]`.
			// What is the max number of 1s than can appear in `root_ids[]`?
      ts[ts_end++] = (t_ref_pair.first + t_ref_pair.second)*0.5;
		}
	}
	*output_ts_len = ts_end;

	//% Each root is now ts(i), plus minus t_stddev.
	//% Now get rho1(t):

	T (&rhos1)[ROOT_IDS_LEN] = (*output)[1];
	T (&rhos2)[ROOT_IDS_LEN] = (*output)[2];

  const T alpha_times_2 = 2.*alpha;
	for (int i = 0; i < ts_end; i++) {
		const T ts_new = ts[i],
    x2 = ts_new * ts_new,
    ts_den = 1. + x2,
    alpha_ts_new2 = alpha_times_2 * ts_new,
    beta_1_minus_x2 = beta * (1. - x2);
    
		rhos1[i] = (( alpha_ts_new2 * cth) + (beta_1_minus_x2 * sth)) / ts_den;
		rhos2[i] = ((-alpha_ts_new2 * sth) + (beta_1_minus_x2 * cth)) / ts_den;
	}

	// Not used?
	#if 0
	static T t_stddev;
	t_stddev = t_vector[1] - t_vector[0];

	static T rhos1_minus, rhos1_plus;
	static T rhos2_minus, rhos2_plus;

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
	#endif

}

}

#endif // !rhos_from_root_ids_hxx_

