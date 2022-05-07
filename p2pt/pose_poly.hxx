#ifndef pose_poly_hxx_
#define pose_poly_hxx_

#include "pose_poly.h"

namespace P2Pt {
  
template <typename T>
void
pose_poly<T>::
rhos_from_root_ids(
	const T (&root_ids)[ROOT_IDS_LEN],
	T (*output)[3][ROOT_IDS_LEN],
	int *output_ts_len
)
{
	T (&ts)[ROOT_IDS_LEN] = (*output)[0];

	int &ts_end = *output_ts_len; ts_end = 0;
	for (unsigned i = 0; i < ROOT_IDS_LEN; i++) {
		if (!root_ids[i]) continue;
    T t0 = t_vec(i), t1 = t_vec(i+1), &t2 = ts[ts_end++];
    T f0 = fn_t(t_vec(i)), f1 = fn_t(t_vec(i+1));
    for (unsigned k = 0; k < 3; ++k) {
      t2 = t1 - f1*(t1-t0)/(f1-f0); t0 = t1; t1 = t2;
      f0 = f1; if (k + 1 < 3) f1 = fn_t(t2);
    }
	}

	//% Each root is now ts(i), plus minus t_stddev.
	//% Now get rho1(t):

	T (&rhos1)[ROOT_IDS_LEN] = (*output)[1]; T (&rhos2)[ROOT_IDS_LEN] = (*output)[2];

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
}

} // ! namespace P2Pt

#endif // !pose_poly_hxx_

