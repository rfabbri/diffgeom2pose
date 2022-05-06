#ifndef pose_poly_hxx_
#define pose_poly_hxx_

#include "pose_poly.h"

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
			unsigned long max_iter = 7; // INFO: 7 seems to be the ideal # of iters

      std::cout << "started iterative method --------------------------------------\n";
      std::cerr << "current eval count " << std::endl;
      // secant max_iter = 3: 16 evals total. slower than boost,
      // I geuss since boost uses quadratic interpolation after 1st iteration
      /*
      double t0 = t_vector[i], t1 = t_vector[i+1], t2;
      double f0 = fn_t(t_vector[i]), f1 = fn_t(t_vector[i+1]);
      std::cerr << "XXXXXXXX  f0, f1 " << f0 << " " << f1 << std::endl;
      for (unsigned k=0; k < max_iter; ++k) {
        t2 = t1 - f1 * (t1 - t0) / (f1 - f0);
        t0 = t1; t1 = t2;
        f0 = f1; if (k + 1 < max_iter) f1 = fn_t(t2);
        std::cerr << "XXXXXXXX  f(x2)" << f1 << std::endl;
      }

      ts[ts_end++] = t2;
      */
      
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
}

} // ! namespace P2Pt

#endif // !pose_poly_hxx_
