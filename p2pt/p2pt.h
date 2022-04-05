#ifndef p2pt_h_
#define p2pt_h_

namespace P2Pt {
  
template <typename T=double>
class p2pt { // fully static, not to be instantiated - just used for templating
	public:
	static void hello();
	static void rf_find_bounded_root_intervals(const T t_vector[2001], T root_ids[2001]);
	static void rf_sample_pose_poly(
		const T t[2001], 
		T A[2001], T B[2001], T C[2001], T E[2001], T F[2001],
		T G[2001], T H[2001], T J[2001], T K[2001], T L[2001], T fvalue[2001]
	);
};
  
} // namespace minus

#endif  // p2pt_h_
