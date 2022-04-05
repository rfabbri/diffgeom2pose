namespace P2Pt {

template<typename T>
void 
pose_poly<T>::
rf_find_bounded_root_intervals(const T t_vector[2001], T root_ids[2001])
{
	static T A[2001];
	static T B[2001];
	static T C[2001];
	// TODO: check why `D` is missing
	static T E[2001];
	static T F[2001];
	static T G[2001];
	static T H[2001];
	// TODO: check why `I` is missing
	static T J[2001];
	static T K[2001];
	static T L[2001];
	static T fvalue[2001] = {0};

	static T *sampled_poly = fvalue;

	// TODO: use this separation of arguments in more functions instead of an aggregate array
	rf_sample_pose_poly(t_vector, A, B, C, E, F, G, H, J, K, L, fvalue);

	T curr_val = sampled_poly[0];
	for (int i = 0; i < 2001 /* sizeof(root_ids) does not work; returns pointer size */; i++) {
		T nxt_val = sampled_poly[i + 1];
		root_ids[i] = (curr_val * nxt_val) < 0;
		curr_val = nxt_val;
	}
}

}
