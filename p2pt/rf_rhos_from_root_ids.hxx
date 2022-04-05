template <typename T>
void
rf_rhos_from_root_ids(
	T t_vector[2001], T root_ids[2001], 
	T output[7] /* = {rhos1, rhos1_minus, rhos1_plus, rhos2, rhos2_minus, rhos2_plus, ts} */
)
{
	extern double alpha, beta, theta;

	static T ts[2001];
	// TODO: Check for out-of-bounds access on `i+1`
	for (int i = 0; i < 2001; i++) {
		if (root_ids[i] == 1) {
			rf_pose_from_point_tangents_2_fn_t_for_root(t_vector[i]);
			rf_pose_from_point_tangents_2_fn_t_for_root(t_vector[i+1]); 

			// TODO: check implementation of `fzero()`
			T t_ref = fzero(@rf_pose_from_point_tangents_2_fn_t_for_root, [t_vector(i) t_vector(i+1)]);;
			rf_pose_from_point_tangents_2_fn_t_for_root(t_ref);
			// TODO: Check total size/appendng of elements for `ts[]` as in MATLAB
			ts[i] = t_ref;
		}
	}
}