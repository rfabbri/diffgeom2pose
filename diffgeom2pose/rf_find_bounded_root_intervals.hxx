
template<typename T>
T **rf_find_bounded_root_intervals(T *t_vector, int t_vector_size)
{

	static T root_ids[t_vector_size];
	T *sampled_poly = rf_sample_pose_poly(t_vector);

	T curr_val = sampled_poly[0];
	for (int i = 0; i < sizeof(root_ids); i++) {
		T nxt_val = sampled_poly[i + 1];
		root_ids[i] = (curr_val * nxt_val) < 0;
		curr_val = nxt_val;
	}

	// [root_ids, sampled_poly]
	static T *output[] = {root_ids, sampled_poly};
	return output;
}