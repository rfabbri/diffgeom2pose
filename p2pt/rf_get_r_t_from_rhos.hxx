// % to be called from rf_pose_from_point_tangents_root_find_function_any.m

template<typename T>
void
rf_get_r_from_rhos()
{
	//% to be called from rf_pose_from_point_tangents_root_find_function_any.m

	// TODO: Check what the maximum number of elements can be for a cell array
	// TODO: Check if this can be implemented as a vector for pushback

	static constexpr int max_array_len = t_vector_len;
	//% Lambdas:
	static T lambdas1[max_array_len][max_array_len];
	static T lambdas2[max_array_len][max_array_len];

	static int end1 = 0;
	static int end2 = 0;

	for (int i = 0; i < t_vector_len; i++) {
		// TODO: implement `isempty`. Nullptr?
		if (isempty(sigmas1[i])) {
			lambdas1[i] = {};
			lambdas2[i] = {};
		}
		lambdas1[i] =
		lambdas2[i] =
	}

}
