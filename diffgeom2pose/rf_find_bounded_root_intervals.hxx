#include <assert.h>

template <typename T>
void rf_find_bounded_root_intervals(T* t_vector, T** output, int output_size)
{
	// Output array MUST have size 2
	assert(output_size == 2);

	T root_ids[sizeof(t_vector)] = { 0 };
	T sampled_poly[11] = rf_sample_pose_poly(t_vector); // TODO: implement `rf_sample_pose_poly()`

	T curr_val = sampled_poly[0];
	for (int i = 0; i < sizeof(root_ids); i++) {
		T nxt_val = sampled_poly[i + 1];
		root_ids[i] = (curr_val * nxt_val) < 0;
		curr_val = nxt_val;
	}

	// TODO: FIX ME!!! Local variable gets out-of-scope on return
	output[0] = &root_ids;
	output[1] = &sampled_poly;
}