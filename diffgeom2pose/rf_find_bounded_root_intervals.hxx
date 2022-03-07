#include <assert.h>
#include <stdlib.h>
#include <string.h>

template<typename T>
T** rf_find_bounded_root_intervals(T* t_vector)
{
	// Output array MUST have size 2
	assert(output_size == 2);

	T* root_ids = (T*)malloc(sizeof(t_vector));      // TODO: check if calloc is needed in place of malloc
	T* sampled_poly = rf_sample_pose_poly(t_vector);

	T curr_val = sampled_poly[0];
	for (int i = 0; i < sizeof(root_ids); i++) {
		T nxt_val = sampled_poly[i + 1];
		root_ids[i] = (curr_val * nxt_val) < 0;
		curr_val = nxt_val;
	}

	// [root_ids, sampled_poly]
	T** output = (T**)malloc(sizeof(2 * T*));
	output[0] = root_ids;
	output[1] = sampled_poly;
	return output;
}