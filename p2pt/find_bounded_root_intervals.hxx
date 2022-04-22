#include "poly.h"  // TODO: Check if this include is really needed
#include "sample_pose_poly.hxx"

namespace P2Pt {

template<typename T>
void
pose_poly<T>::
find_bounded_root_intervals(const T t_vector[T_VECTOR_LEN], T root_ids[ROOT_IDS_LEN])
{
	// fvalue, A, B, C, E, F, G, H, J, K, L
	static T output[11][T_VECTOR_LEN] = { {0} };
	static T *sampled_poly = output[0];

	sample_pose_poly(t_vector, output);

	T curr_val = sampled_poly[0];
	for (int i = 0; i < ROOT_IDS_LEN; i++) {
		T nxt_val = sampled_poly[i + 1];
		root_ids[i] = (curr_val * nxt_val) < 0;
		curr_val = nxt_val;
	}
}

}