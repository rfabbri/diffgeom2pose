#ifndef find_bounded_root_intervals_hxx_
#define find_bounded_root_invervals_hxx_

#include "sample_pose_poly.hxx"

namespace P2Pt {

template<typename T>
void
pose_poly<T>::
find_bounded_root_intervals(const T (&t_vector)[T_VECTOR_LEN], T (*root_ids_output)[ROOT_IDS_LEN])
{
	// fvalue, A, B, C, E, F, G, H, J, K, L
	static T sampled_poly[T_VECTOR_LEN] = {0};

	sample_pose_poly(t_vector, &sampled_poly);

	static T curr_val, nxt_val;

	curr_val = sampled_poly[0];
	for (int i = 0; i < ROOT_IDS_LEN; i++) {
		nxt_val = sampled_poly[i+1];
		(*root_ids_output)[i] = (curr_val * nxt_val) < 0;
		curr_val = nxt_val;
	}
}

}

#endif // !find_bounded_root_intervals_hxx_

