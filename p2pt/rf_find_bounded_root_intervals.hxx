#include "poly.h"  // TODO: Check if this include is really needed
#include "rf_sample_pose_poly.hxx"

namespace P2Pt {

template<typename T>
void
pose_poly<T>::
rf_find_bounded_root_intervals(const T t_vector[t_vector_len], T root_ids[root_ids_len])
{
	// fvalue, A, B, C, E, F, G, H, J, K, L
	static T output[11][t_vector_len] = { {0} };
	static T *sampled_poly = output[0];

	rf_sample_pose_poly(t_vector, output);

	T curr_val = sampled_poly[0];
	for (int i = 0; i < root_ids_len; i++) {
		T nxt_val = sampled_poly[i + 1];
		root_ids[i] = (curr_val * nxt_val) < 0;
		curr_val = nxt_val;
	}
}

}
