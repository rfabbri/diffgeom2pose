#ifndef find_bounded_root_intervals_hxx_
#define find_bounded_root_intervals_hxx_

//#include "sample_pose_poly.hxx"
#include "pose_from_point_tangents_2_fn_t.hxx"

namespace P2Pt {

template<typename T>
void
pose_poly<T>::
find_bounded_root_intervals(const T (&t_vector)[T_VECTOR_LEN], T (*root_ids_output)[ROOT_IDS_LEN])
{
	//static T sampled_poly[T_VECTOR_LEN] __attribute__((aligned (16))) = {0};

	//sample_pose_poly(t_vector, &sampled_poly);

	static T curr_val, next_val;

	curr_val = pose_from_point_tangents_2_fn_t(t_vector[0]);
	for (int i = 0; i < ROOT_IDS_LEN; i++) {
		next_val = pose_from_point_tangents_2_fn_t(t_vector[i+1]);
		(*root_ids_output)[i] = (curr_val * next_val) < 0;
		curr_val = next_val;
	}

	//curr_val = sampled_poly[0];
	//for (int i = 0; i < ROOT_IDS_LEN; i++) {
	//	next_val = sampled_poly[i+1];
	//	(*root_ids_output)[i] = (curr_val * next_val) < 0;
	//	curr_val = next_val;
	//}
}

}

#endif // !find_bounded_root_intervals_hxx_

