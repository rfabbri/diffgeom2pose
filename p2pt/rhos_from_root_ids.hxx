#include "poly.h"
#include "pose_from_point_tangents_2_fn_t_for_root.hxx"

namespace P2Pt {

template <typename T>
void
pose_poly<T>::
rhos_from_root_ids(
	const T t_vector[t_vector_len], const T root_ids[t_vector_len],
	T output[7][t_vector_len] /* = {rhos1, rhos1_minus, rhos1_plus, rhos2, rhos2_minus, rhos2_plus, ts} */
)
{
	T* rhos1       = output[0];
	T* rhos1_minus = output[1];
	T* rhos1_plus  = output[2];
	T* rhos2       = output[3];
	T* rhos2_minus = output[4];
	T* rhos2_plus  = output[5];
	T* ts          = output[6];

	int ts_end = 0;
	for (int i = 0; i < root_ids_len; i++) {
		if (root_ids[i] == 1) {
			// TODO: implement fzero using (t_vector(i) + t_vector(i+1))/2; maybe use Newton's method later
			static T t_ref_arr[11];
			static T t_ref;

			pose_from_point_tangents_2_fn_t_for_root((t_vector[i] + t_vector[i+1])/2, t_ref_arr);
			t_ref = t_ref_arr[0];

			//pose_from_point_tangents_2_fn_t_for_root(t_ref);

			// TODO: Possibly optimize the size of `ts[]`.
			// What is the max number of 1s than can appear in `root_ids[]`?
			ts[ts_end++] = t_ref;
		}
	}

	T t_stddev = t_vector[1] - t_vector[0];

	#define ALPHA_TS_COS(x) (2 * alpha * (x) * cos(theta))
	#define ALPHA_TS_SIN(x) (-2 * alpha * (x) * sin(theta))
	#define BETA_TS_SIN(x) (beta * (1 - (x) * (x)) * sin(theta))
	#define BETA_TS_COS(x) (beta * (1 - (x) * (x)) * cos(theta))
	#define TS_DEN(x) (1 + (x) * (x))

	//% Each root is now ts(i), plus minus t_stddev.
	//% Now get rho1(t):

	for (int i = 0; i < t_vector_len; i++) {
		T ts_new = ts[i];
		rhos1[i] = ALPHA_TS_COS(ts_new) + BETA_TS_SIN(ts_new);
		rhos1[i] /= TS_DEN(ts_new);
		rhos2[i] = ALPHA_TS_SIN(ts_new) + BETA_TS_COS(ts_new);
		rhos2[i] /= TS_DEN(ts_new);
	}
	for (int i = 0; i < t_vector_len; i++) {
		T ts_new = ts[i] - t_stddev;
		rhos1_minus[i] = ALPHA_TS_COS(ts_new) + BETA_TS_SIN(ts_new);
		rhos1_minus[i] /= TS_DEN(ts_new);
		rhos2_minus[i] = ALPHA_TS_SIN(ts_new) + BETA_TS_COS(ts_new);
		rhos2_minus[i] /= TS_DEN(ts_new);
	}
	for (int i = 0; i < t_vector_len; i++) {
		T ts_new = ts[i] + 2 * t_stddev;
		rhos1_plus[i] = ALPHA_TS_COS(ts_new) + BETA_TS_SIN(ts_new);
		rhos1_plus[i] /= TS_DEN(ts_new);
		rhos2_plus[i] = ALPHA_TS_SIN(ts_new) + BETA_TS_COS(ts_new);
		rhos2_plus[i] /= TS_DEN(ts_new);
	}

}

}


