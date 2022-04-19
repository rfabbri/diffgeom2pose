#include "rf_pose_from_point_tangents_2_fn_t_for_root.hxx"

template <typename T>
void
rf_rhos_from_root_ids(
	T t_vector[t_vector_len], T root_ids[t_vector_len],
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
	// TODO: Create element counter for ts

	// TODO: Make member of poly? Needs access to `alpha`, `beta`, and `theta`
	extern T alpha, beta, theta;

	for (int i = 0; i < root_ids_len; i++) {
		if (root_ids[i] == 1) {
			// TODO: implement fzero using (t_vector(i) + t_vector(i+1))/2; maybe use Newton's method later
			static T t_ref_arr[11];
			static T t_ref;

			// First round
			rf_pose_from_point_tangents_2_fn_t_for_root((t_vector[i] + t_vector[i])/2, t_ref_arr); 
			t_ref = t_ref_arr[0];

			// Second round
			rf_pose_from_point_tangents_2_fn_t_for_root(t_ref/2, t_ref_arr); 
			t_ref = t_ref_arr[0];



			//rf_pose_from_point_tangents_2_fn_t_for_root(t_ref);

			// TODO: Check total size/appending of elements for `ts[]` as in MATLAB
			// TODO: Possibly optimize the size of `ts[]`.
			// What is the max number of 1s than can appear in `root_ids[]`?
			ts[i] = t_ref;
		}
	}

	T t_stddev = t_vector[1] - t_vector[0];

	// TODO: Check if inline functions are a better fit
	template<typename T>
	constexpr auto ALPHA_TS_COS(T x) { return (2 * alpha * (x) * cos(theta)); }
	template<typename T>
	constexpr auto ALPHA_TS_SIN(T x) { return (-2 * alpha * (x) * sin(theta)); }
	template<typename T>
	constexpr auto BETA_TS_SIN(T x) { return (beta * (1 - (x) * (x)) * sin(theta)); }
	template<typename T>
	constexpr auto BETA_TS_COS(T x) { return (beta * (1 - (x) * (x)) * cos(theta)); }
	template<typename T>
	constexpr auto TS_DEN(T x) { return (1 + (x) * (x)); }

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

