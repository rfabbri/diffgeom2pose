#include <common.hxx>

template<typename T>
T** rf_pose_from_point_tangents_root_find_function_any(
	// TODO: GamaX are vectors; check if other values are as well.
	T gama1, T tgt1,
	T gama2, T tgt2,
	T Gama1, T Tgt1,
	T Gama2, T Tgt2
)
{
	// % This is the main routine to find roots. Can be used with any input.

	// Output array MUST have size 3
	assert(output_size == 3);

	T* Rots    = nullptr;
	T* Transls = nullptr;
	T degen;

	T** output = (T**)malloc(sizeof(3 * T*));

	// % test for geometric degeneracy -------------------------------

	// TODO: get length of DGama vector
	T DGama = Gama1 - Gama2;
	DGama = Dgama / common::norm(DGama, DGama_len);

	// Array for degeneracy calculation
	T degen_array[] = { DGama, Tgt1, Tgt2 };
	degen = common::det(degen_array);

	if (abs(degen) < 1.0e-3) {
		puts("data point not reliable");
		output[0] = Rots;
		output[1] = Transls;
		output[2] = degen;
		return;
	}

	// % compute roots -------------------------------

	T* t_vector = common::colon(-1, 0.001, 1);
	#include rf_pose_from_point_tangents_2

	T** root_ids = rf_find_bounded_root_intervals(t_vector);

	// % compute rhos, r, t --------------------------
	rf_rhos_from_root_ids(t_vector, root_ids); // TODO: implement `rf_rhos_from_root_ids()`

	#include rf_get_sigmas.hxx
	#include rf_get_r_t_from_rhos.hxx
}
