template <typename T>
void rf_pose_from_point_tangents_root_find_function_any(
	T gama1, T tgt1,
	T gama2, T tgt2,
	T Gama1, T Tgt1,
	T Gama2, T Tgt2
	T** output,
	int output_size
)
{
	// % This is the main routine to find roots. Can be used with any input.

	// Output array MUST have size 3
	assert(output_size == 3);

	// INFO: `Rots[]` and `Transls[]` should be initialized to null
	// INFO: FIX ME!!!
	T Rots[]    = {};
	T Transls[] = {};
	T degen;

	// % test for geometric degeneracy -------------------------------

	T DGama = Gama1 - Gama2;
	DGama = Dgama / norm(DGama); // TODO: implement `norm()`

	// Array for degeneracy calculation
	T degen_array[] = { DGama, Tgt1, Tgt2 };

	degen = det(degen_array); // TODO: implement `det()`

	// TODO: FIX ME!!! Local variable gets out-of-scope on return
	output[0] = &Rots;
	output[1] = &Transls;
	//output[2] = 

	if (abs(degen) < 1.0e-3) { // TODO: implement `abs()`
		puts("data point not reliable");
	}

	// % compute roots -------------------------------

	T* t_vector = linspace(-1, 0.001, 1); // TODO: implment `linspace()`
	#include rf_pose_from_point_tangents_2

	root_ids = (T**)malloc(2 * sizeof(T*));
	rf_find_bounded_root_intervals(t_vector, root_ids);

	// % compute rhos, r, t --------------------------
	rf_rhos_from_root_ids(t_vector, root_ids); // TODO: implement `rf_rhos_from_root_ids()`

	#include rf_get_sigmas.hxx
	#include rf_get_r_t_from_rhos.hxx
}
