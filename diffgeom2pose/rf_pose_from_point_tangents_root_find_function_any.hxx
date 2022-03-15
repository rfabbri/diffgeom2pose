#include "common.hxx"

#define GAMA_SIZE 3 // TODO: change this later!

template<typename T>
T **rf_pose_from_point_tangents_root_find_function_any(
	T *gama1, T *tgt1,
	T *gama2, T *tgt2,
	T *Gama1, T *Tgt1,
	T *Gama2, T *Tgt2
)
{
	// % This is the main routine to find roots. Can be used with any input.

	// TODO: get number of "solutions" `n`??? to initialize fixed size array 
	static T Rots[]    = {};
	static T Transls[] = {};
	static T degen[1];

	static T *output[3];

	// % test for geometric degeneracy -------------------------------

	// TODO: get length of DGama vector
	T DGama = Gama1 - Gama2;
	T* DGama = common::vec_subtract(Gama1, GAMA_SIZE, Gama2, GAMA_SIZE);
	DGama = Dgama / common::norm(DGama, DGama_len);

	// Matrix for degeneracy calculation
	T degen_matrix[] = { DGama, Tgt1, Tgt2 };
	degen[0] = common::det3x3(degen_matrix);

	if (abs(degen) < 1.0e-3) {
		std::cout << "data point not reliable" << std::endl;
		output[0] = nullptr;
		output[1] = nullptr;
		output[2] = nullptr;
		return output;
	}

	// % compute roots -------------------------------

	T *t_vector = common::colon(-1, 0.001, 1);
	#include rf_pose_from_point_tangents_2

	T **root_ids = rf_find_bounded_root_intervals(t_vector);

	// % compute rhos, r, t --------------------------
	rf_rhos_from_root_ids(t_vector, root_ids); // TODO: implement `rf_rhos_from_root_ids()`

	#include rf_get_sigmas.hxx
	#include rf_get_r_t_from_rhos.hxx
}
