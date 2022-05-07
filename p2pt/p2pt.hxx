#ifndef p2pt_hxx_
#define p2pt_hxx_

#include <iostream>
#include "p2pt.h"
#include "pose_poly.hxx"
// thesse are pose_poly.hxx actually
#include "get_r_t_from_rhos.hxx"
#include "get_sigmas.hxx"
#include "pose_from_point_tangents_2.hxx"


namespace P2Pt {
  
// This is the main routine
template <typename T>
void p2pt<T>::
pose_from_point_tangents(
	const T (&gama1)[3], const T (&tgt1)[3],
	const T (&gama2)[3], const T (&tgt2)[3],
	const T (&Gama1)[3], const T (&Tgt1)[3],
	const T (&Gama2)[3], const T (&Tgt2)[3],
	T (*output_RT)[RT_MAX_LEN][4][3], int *output_RT_len,
	T *output_degen
)
{
	T DGama[3] = { Gama1[0] - Gama2[0], Gama1[1] - Gama2[1], Gama1[2] - Gama2[2] };
  { // % test for geometric degeneracy -------------------------------
  const T norm = sqrt(DGama[0]*DGama[0] + DGama[1]*DGama[1] + DGama[2]*DGama[2]);
  
	// Matrix for degeneracy calculation
	const T d[3][3] = {
		DGama[0]/norm, Tgt1[0], Tgt2[0],
		DGama[1]/norm, Tgt1[1], Tgt2[1],
		DGama[2]/norm, Tgt1[2], Tgt2[2]
	};
	T &degen = *output_degen;
	degen = (d[0][0]*d[1][1]*d[2][2]+d[0][1]*d[1][2]*d[2][0]+d[0][2]*d[1][0]*d[2][1]) // det(d)
		     -(d[2][0]*d[1][1]*d[0][2]+d[2][1]*d[1][2]*d[0][0]+d[2][2]*d[1][0]*d[0][1]);

	if (std::abs(degen) < 1.0e-3) {
    //TODO(openmvg) GUARD this print
		std::cerr << "data point not reliable, please skip this solve within RANSAC" << std::endl;
		output_RT     = nullptr;
		output_RT_len = nullptr;
		output_degen  = nullptr;
		return;
	}
  }

	// % compute roots -------------------------------
	pose_poly<T> p;
	p.pose_from_point_tangents_2( gama1, tgt1, gama2, tgt2, Gama1, Tgt1, Gama2, Tgt2);

	T root_ids[ROOT_IDS_LEN] __attribute__((aligned (16)));
	p.find_bounded_root_intervals(&root_ids);

	// % compute rhos, r, t --------------------------
	T rhos[3][ROOT_IDS_LEN];
	int ts_len;

	p.rhos_from_root_ids(root_ids, &rhos, &ts_len);

	const T (&ts)[ROOT_IDS_LEN]    = rhos[0];
	const T (&rhos1)[ROOT_IDS_LEN] = rhos[1];
	const T (&rhos2)[ROOT_IDS_LEN] = rhos[2];

	T sigmas[2][TS_MAX_LEN][TS_MAX_LEN];
	int sigmas_len[2][TS_MAX_LEN];

	p.get_sigmas(ts_len, ts, &sigmas, &sigmas_len);

	T (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN] = sigmas[0];
	T (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN] = sigmas[1];
	const int (&sigmas1_len)[TS_MAX_LEN] = sigmas_len[0];
	const int (&sigmas2_len)[TS_MAX_LEN] = sigmas_len[1];

	T (&RT)[RT_MAX_LEN][4][3] = *output_RT;
	int &RT_len               = *output_RT_len;

	p.get_r_t_from_rhos(
		ts_len,
		sigmas1, sigmas1_len, sigmas2, sigmas2_len,
		rhos1, rhos2,
		gama1, tgt1, gama2, tgt2,
		Gama1, Tgt1, Gama2, Tgt2,
		&RT, &RT_len);
}

} // namespace p2pt

#endif // !p2pt_hxx_

