#ifndef get_r_t_from_rhos_hxx_
#define get_r_t_from_rhos_hxx_

namespace P2Pt {

template<typename T>
void
pose_poly<T>::
get_r_t_from_rhos(
	const int ts_len,
	const T (&sigmas1)[TS_MAX_LEN][TS_MAX_LEN], const int (&sigmas1_len)[TS_MAX_LEN],
	const T (&sigmas2)[TS_MAX_LEN][TS_MAX_LEN], const int (&sigmas2_len)[TS_MAX_LEN],
	const T (&rhos1)[ROOT_IDS_LEN], const T (&rhos2)[ROOT_IDS_LEN],
	const T (&gama1)[3], const T (&tgt1)[3],
	const T (&gama2)[3], const T (&tgt2)[3],
	const T (&Gama1)[3], const T (&Tgt1)[3],
	const T (&Gama2)[3], const T (&Tgt2)[3],
	T (*output)[RT_MAX_LEN][4][3], int *output_len
)
{
	T lambdas1[TS_MAX_LEN][TS_MAX_LEN];
	T lambdas2[TS_MAX_LEN][TS_MAX_LEN];
	const T DGama[3] = {
    Gama1[0] - Gama2[0],
    Gama1[1] - Gama2[1],
    Gama1[2] - Gama2[2]
  };
  
	for (int i = 0; i < ts_len; i++) {
    assert(sigmas1_len[i] == sigmas2_len[i]);
    const T dgamas_rhos[3] = {
     rhos1[i]*gama1[0] - rhos2[i]*gama2[0],
     rhos1[i]*gama1[1] - rhos2[i]*gama2[1],
     rhos1[i]*gama1[2] - rhos2[i]*gama2[2]
    };
		for (int j = 0; j < sigmas1_len[i]; j++) {
			lambdas1[i][j] = 
        (DGama[0]*Tgt1[0]+DGama[1]*Tgt1[1] + DGama[2]*Tgt1[2]) / 
        (dgamas_rhos[0]*(rhos1[i]*tgt1[0] + sigmas1[i][j]*gama1[0]) + 
        dgamas_rhos[1]*(rhos1[i]*tgt1[1] + sigmas1[i][j]*gama1[1]) +
        dgamas_rhos[2]*(rhos1[i]*tgt1[2] + sigmas1[i][j]*gama1[2]));
      
			lambdas2[i][j] = 
        (DGama[0]*Tgt2[0]+DGama[1]*Tgt2[1] + DGama[2]*Tgt2[2]) /
        (dgamas_rhos[0]*(rhos2[i]*tgt2[0] + sigmas2[i][j]*gama2[0]) + 
        dgamas_rhos[1]*(rhos2[i]*tgt2[1] + sigmas2[i][j]*gama2[1]) +
        dgamas_rhos[2]*(rhos2[i]*tgt2[2] + sigmas2[i][j]*gama2[2]));
		}
	}

	//% Rotation:
	const T A[3][3] = {
		DGama[0], Tgt1[0], Tgt2[0],
		DGama[1], Tgt1[1], Tgt2[1],
		DGama[2], Tgt1[2], Tgt2[2],
	};

	T inv_A[3][3]; common::invm3x3(A, inv_A);

	// Matrix containing Rotations and Translations
	T (&RT)[RT_MAX_LEN][4][3] = *output;
	int &RT_len               = *output_len;

	RT_len = 0;
	for (int i = 0; i < ts_len; i++) {
		for (int j = 0; j < sigmas1_len[i]; j++, RT_len++) {
			T (&Rots)[4][3] = RT[RT_len];
			T (&Transls)[3] = RT[RT_len][3];

			#define B_row(r) \
				rhos1[i]*gama1[(r)] - rhos2[i]*gama2[(r)], \
				lambdas1[i][j]*(rhos1[i]*tgt1[(r)] + sigmas1[i][j]*gama1[(r)]), \
				lambdas2[i][j]*(rhos2[i]*tgt2[(r)] + sigmas2[i][j]*gama2[(r)])

			const T B[3][3] = {
				B_row(0),
				B_row(1),
				B_row(2)
			};

			common::multm3x3(B, inv_A, Rots);

			T buff2[3];

			// Transls{end+1} = rhos1(i)*gama1 - Rots{end}*Gama1;
			common::multm_3x3_3x1(Rots, Gama1, buff2);
      Transls[0] = rhos1[i]*gama1[0] - buff2[0];
      Transls[1] = rhos1[i]*gama1[1] - buff2[1];
      Transls[2] = rhos1[i]*gama1[2] - buff2[2];
		}
	}
}
}

#endif // !get_r_t_from_rhos_hxx_

