//
// \author Ricardo Fabbri based on original code by Anton Leykin
// \date February 2019
//
// Tests more comprehensive runs of p2pt using the public interface
//
#include <cstring>
#include <iostream>
#include <testlib/testlib_test.h>
#include <p2pt/p2pt.h>
#include <p2pt/poly.h>

#include <p2pt/pose_from_point_tangents_root_find_function_any.hxx>
#include <p2pt/rhos_from_root_ids.hxx>
#include <p2pt/get_sigmas.hxx>

#include <tests/test-p2pt-constants.hxx>

using namespace P2Pt;

static void
test_hello()
{
	p2pt<double>::hello();
	TEST("Rodou hello? ", true, true);
}

static void
test_pose_from_point_tangents_root_find_function_any()
{
	double output[2][RT_MAX_LEN + 1][4][3] = {0};

	pose_from_point_tangents_root_find_function_any(
		sample_gama1, sample_tgt1,
		sample_gama2, sample_tgt2,
		sample_Gama1, sample_Tgt1,
		sample_Gama2, sample_Tgt2,
		&output
	);

	double (&RT)[RT_MAX_LEN + 1][4][3] = output[0];
	double(&RT_len)    = RT[0][0][0];
	double &test_degen = output[1][0][0][0];
	char indexstr[128];

	// HACK: `i+1` used to skip first value of `RT` which is `RT_len`
	for (int i = 0; i < sample_RT_len; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				double (&test_Rots)[4][3] = RT[i+1];
				snprintf(indexstr, 128, "RT[%d]: Rots[%d][%d]", i+1, j, k);
				TEST_NEAR_REL(indexstr, test_Rots[j][k], sample_Rots[i][j][k], eps);
			}
		}
		for (int j = 0; j < 3; j++) {
			double (&test_Transls)[3] = RT[i+1][3];
			snprintf(indexstr, 128, "RT[%d]: Transls[%d]", i+1, j);
			TEST_NEAR_REL(indexstr, test_Transls[j], sample_Transls[i][j], eps);
		}
		std::cout << std::endl;
	}
	TEST_NEAR_REL("degen", test_degen, sample_degen, eps);
}

static void
test_pose_from_point_tangents_2()
{
	pose_poly<double> p;

	p.pose_from_point_tangents_2(
		sample_gama1, sample_tgt1,
		sample_gama2, sample_tgt2,
		sample_Gama1, sample_Tgt1,
		sample_Gama2, sample_Tgt2
	);

	TEST_NEAR_REL("A0", p.A0, sample_A0, eps);
	TEST_NEAR_REL("A1", p.A1, sample_A1, eps);
	TEST_NEAR_REL("A2", p.A2, sample_A2, eps);

	TEST_NEAR_REL("B0", p.B0, sample_B0, eps);
	TEST_NEAR_REL("B1", p.B1, sample_B1, eps);
	TEST_NEAR_REL("B2", p.B2, sample_B2, eps);
	TEST_NEAR_REL("B3", p.B3, sample_B3, eps);

	TEST_NEAR_REL("C0", p.C0, sample_C0, eps);
	TEST_NEAR_REL("C1", p.C1, sample_C1, eps);
	TEST_NEAR_REL("C2", p.C2, sample_C2, eps);
	TEST_NEAR_REL("C3", p.C3, sample_C3, eps);
	TEST_NEAR_REL("C4", p.C4, sample_C4, eps);

	TEST_NEAR_REL("E0", p.E0, sample_E0, eps);
	TEST_NEAR_REL("E1", p.E1, sample_E1, eps);
	TEST_NEAR_REL("E2", p.E2, sample_E2, eps);

	TEST_NEAR_REL("F0", p.F0, sample_F0, eps);
	TEST_NEAR_REL("F1", p.F1, sample_F1, eps);
	TEST_NEAR_REL("F2", p.F2, sample_F2, eps);
	TEST_NEAR_REL("F3", p.F3, sample_F3, eps);

	TEST_NEAR_REL("G0", p.G0, sample_G0, eps);
	TEST_NEAR_REL("G1", p.G1, sample_G1, eps);
	TEST_NEAR_REL("G2", p.G2, sample_G2, eps);
	TEST_NEAR_REL("G3", p.G3, sample_G3, eps);
	TEST_NEAR_REL("G4", p.G4, sample_G4, eps);

	TEST_NEAR_REL("H0", p.H0, sample_H0, eps);
	TEST_NEAR_REL("H1", p.H1, sample_H1, eps);
	TEST_NEAR_REL("H2", p.H2, sample_H2, eps);
	TEST_NEAR_REL("H3", p.H3, sample_H3, eps);
	TEST_NEAR_REL("H4", p.H4, sample_H4, eps);

	TEST_NEAR_REL("J0", p.J0, sample_J0, eps);
	TEST_NEAR_REL("J1", p.J1, sample_J1, eps);
	TEST_NEAR_REL("J2", p.J2, sample_J2, eps);
	TEST_NEAR_REL("J3", p.J3, sample_J3, eps);

	TEST_NEAR_REL("K0", p.K0, sample_K0, eps);
	TEST_NEAR_REL("K1", p.K1, sample_K1, eps);
	TEST_NEAR_REL("K2", p.K2, sample_K2, eps);
	TEST_NEAR_REL("K3", p.K3, sample_K3, eps);

	TEST_NEAR_REL("L0", p.L0, sample_L0, eps);
	TEST_NEAR_REL("L1", p.L1, sample_L1, eps);
	TEST_NEAR_REL("L2", p.L2, sample_L2, eps);
}

static void
test_find_bounded_root_intervals()
{
	static double test_root_ids[ROOT_IDS_LEN];

	pose_poly<double> p = {
		sample_A0, sample_A1, sample_A2,
		sample_B0, sample_B1, sample_B2, sample_B3,
		sample_C0, sample_C1, sample_C2, sample_C3, sample_C4,
		sample_E0, sample_E1, sample_E2,
		sample_F0, sample_F1, sample_F2, sample_F3,
		sample_G0, sample_G1, sample_G2, sample_G3, sample_G4,
		sample_H0, sample_H1, sample_H2, sample_H3, sample_H4,
		sample_J0, sample_J1, sample_J2, sample_J3,
		sample_K0, sample_K1, sample_K2, sample_K3,
		sample_L0, sample_L1, sample_L2
	};

	p.find_bounded_root_intervals(sample_t_vector, &test_root_ids);

	for (int i = 0; i < ROOT_IDS_LEN; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "root_ids[%d]", i);
		TEST_NEAR_REL(indexstr, test_root_ids[i], sample_root_ids[i], eps);
	}
}

static void
test_sample_pose_poly()
{
	pose_poly<double> p = {
		sample_A0, sample_A1, sample_A2,
		sample_B0, sample_B1, sample_B2, sample_B3,
		sample_C0, sample_C1, sample_C2, sample_C3, sample_C4,
		sample_E0, sample_E1, sample_E2,
		sample_F0, sample_F1, sample_F2, sample_F3,
		sample_G0, sample_G1, sample_G2, sample_G3, sample_G4,
		sample_H0, sample_H1, sample_H2, sample_H3, sample_H4,
		sample_J0, sample_J1, sample_J2, sample_J3,
		sample_K0, sample_K1, sample_K2, sample_K3,
		sample_L0, sample_L1, sample_L2
	};

	double output[11][T_VECTOR_LEN];

	p.sample_pose_poly(sample_t_vector, &output);

	double (&test_fvalue)[T_VECTOR_LEN] = output[0];

	for (int i = 0; i < T_VECTOR_LEN; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "fvalue[%d]", i);
		TEST_NEAR_REL(indexstr, test_fvalue[i], sample_fvalue_pose_poly[i], eps);
	}
}

static void
test_rhos_from_root_ids()
{
	pose_poly<double> p;
	p.alpha = sample_alpha;
	p.beta  = sample_beta;
	p.theta = sample_theta;

	double output[8][ROOT_IDS_LEN];
	char indexstr[128];

	p.rhos_from_root_ids(sample_t_vector, sample_root_ids, &output);

	double (&test_rhos1)[ROOT_IDS_LEN]       = output[0];
	double (&test_rhos1_minus)[ROOT_IDS_LEN] = output[1];
	double (&test_rhos1_plus)[ROOT_IDS_LEN]  = output[2];
	double (&test_rhos2)[ROOT_IDS_LEN]       = output[3];
	double (&test_rhos2_minus)[ROOT_IDS_LEN] = output[4];
	double (&test_rhos2_plus)[ROOT_IDS_LEN]  = output[5];
	double (&test_ts)[ROOT_IDS_LEN]          = output[6];
	double (&test_ts_len)                    = output[7][0];

	for (int i = 0; i < sample_ts_len; i++) {
		snprintf(indexstr, 128, "test_rhos1[%d]", i);
		TEST_NEAR_REL(indexstr, test_rhos1[i], sample_rhos1[i], eps);
	}
	for (int i = 0; i < sample_ts_len; i++) {
		snprintf(indexstr, 128, "test_rhos1_minus[%d]", i);
		TEST_NEAR_REL(indexstr, test_rhos1_minus[i], sample_rhos1_minus[i], eps);
	}
	for (int i = 0; i < sample_ts_len; i++) {
		snprintf(indexstr, 128, "test_rhos1_plus[%d]", i);
		TEST_NEAR_REL(indexstr, test_rhos1_plus[i] , sample_rhos1_plus[i], eps);
	}
	for (int i = 0; i < sample_ts_len; i++) {
		snprintf(indexstr, 128, "test_rhos2[%d]", i);
		TEST_NEAR_REL(indexstr, test_rhos2[i], sample_rhos2[i], eps);
	}
	for (int i = 0; i < sample_ts_len; i++) {
		snprintf(indexstr, 128, "test_rhos2_minus[%d]", i);
		TEST_NEAR_REL(indexstr, test_rhos2_minus[i], sample_rhos2_minus[i], eps);
	}
	for (int i = 0; i < sample_ts_len; i++) {
		snprintf(indexstr, 128, "test_rhos2_plus[%d]", i);
		TEST_NEAR_REL(indexstr, test_rhos2_plus[i], sample_rhos2_plus[i], eps);
	}
	for (int i = 0; i < sample_ts_len; i++) {
		snprintf(indexstr, 128, "test_ts[%d]", i);
		TEST_NEAR_REL(indexstr, test_ts[i], sample_ts[i], eps);
	}
	TEST_EQUAL("ts_len", test_ts_len, sample_ts_len);
}

static void
test_pose_from_point_tangents_2_fn_t_for_root()
{
	pose_poly<double> p = {
		sample_A0, sample_A1, sample_A2,
		sample_B0, sample_B1, sample_B2, sample_B3,
		sample_C0, sample_C1, sample_C2, sample_C3, sample_C4,
		sample_E0, sample_E1, sample_E2,
		sample_F0, sample_F1, sample_F2, sample_F3,
		sample_G0, sample_G1, sample_G2, sample_G3, sample_G4,
		sample_H0, sample_H1, sample_H2, sample_H3, sample_H4,
		sample_J0, sample_J1, sample_J2, sample_J3,
		sample_K0, sample_K1, sample_K2, sample_K3,
		sample_L0, sample_L1, sample_L2
	};

	double output[T_VECTOR_LEN][11];
	for (int i = 0; i < T_VECTOR_LEN; i++) {
		p.pose_from_point_tangents_2_fn_t_for_root(sample_t_vector[i], &output[i]);
	}

	for (int i = 0; i < T_VECTOR_LEN; i++) {
		double test_fvalue = output[i][0];
		char indexstr[128];
		snprintf(indexstr, 128, "fvalue[%d]", i);
		TEST_NEAR_REL(indexstr, test_fvalue, sample_fvalue_pose_poly[i], eps);
	}
}

static void
test_get_sigmas()
{
	pose_poly<double> p = {
		sample_A0, sample_A1, sample_A2,
		sample_B0, sample_B1, sample_B2, sample_B3,
		sample_C0, sample_C1, sample_C2, sample_C3, sample_C4,
		sample_E0, sample_E1, sample_E2,
		sample_F0, sample_F1, sample_F2, sample_F3,
		sample_G0, sample_G1, sample_G2, sample_G3, sample_G4,
		sample_H0, sample_H1, sample_H2, sample_H3, sample_H4,
		sample_J0, sample_J1, sample_J2, sample_J3,
		sample_K0, sample_K1, sample_K2, sample_K3,
		sample_L0, sample_L1, sample_L2
	};

	bool pass;
	double output[4][TS_MAX_LEN][TS_MAX_LEN];
	char indexstr[128];

	double (&test_sigmas1)[TS_MAX_LEN][TS_MAX_LEN] = output[0];
	double (&test_sigmas2)[TS_MAX_LEN][TS_MAX_LEN] = output[1];
	double (&test_sigmas1_end)[TS_MAX_LEN]        = output[2][0];
	double (&test_sigmas2_end)[TS_MAX_LEN]        = output[3][0];

	p.get_sigmas(sample_ts_len, sample_ts, &output);

	pass = true;
	for (int i = 0; i < sample_ts_len; i++) {
		snprintf(indexstr, 128, "test_sigmas1_end[%d]", i);
		TEST_EQUAL(indexstr, test_sigmas1_end[i], sample_sigmas1_end[i]);

		pass = test_sigmas1_end[i] == sample_sigmas1_end[i] && pass;
		for (int j = 0; pass && j < sample_sigmas1_end[i]; j++) {
			snprintf(indexstr, 128, "test_sigmas1[%d][%d]", i, j);
			TEST_NEAR_REL(indexstr, test_sigmas1[i][j], sample_sigmas1[i][j], eps);
		}
		std::cout << std::endl;
	}

	pass = true;
	for (int i = 0; i < sample_ts_len; i++) {
		snprintf(indexstr, 128, "test_sigmas2_end[%d]", i);
		TEST_EQUAL(indexstr, test_sigmas2_end[i], sample_sigmas2_end[i]);

		pass = test_sigmas2_end[i] == sample_sigmas2_end[i] && pass;
		for (int j = 0; pass && j < sample_sigmas2_end[i]; j++) {
			snprintf(indexstr, 128, "test_sigmas2[%d][%d]", i, j);
			TEST_NEAR_REL(indexstr, test_sigmas2[i][j], sample_sigmas2[i][j], eps);
		}
		std::cout << std::endl;
	}
}

static void
test_pose_from_point_tangents_2_fn_t()
{
	pose_poly<double> p = {
		sample_A0, sample_A1, sample_A2,
		sample_B0, sample_B1, sample_B2, sample_B3,
		sample_C0, sample_C1, sample_C2, sample_C3, sample_C4,
		sample_E0, sample_E1, sample_E2,
		sample_F0, sample_F1, sample_F2, sample_F3,
		sample_G0, sample_G1, sample_G2, sample_G3, sample_G4,
		sample_H0, sample_H1, sample_H2, sample_H3, sample_H4,
		sample_J0, sample_J1, sample_J2, sample_J3,
		sample_K0, sample_K1, sample_K2, sample_K3,
		sample_L0, sample_L1, sample_L2
	};

	double output[T_VECTOR_LEN][11];
	for (int i = 0; i < T_VECTOR_LEN; i++) {
		p.pose_from_point_tangents_2_fn_t(sample_t_vector[i], &output[i]);
	}

	for (int i = 0; i < T_VECTOR_LEN; i++) {
		double test_fvalue = output[i][0];
		char indexstr[128];
		snprintf(indexstr, 128, "fvalue[%d]", i);
		TEST_NEAR_REL(indexstr, test_fvalue, sample_fvalue_pose_poly[i], eps);
	}
}

static void
test_get_r_t_from_rhos()
{
	pose_poly<double> p;

	double RT[RT_MAX_LEN + 1][4][3];
	char indexstr[128];

	p.get_r_t_from_rhos(
		sample_ts_len,
		sample_sigmas1, sample_sigmas1_end,
		sample_sigmas2, sample_sigmas2_end,
		sample_rhos1, sample_rhos2,
		sample_gama1, sample_tgt1,
		sample_gama2, sample_tgt2,
		sample_Gama1, sample_Tgt1,
		sample_Gama2, sample_Tgt2,
		&RT
	);

	double (&RT_len) = RT[0][0][0];
	for (int i = 0; i < RT_len; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				double (&test_Rots)[4][3] = RT[i+1];
				snprintf(indexstr, 128, "RT[%d]: Rots[%d][%d]", i+1, j, k);
				TEST_NEAR_REL(indexstr, test_Rots[j][k], sample_Rots[i][j][k], eps);
			}
		}
		for (int j = 0; j < 3; j++) {
			double (&test_Transls)[3] = RT[i+1][3];
			snprintf(indexstr, 128, "RT[%d]: Transls[%d]", i+1, j);
			TEST_NEAR_REL(indexstr, test_Transls[j], sample_Transls[i][j], eps);
		}
		std::cout << std::endl;
	}

}

// main test function - place all the tests to be run here
void
test_p2pt()
{
	std::cout << "*** USING eps = " << std::scientific << eps << " FOR ALL TESTS ***" << std::endl;

	std::cout << "\nTEST #1 - test_hello" << std::endl;
	test_hello();

	std::cout << "\nTEST #2 - test_pose_from_tangents_root_find_function_any" << std::endl;
	test_pose_from_point_tangents_root_find_function_any();

	std::cout << "\nTEST #3 - test_pose_from_point_tangents_2" << std::endl;
	test_pose_from_point_tangents_2();

	std::cout << "\nTEST #4 - test_find_bounded_root_intervals" << std::endl;
	test_find_bounded_root_intervals();

	std::cout << "\nTEST #5 - test_sample_pose_poly" << std::endl;
	test_sample_pose_poly();

	std::cout << "\nTEST #6 - test_rhos_from_root_ids" << std::endl;
	test_rhos_from_root_ids();

	std::cout << "\nTEST #7 - test_pose_from_point_tangents_2_fn_t_for_root" << std::endl;
	test_pose_from_point_tangents_2_fn_t_for_root();

	std::cout << "\nTEST #8 - test_get_sigmas" << std::endl;
	test_get_sigmas();

	std::cout << "\nTEST #9 - test_pose_from_point_tangents_2_fn_t" << std::endl;
	test_pose_from_point_tangents_2_fn_t();

	std::cout << "\nTEST #10 - test_get_r_t_from_rhos" << std::endl;
	test_get_r_t_from_rhos();
}

TESTMAIN(test_p2pt);

