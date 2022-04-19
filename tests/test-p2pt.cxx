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

#include <p2pt/common.hxx>
#include <p2pt/rf_pose_from_point_tangents_2.hxx>
#include <p2pt/rf_find_bounded_root_intervals.hxx>
#include <p2pt/rf_rhos_from_root_ids.hxx>
#include <p2pt/rf_get_sigmas.hxx>

#include <tests/test-p2pt-constants.hxx>

// TODO: Make this global across files
//constexpr int t_vector_len = 2001;
//constexpr int root_ids_len = t_vector_len - 1;

using namespace P2Pt;

static void
test_hello()
{
	p2pt<double>::hello();
	TEST("Rodou hello? ", true, true);
}

static void
test_rf_pose_from_point_tangents_2()
{
	pose_poly<double> p;

	p.rf_pose_from_point_tangents_2(
		sample_gama1, sample_tgt1,
		sample_gama2, sample_tgt2,
		sample_Gama1, sample_Tgt1,
		sample_Gama2, sample_Tgt2
	);

	TEST_NEAR("A0", p.A0, sample_A0, eps);
	TEST_NEAR("A1", p.A1, sample_A1, eps);
	TEST_NEAR("A2", p.A2, sample_A2, eps);

	TEST_NEAR("B0", p.B0, sample_B0, eps);
	TEST_NEAR("B1", p.B1, sample_B1, eps);
	TEST_NEAR("B2", p.B2, sample_B2, eps);
	TEST_NEAR("B3", p.B3, sample_B3, eps);

	TEST_NEAR("C0", p.C0, sample_C0, eps);
	TEST_NEAR("C1", p.C1, sample_C1, eps);
	TEST_NEAR("C2", p.C2, sample_C2, eps);
	TEST_NEAR("C3", p.C3, sample_C3, eps);
	TEST_NEAR("C4", p.C4, sample_C4, eps);

	TEST_NEAR("E0", p.E0, sample_E0, eps);
	TEST_NEAR("E1", p.E1, sample_E1, eps);
	TEST_NEAR("E2", p.E2, sample_E2, eps);

	TEST_NEAR("F0", p.F0, sample_F0, eps);
	TEST_NEAR("F1", p.F1, sample_F1, eps);
	TEST_NEAR("F2", p.F2, sample_F2, eps);
	TEST_NEAR("F3", p.F3, sample_F3, eps);

	TEST_NEAR("G0", p.G0, sample_G0, eps);
	TEST_NEAR("G1", p.G1, sample_G1, eps);
	TEST_NEAR("G2", p.G2, sample_G2, eps);
	TEST_NEAR("G3", p.G3, sample_G3, eps);
	TEST_NEAR("G4", p.G4, sample_G4, eps);

	TEST_NEAR("H0", p.H0, sample_H0, eps);
	TEST_NEAR("H1", p.H1, sample_H1, eps);
	TEST_NEAR("H2", p.H2, sample_H2, eps);
	TEST_NEAR("H3", p.H3, sample_H3, eps);
	TEST_NEAR("H4", p.H4, sample_H4, eps);

	TEST_NEAR("J0", p.J0, sample_J0, eps);
	TEST_NEAR("J1", p.J1, sample_J1, eps);
	TEST_NEAR("J2", p.J2, sample_J2, eps);
	TEST_NEAR("J3", p.J3, sample_J3, eps);

	TEST_NEAR("K0", p.K0, sample_K0, eps);
	TEST_NEAR("K1", p.K1, sample_K1, eps);
	TEST_NEAR("K2", p.K2, sample_K2, eps);
	TEST_NEAR("K3", p.K3, sample_K3, eps);

	TEST_NEAR("L0", p.L0, sample_L0, eps);
	TEST_NEAR("L1", p.L1, sample_L1, eps);
	TEST_NEAR("L2", p.L2, sample_L2, eps);
}

static void
test_rf_find_bounded_root_intervals()
{
	static double test_root_ids[root_ids_len];

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

	p.rf_find_bounded_root_intervals(sample_t_vector, test_root_ids);

	for (int i = 0; i < root_ids_len; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "root_ids[%d]", i);
		TEST_NEAR(indexstr, test_root_ids[i], sample_root_ids[i], eps);
	}
}

static void
test_rf_sample_pose_poly()
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

	double output[11][t_vector_len];

	p.rf_sample_pose_poly(sample_t_vector, output);

	double *test_fvalue = output[0];

	for (int i = 0; i < t_vector_len; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "fvalue[%d]", i);
		TEST_NEAR(indexstr, test_fvalue[i], sample_fvalue_pose_poly[i], eps);
	}
}

static void
test_rf_rhos_from_root_ids()
{
	pose_poly<double> p;
	p.alpha = sample_alpha;
	p.beta  = sample_beta;
	p.theta = sample_theta;

	double output[7][t_vector_len];

	p.rf_rhos_from_root_ids(sample_t_vector, sample_root_ids, output);

	double* test_rhos1       = output[0];
	double* test_rhos1_minus = output[1];
	double* test_rhos1_plus  = output[2];
	double* test_rhos2       = output[3];
	double* test_rhos2_minus = output[4];
	double* test_rhos2_plus  = output[5];
	double* test_ts          = output[6];

	for (int i = 0; i < 4; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "test_rhos1[%d]", i);
		TEST_NEAR(indexstr, test_rhos1[i], sample_rhos1[i], eps);
	}
	for (int i = 0; i < 4; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "test_rhos1_minus[%d]", i);
		TEST_NEAR(indexstr, test_rhos1_minus[i], sample_rhos1_minus[i], eps);
	}
	for (int i = 0; i < 4; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "test_rhos1_plus[%d]", i);
		TEST_NEAR(indexstr, test_rhos1_plus[i] , sample_rhos1_plus[i], eps);
	}
	for (int i = 0; i < 4; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "test_rhos2[%d]", i);
		TEST_NEAR(indexstr, test_rhos2[i], sample_rhos2[i], eps);
	}
	for (int i = 0; i < 4; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "test_rhos2_minus[%d]", i);
		TEST_NEAR(indexstr, test_rhos2_minus[i], sample_rhos2_minus[i], eps);
	}
	for (int i = 0; i < 4; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "test_rhos2_plus[%d]", i);
		TEST_NEAR(indexstr, test_rhos2_plus[i], sample_rhos2_plus[i], eps);
	}
	for (int i = 0; i < 4; i++) {
		char indexstr[128];
		snprintf(indexstr, 128, "test_ts[%d]", i);
		TEST_NEAR(indexstr, test_ts[i], sample_ts[i], eps);
	}
}

static void
test_rf_pose_from_point_tangents_2_fn_t_for_root()
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

	double output[t_vector_len][11];
	for (int i = 0; i < t_vector_len; i++) {
		p.rf_pose_from_point_tangents_2_fn_t_for_root(sample_t_vector[i], output[i]);
	}

	for (int i = 0; i < t_vector_len; i++) {
		double test_fvalue = output[i][0];
		char indexstr[128];
		snprintf(indexstr, 128, "fvalue[%d]", i);
		TEST_NEAR(indexstr, test_fvalue, sample_fvalue_pose_poly[i], eps);
	}

	//double test_fvalue;
	//char indexstr[128];

	//p.rf_pose_from_point_tangents_2_fn_t_for_root(0, output[0]);
	//test_fvalue = output[0][0];
	//snprintf(indexstr, 128, "fvalue[%d]", 0);
	//TEST_NEAR(indexstr, test_fvalue, 5159552357.01911, eps);

	//p.rf_pose_from_point_tangents_2_fn_t_for_root(1, output[0]);
	//test_fvalue = output[0][0];
	//snprintf(indexstr, 128, "fvalue[%d]", 0);
	//TEST_NEAR(indexstr, test_fvalue, 1599787.71374685, eps);
}

static void
test_rf_get_sigmas()
{

}

static void
test_rf_pose_from_point_tangents_2_fn_t()
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

	double output[t_vector_len][11];
	for (int i = 0; i < t_vector_len; i++) {
		p.rf_pose_from_point_tangents_2_fn_t(sample_t_vector[i], output[i]);
	}

	for (int i = 0; i < t_vector_len; i++) {
		double test_fvalue = output[i][0];
		char indexstr[128];
		snprintf(indexstr, 128, "fvalue[%d]", i);
		TEST_NEAR(indexstr, test_fvalue, sample_fvalue_pose_poly[i], eps);
	}
}

static void
test_rf_get_r_t_from_rhos()
{

}

// main test function - place all the tests to be run here
void
test_p2pt()
{
	std::cout << "*** USING eps = " << std::scientific << eps << " FOR ALL TESTS ***" << std::endl;

	std::cout << "\nTEST #1 - test_hello" << std::endl;
	test_hello();

	std::cout << "\nTEST #2 - test_rf_pose_from_point_tangents_2" << std::endl;
	test_rf_pose_from_point_tangents_2();

	std::cout << "\nTEST #3 - test_rf_find_bounded_root_intervals" << std::endl;
	test_rf_find_bounded_root_intervals();

	std::cout << "\nTEST #4 - test_rf_sample_pose_poly" << std::endl;
	test_rf_sample_pose_poly();

	std::cout << "\nTEST #5 - test_rf_rhos_from_root_ids" << std::endl;
	test_rf_rhos_from_root_ids();

	std::cout << "\nTEST #6 - test_rf_pose_from_point_tangents_2_fn_t_for_root" << std::endl;
	test_rf_pose_from_point_tangents_2_fn_t_for_root();

	std::cout << "\nTEST #7 - test_rf_get_sigmas" << std::endl;
	test_rf_get_sigmas();

	std::cout << "\nTEST #8 - test_rf_pose_from_point_tangents_2_fn_t" << std::endl;
	test_rf_pose_from_point_tangents_2_fn_t();

	std::cout << "\nTEST #9 - test_rf_get_r_t_from_rhos" << std::endl;
	test_rf_get_r_t_from_rhos();
}

TESTMAIN(test_p2pt);

