// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019
//
// Tests more comprehensive runs of p2pt using the public interface
// 
#include <cstring>
#include <iostream>
#include <p2pt/p2pt.h>
#include <testlib/testlib_test.h>
#include <p2pt/poly.h>

#include <p2pt/common.hxx>
#include <p2pt/rf_pose_from_point_tangents_2.hxx>
#include <p2pt/rf_find_bounded_root_intervals.hxx>

#include <tests/test-p2pt-constants.hxx>

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
		truth_gama1, truth_tgt1,
		truth_gama2, truth_tgt2,
		truth_Gama1, truth_Tgt1,
		truth_Gama2, truth_Tgt2
	);

	TEST_NEAR("A0", p.A0, truth_A0, eps);
	TEST_NEAR("A1", p.A1, truth_A1, eps);
	TEST_NEAR("A2", p.A2, truth_A2, eps);

	TEST_NEAR("B0", p.B0, truth_B0, eps);
	TEST_NEAR("B1", p.B1, truth_B1, eps);
	TEST_NEAR("B2", p.B2, truth_B2, eps);
	TEST_NEAR("B3", p.B3, truth_B3, eps);

	TEST_NEAR("C0", p.C0, truth_C0, eps);
	TEST_NEAR("C1", p.C1, truth_C1, eps);
	TEST_NEAR("C2", p.C2, truth_C2, eps);
	TEST_NEAR("C3", p.C3, truth_C3, eps);
	TEST_NEAR("C4", p.C4, truth_C4, eps);

	TEST_NEAR("E0", p.E0, truth_E0, eps);
	TEST_NEAR("E1", p.E1, truth_E1, eps);
	TEST_NEAR("E2", p.E2, truth_E2, eps);

	TEST_NEAR("F0", p.F0, truth_F0, eps);
	TEST_NEAR("F1", p.F1, truth_F1, eps);
	TEST_NEAR("F2", p.F2, truth_F2, eps);
	TEST_NEAR("F3", p.F3, truth_F3, eps);

	TEST_NEAR("G0", p.G0, truth_G0, eps);
	TEST_NEAR("G1", p.G1, truth_G1, eps);
	TEST_NEAR("G2", p.G2, truth_G2, eps);
	TEST_NEAR("G3", p.G3, truth_G3, eps);
	TEST_NEAR("G4", p.G4, truth_G4, eps);

	TEST_NEAR("H0", p.H0, truth_H0, eps);
	TEST_NEAR("H1", p.H1, truth_H1, eps);
	TEST_NEAR("H2", p.H2, truth_H2, eps);
	TEST_NEAR("H3", p.H3, truth_H3, eps);
	TEST_NEAR("H4", p.H4, truth_H4, eps);

	TEST_NEAR("J0", p.J0, truth_J0, eps);
	TEST_NEAR("J1", p.J1, truth_J1, eps);
	TEST_NEAR("J2", p.J2, truth_J2, eps);
	TEST_NEAR("J3", p.J3, truth_J3, eps);

	TEST_NEAR("K0", p.K0, truth_K0, eps);
	TEST_NEAR("K1", p.K1, truth_K1, eps);
	TEST_NEAR("K2", p.K2, truth_K2, eps);
	TEST_NEAR("K3", p.K3, truth_K3, eps);

	TEST_NEAR("L0", p.L0, truth_L0, eps);
	TEST_NEAR("L1", p.L1, truth_L1, eps);
	TEST_NEAR("L2", p.L2, truth_L2, eps);
}

static void
test_rf_find_bounded_root_intervals()
{
	static double test_root_ids[2000];

	pose_poly<double> p = {
		truth_A0, truth_A1, truth_A2,
		truth_B0, truth_B1, truth_B2, truth_B3,
		truth_C0, truth_C1, truth_C2, truth_C3, truth_C4,
		truth_E0, truth_E1, truth_E2,
		truth_F0, truth_F1, truth_F2, truth_F3,
		truth_G0, truth_G1, truth_G2, truth_G3, truth_G4,
		truth_H0, truth_H1, truth_H2, truth_H3, truth_H4,
		truth_J0, truth_J1, truth_J2, truth_J3,
		truth_K0, truth_K1, truth_K2, truth_K3,
		truth_L0, truth_L1, truth_L2
	};

	p.rf_find_bounded_root_intervals(truth_t_vector, test_root_ids);
	
	for (int i = 0; i < 2000; i++) {
		char indexstr[20];
		snprintf(indexstr, 20, "root_ids[%d]", i);
		TEST_NEAR(indexstr, test_root_ids[i], truth_root_ids[i], eps);
	}
}

// main test function - place all the tests to be run here
void
test_p2pt()
{
	test_hello();
	test_rf_pose_from_point_tangents_2();
	test_rf_find_bounded_root_intervals();
}

TESTMAIN(test_p2pt);

