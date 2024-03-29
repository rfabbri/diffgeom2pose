#include <iostream>
#include <chrono>

#include "p2pt/p2pt.h"
#include "tests/test-p2pt-constants.hxx"
#include <cstring>

using namespace P2Pt;

void
test_run(const char *type, const char *benchmark = "no")
{
	double test_RT[RT_MAX_LEN][4][3] = {0};
	double test_degen;
	float  test_RT_float[RT_MAX_LEN][4][3] = {0};
	float  test_degen_float;
	unsigned    test_RT_len;

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	if (strcmp(type, "float") == 0)
		p2pt<float>::pose_from_point_tangents(
			sample_gama1_float, sample_tgt1_float,
			sample_gama2_float, sample_tgt2_float,
			sample_Gama1_float, sample_Tgt1_float,
			sample_Gama2_float, sample_Tgt2_float,
			&test_RT_float, &test_RT_len, &test_degen_float
		);
	else
		p2pt<double>::pose_from_point_tangents(
			sample_gama1, sample_tgt1,
			sample_gama2, sample_tgt2,
			sample_Gama1, sample_Tgt1,
			sample_Gama2, sample_Tgt2,
			&test_RT, &test_RT_len, &test_degen
		);
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

	double diff;
	char indexstr[128];
	bool is_benchmark = strcmp(benchmark, "benchmark") == 0 ? true : false;

	// RT
	for (unsigned i = 0; i < sample_RT_len; i++) {
		// Rots
		for (unsigned j = 0; j < 3; j++) {
			for (unsigned k = 0; k < 3; k++) {
				double (&test_Rots)[4][3] = test_RT[i];
				if (!is_benchmark) {
          diff = std::abs(test_Rots[j][k] - sample_Rots[i][j][k]);
					snprintf(
						indexstr, 128,
						"RT[%d]: Rots[%d][%d] = %.15g,   \tSample = %.15g,\tdiff = %e\t| %s",
						i, j, k, test_Rots[j][k], sample_Rots[i][j][k],
						diff,
						diff < eps ? "OK" : "FAILED"
					);
					std::cout << indexstr << std::endl;
				}
			}
		}
		// Transls
		for (unsigned j = 0; j < 3; j++) {
			double (&test_Transls)[3] = test_RT[i][3];
			if (!is_benchmark) {
			  diff = std::abs(test_Transls[j] - sample_Transls[i][j]);
				snprintf(
					indexstr, 128,
					"RT[%d]: Transls[%d] = %.15g,   \tSample = %.15g,\tdiff = %e\t| %s",
					i, j, test_Transls[j], sample_Transls[i][j],
					diff,
					diff < eps ? "OK" : "FAILED"
				);
				std::cout << indexstr << std::endl;
			}
		}
		if (!is_benchmark) std::cout << std::endl;
	}
	// Degen
	if (!is_benchmark) {
	  diff = std::abs(test_degen - sample_degen);
		snprintf(
			indexstr, 128,
			"degen = %.15g,   \tSample = %.15g,\tdiff = %e\t| %s",
			test_degen, sample_degen,
			diff,
			diff < eps ? "OK" : "FAILED"
		);
		std::cout << indexstr << std::endl;
	}
  std::cout << "Time of solver: " << duration << "us" << std::endl;
}

int main(int argc, const char* argv[])
{
	switch (argc) {
	case 1:
		test_run("double");
		return 0;
	case 2:
		if      (strcmp(argv[1], "double")     == 0) test_run("double");
		else if (strcmp(argv[1], "float")      == 0) test_run("float");
		else if (strcmp(argv[1], "-benchmark") == 0) test_run("double", "benchmark");
		else break;
		return 0;
	case 3:
		if      ((strcmp(argv[1], "double") == 0) && (strcmp(argv[2], "-benchmark") == 0))
			test_run("double", "benchmark");
		else if ((strcmp(argv[1], "float")  == 0) && (strcmp(argv[2], "-benchmark") == 0))
			test_run("float", "benchmark");
		else break;
		return 0;
	}
	std::cout << "usage: ./cmd [-benchmark]\n"
				 "       ./cmd [double | float] [-benchmark]" << std::endl;
    return 1;
}
