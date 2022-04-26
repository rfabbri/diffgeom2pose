#include <iostream>
#include <chrono>

#include "p2pt/p2pt.h"
#include "p2pt/poly.h"
#include "tests/test-p2pt-constants.hxx"
#include "p2pt/pose_from_point_tangents_root_find_function_any.hxx"

#define TEST // Disable stdout for testing speed

int main()
{
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	double test_RT[RT_MAX_LEN][4][3] = {0};
	int test_RT_len;
	double test_degen;

	P2Pt::pose_from_point_tangents_root_find_function_any(
		sample_gama1, sample_tgt1,
		sample_gama2, sample_tgt2,
		sample_Gama1, sample_Tgt1,
		sample_Gama2, sample_Tgt2,
		&test_RT, &test_RT_len, &test_degen
	);

	double diff;
	char indexstr[128];
	for (int i = 0; i < sample_RT_len; i++) {
		// Rots
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				double (&test_Rots)[4][3] = test_RT[i];
#ifndef TEST
				snprintf(
					indexstr, 128,
					"RT[%d]: Rots[%d][%d] = %.15g,   \tSample = %.15g,\tdiff = %e\t| %s",
					i, j, k, test_Rots[j][k], sample_Rots[i][j][k],
					(diff = std::abs(test_Rots[j][k] - sample_Rots[i][j][k])),
					diff < eps ? "OK" : "FAILED"
				);
				std::cout << indexstr << std::endl;
#endif
			}
		}
		// Transls
		for (int j = 0; j < 3; j++) {
			double (&test_Transls)[3] = test_RT[i][3];
#ifndef TEST
			snprintf(
				indexstr, 128,
				"RT[%d]: Transls[%d] = %.15g,   \tSample = %.15g,\tdiff = %e\t| %s",
				i, j, test_Transls[j], sample_Transls[i][j],
				(diff = std::abs(test_Transls[j] - sample_Transls[i][j])),
				diff < eps ? "OK" : "FAILED"
			);
			std::cout << indexstr << std::endl;
#endif
		}
#ifndef TEST
		std::cout << std::endl;
#endif
	}
	// Degen
#ifndef TEST
	snprintf(
		indexstr, 128,
		"degen = %.15g,   \tSample = %.15g,\tdiff = %e\t| %s",
		test_degen, sample_degen,
		(diff = std::abs(test_degen - sample_degen)),
		diff < eps ? "OK" : "FAILED"
	);
	std::cout << indexstr << std::endl;
#endif
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "Time of solver: " << duration << "us" << std::endl;
    return 0;
}
