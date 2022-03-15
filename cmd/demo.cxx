// %% Check 2 x 1 1 / 2 points SRS, following Fabbri et al. 2019
// %% Special thanks to anonymous reviewer for this short demo
// %% Simply runs on completely random data.

#include <iostream>
#include "../diffgeom2pose/common.hxx"
#include "../tests/test.cxx"
#include "../diffgeom2pose/rf_pose_from_point_tangents_root_find_function_any.hxx"

#define MAGIC_NUM 100

int main()
{
	std::cout
		<< "--------------------------------------" << '\n'
		<< "--- demo point pairs with tangents ---" << '\n'
		<< "--------------------------------------" << std::endl;

	double Rots[MAGIC_NUM];
	double Transls[MAGIC_NUM];
	double degen[MAGIC_NUM];

	double* result[] = {Rots, Transls, degen};

	rf_pose_from_point_tangents_root_find_function_any(gama1, tgt1, gama2, tgt2, Gama1, Tgt1, Gama2, Tgt2);

	return 0;
}