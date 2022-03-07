#include <iostream>
#include "../diffgeom2pose/common.hxx"

// %% Check 2 x 1 1 / 2 points SRS, following Fabbri et al. 2019
// %% Special thanks to anonymous reviewer for this short demo
// %% Simply runs on completely random data.


// %% Functions
// From vector to skew matrix
template<typename T>
T** v2skew(T* v)
{
	T** output = (T**)malloc(9 * sizeof(T*));
	output[0] =  0;
	output[1] = -v[2];
	output[2] =  v[1];
	output[3] =  v[2];
	output[4] =  0;
	output[5] = -v[0];
	output[6] = -v[1];
	output[7] =  v[0];
	output[8] =  0;
	return output;
}

// From vector to rotation
template<typename T>
T** v2rot(T* v)
{
	return common::expm(v2skew(v)); // TODO: implement `expm()`
}

// From skew matrix to vector
template<typename T>
T* skew2v(T** S)
{
	T* output = (T*)malloc(3 * sizeof(T));
	output[0] = S(2, 1);
	output[1] = S(0, 2);
	output[2] = S(1, 0);
	return output;
}

int main()
{
	std::cout
		<< "--------------------------------------" << '\n'
		<< "--- demo point pairs with tangents ---" << '\n'
		<< "--------------------------------------" << std::endl;

	double** R_tilde = v2rot(common::randn(3,1));
	double** T_tilde;
	return 0;
}