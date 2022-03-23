#include <iostream>
#include "../diffgeom2pose/rf_pose_from_point_tangents_root_find_function_any.hxx"
#include "../diffgeom2pose/common.hxx"

#define MAGIC_NUM 100

/*****************************************************************************
 * Test variables 
 * 
 * These values were taken from a single MATLAB run and are generated randomly
 *****************************************************************************/

// TODO: Check if `T_tilde`, `Gama1`, `Gama2`, `Tgt1`, `Tgt2`, `gama1`, `gama2`, `tgt1`, and `tgt2` need to be transposed

#pragma region variables
static double R_tilde[3][3] = { 
	{ -0.917851347078650,   0.362140358287978,  -0.162490816863482 },
	{  0.084622149462850,  -0.221430418224087,  -0.971497638548542 },
	{ -0.387798912435548,  -0.905440738416469,   0.172595112125592 } 
};
static double T_tilde[3] = { 0.862173320368121, 0.318765239858981, 8.692311703694727 };

static double X[2][3] = {
	{ -0.433592022305684,   3.578396939725760,  -1.349886940156521 },
	{  0.342624466538650,   2.769437029884877,   3.034923466331855 }
};
static double T[2][3] = {
	{  0.707087314733721,   0.696695198637593,  -0.121009625807131 },
	{ -0.041895437077509,  -0.136185233022273,   0.989797128031171 }
};

static double Gama1[3] = { -0.433592022305684, 3.578396939725760, -1.349886940156521 };
static double Gama2[3] = { 0.342624466538650, 2.769437029884877, 3.034923466331855 };
static double Tgt1[3] = { 0.707087314733721, 0.696695198637593, -0.121009625807131 };
static double Tgt2[3] = { -0.041895437077509, -0.136185233022273, 0.989797128031171 };
static double espi = 0.1;

static double Y[4][3] = {
	{ -0.433592022305684,   3.578396939725760,  -1.349886940156521 },
	{  0.342624466538650,   2.769437029884877,   3.034923466331855 },
	{ -0.362883290832312,   3.648066459589520,  -1.361987902737234 },
	{  0.338434922830899,   2.755818506582650,   3.133903179134971 }
};
static double y[4][3] = {
	{  2.775372523552219,   0.801119794675837,   5.387447963139889 },
	{  1.057472717309608,  -3.213899721813590,   6.575694154824491 },
	{  2.737668859298950,   0.803432425169329,   5.294857005470810 },
	{  1.040302984351099,  -3.307397251949443,   6.606733035742282 }

};
static double gama[4][3] = {
	{  0.515155328188950,   0.148701166147122,   1.000000000000000 },
	{  0.160815374378955,  -0.488754441149852,   1.000000000000000 },
	{  0.517043020514115,   0.151738266838783,   1.000000000000000 },
	{  0.157461029335238,  -0.500610094892059,   1.000000000000000 }
};

static double gama1[3] = { 0.515155328188950, 0.148701166147122, 1.000000000000000 };
static double gama2[3] = { 0.160815374378955, -0.488754441149852, 1.000000000000000 };
static double tgt1[3] = { 0.527886693031222, 0.849314805782026, 0 };
static double tgt2[3] = { -0.272245168540452, -0.962227919053683, 0 };
#pragma endregion

int main()
{
	std::cout
		<< "--------------------------------------" << '\n'
		<< "--- demo point pairs with tangents ---" << '\n'
		<< "--------------------------------------" << std::endl;

	static double Rots[MAGIC_NUM];
	static double Transls[MAGIC_NUM];
	static double degen[MAGIC_NUM];

	static double* result[3] = {Rots, Transls, degen};

	rf_pose_from_point_tangents_root_find_function_any(
		gama1, tgt1,
		gama2, tgt2,
		Gama1, Tgt1,
		Gama2, Tgt2,
		result
	);

	return 0;
}