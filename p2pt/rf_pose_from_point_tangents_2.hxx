//% This is usually called from rf_pose_from_point_tangents_root_find_function_any.m
//%
//% Input: gama1, gama2, tgt1, tgt2, Gama1, Gama2, Tgt1, Tgt2
//% ----------------------------------------------------------------------------
//
//% dummy input ---------------------------
//%gama1 = rand(1,3);
//%gama1(3) = 1;
//%gama2 = rand(1,3);
//%gama2(3) = 1;
//%tgt1 = rand(1,3);
//%tgt1(3) = 0;
//%tgt2 = rand(1,3);
//%tgt2(3) = 0;
//
//%Gama1 = rand(1,3);
//%Gama2 = rand(1,3);
//%Tgt1 = rand(1,3);
//%Tgt1 = Tgt1/norm(Tgt1);
//%Tgt2 = rand(1,3);
//%Tgt2 = Tgt2/norm(Tgt2);
//
//%t=1;
//
//% !dummy input ---------------------------
//
//
//% Precomputed, t-independent terms:


constexpr double PI = 3.141592653589793;


// TODO: Check if these values need a template and need to be global
double A0, A1, A2;
double B0, B1, B2, B3;
double C0, C1, C2, C3, C4;
double E0, E1, E2;
double F0, F1, F2, F3;
double G0, G1, G2, G3, G4;
double H0, H1, H2, H3, H4;
double J0, J1, J2, J3;
double K0, K1, K2, K3;
double L0, L1, L2;
double alpha, beta, theta;

double g11 = gama1[0];
double g12 = gama1[1];
double g21 = gama2[0];
double g22 = gama2[1];
double h11 = tgt1[0];
double h12 = tgt1[1];
double h21 = tgt2[0];
double h22 = tgt2[1];

double V[3];
common::vec_sub(Gama1, Gama2, V);


static double buff[3];

vec1vec2_el_wise_mult(V, V, buff));
double a1 = vec_sum(buff);

vec1vec2_el_wise_mult(Tgt1, Tgt1, buff));
double a2 = vec_sum(buff);

vec1vec2_el_wise_mult(Tgt2, Tgt2, buff));
double a3 = vec_sum(buff);

vec1vec2_el_wise_mult(V, Tgt1, buff));
double a4 = vec_sum(buff);

vec1vec2_el_wise_mult(Tgt1, Tgt2, buff));
double a5 = vec_sum(buff);

vec1vec2_el_wise_mult(V, Tgt2, buff));
double a6 = vec_sum(buff);


double t4 = g11 * g11;
double t5 = g12 * g12;
double t6 = g21 * g21;
double t7 = g22 * g22;
double t11 = 2 * (1 + g11 * g21 + g12 * g22) / (t4 + t5 - t6 - t7); 

double theta = 0.5*atan(t11);
if (theta < 0) {
	theta = theta + PI / 2;
}

//% 497-798
//%theta = .7865071740;
//
//% 240-1100
//%theta = .7887953040;
//
//% 101-406: 
//%theta = .7852237735
//
//% the above theta has similar sin(2theta), cos(2theta) as the maple spreadsheet.

double t1 = sin(2*theta);
double t2 = cos(2*theta);
double t5 = g22 * g22;
double t6 = t2 * t2;
double t8 = g21 * g21;
double t14 = g12 * g12;
double t15 = t1 * t1;
double t21 = g11 * g11;

double den1 = 2*t1*(g11*g21 +g12*g22 + 1) + t2*(t21 + t14 - t8 - t5);
double den2 = t21 + t14 + t8 + t5 + 2;

double t25 = -2*a1 / (den1 - den2);

double beta = sqrt(t25);


double t24 = 2*a1 / (den1 + den2);

double alpha = sqrt(t24);


//% Coefficient code adapted from Maple ::: can be further cleaned up but works

double A0 = a4 * a4 * g12 * g12 + a4 * a4 * g11 * g11 + a4 * a4 + 2.0 * a2 * pow(g11, 
3.0) * g21 * beta * beta * sin(theta) * cos(theta) + 2.0 * a2 * g21 * g11 * 
g12 * g12 * beta * beta * sin(theta) * cos(theta) - 2.0 * a2 * g11 * g11 * g12 
* g12 * beta * beta * pow(sin(theta), 2.0) - a2 * pow(g12, 4.0) * beta * 
beta * pow(sin(theta), 2.0) - a2 * g21 * g21 * g11 * g11 * beta * beta * 
pow(cos(theta), 2.0) + 2.0 * a2 * g12 * g12 * beta * beta * sin(theta) * 
cos(theta) + 2.0 * a2 * g11 * g11 * beta * beta * sin(theta) * cos(theta) + 
2.0 * a2 * g11 * g11 * g22 * g12 * beta * beta * sin(theta) * cos(theta) - a2 
* beta * beta * pow(cos(theta), 2.0) + 2.0 * a2 * pow(g12, 3.0) * g22 * 
beta * beta * sin(theta) * cos(theta) - a2 * pow(g11, 4.0) * beta * beta * 
pow(sin(theta), 2.0) - 2.0 * a2 * g11 * g11 * beta * beta * pow(sin(theta), 
2.0) - 2.0 * a2 * g12 * g12 * beta * beta * pow(sin(theta), 2.0) + 2.0 * 
a2 * beta * beta * sin(theta) * cos(theta) - 2.0 * a2 * g21 * g11 * g22 * g12 
* beta * beta * pow(cos(theta), 2.0) - a2 * beta * beta * pow(sin(theta), 
2.0) + 2.0 * a2 * g21 * g11 * beta * beta * sin(theta) * cos(theta) - a2 * 
g22 * g22 * g12 * g12 * beta * beta * pow(cos(theta), 2.0) - 2.0 * a2 * g22 
* g12 * beta * beta * pow(cos(theta), 2.0) - 2.0 * a2 * g21 * g11 * beta * 
beta * pow(cos(theta), 2.0) + 2.0 * a2 * g22 * g12 * beta * beta * 
sin(theta) * cos(theta);

double A1 = 0.4e1 * a2 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g21 * g21 * g11 * g11 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * g12 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.8e1 * a2 * g21 * g11 * alpha * sin(theta) * 
beta * cos(theta) - 0.4e1 * a2 * g12 * g12 * beta * pow(sin(theta), 0.2e1) * 
alpha + 0.4e1 * a2 * g22 * g12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * 
a2 * g22 * g12 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g22 * g22 
* g12 * g12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g11 * 
g12 * g12 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g21 * g11 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.4e1 * a2 * g21 * g11 * g12 * g12 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.4e1 * a2 * g21 * g11 * beta * pow(sin(theta), 0.2e1) * alpha - 
0.8e1 * a2 * g11 * g11 * g12 * g12 * alpha * sin(theta) * beta * cos(theta) - 
0.4e1 * a2 * pow(g11, 0.4e1) * alpha * sin(theta) * beta * cos(theta) - 0.8e1 * 
a2 * g11 * g11 * alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a2 * g21 * g11 
* g22 * g12 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * pow(g12, 
0.3e1) * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * pow(g12, 
0.3e1) * g22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * pow(g12, 
0.4e1) * alpha * sin(theta) * beta * cos(theta) - 0.8e1 * a2 * g12 * g12 * alpha 
* sin(theta) * beta * cos(theta) + 0.8e1 * a2 * g22 * g12 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a2 * pow(g11, 0.3e1) * g21 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.4e1 * a2 * pow(g11, 0.3e1) * g21 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.4e1 * a2 * g11 * g11 * g22 * g12 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * g22 * g12 * beta * pow(sin(theta), 
0.2e1) * alpha;

double A2 =  (2 * a4 * a4 * g12 * g12) +  (2 * a4 * a4 * g11 * g11) + 
 (2 * a4 * a4) + 0.2e1 * a2 *   pow( g12,  
4) * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 *   
pow( g11,  4) * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * 
a2 *  (g11 * g11) * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * 
 (g12 * g12) * beta * beta * pow(sin(theta), 0.2e1) - 0.4e1 * a2 * beta 
* beta * sin(theta) * cos(theta) + 0.2e1 * a2 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a2 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g21 * 
 g11 *  (g12 * g12) * beta * beta * sin(theta) * cos(theta) + 
0.2e1 * a2 * g21 * g21 *  (g11 * g11) * beta * beta * pow(cos(theta), 
0.2e1) - 0.4e1 * a2 *  (g12 * g12) * beta * beta * sin(theta) * 
cos(theta) + 0.4e1 * a2 * g21 *  g11 * beta * beta * pow(cos(theta), 
0.2e1) + 0.2e1 * a2 * g22 * g22 *  (g12 * g12) * beta * beta * 
pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g22 *  g12 * beta * beta * 
sin(theta) * cos(theta) - 0.4e1 * a2 *  (g11 * g11) * beta * beta * 
sin(theta) * cos(theta) - 0.4e1 * a2 * g21 *  g11 * beta * beta * 
sin(theta) * cos(theta) + 0.4e1 * a2 *  (g11 * g11) *  (g12 * 
g12) * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g22 *  g12 * 
beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 *   pow( 
g11,  3) * g21 * beta * beta * sin(theta) * cos(theta) - 0.4e1 * a2 * 
 (g11 * g11) * g22 *  g12 * beta * beta * sin(theta) * 
cos(theta) + 0.4e1 * a2 * g21 *  g11 * g22 *  g12 * beta * beta 
* pow(cos(theta), 0.2e1) - 0.4e1 * a2 *   pow( g12, 
 3) * g22 * beta * beta * sin(theta) * cos(theta) - 0.4e1 * a2 * 
  pow( g11,  4) * alpha * alpha * pow(cos(theta), 
0.2e1) - 0.8e1 * a2 *  (g11 * g11) * alpha * alpha * pow(cos(theta), 
0.2e1) - 0.4e1 * a2 *   pow( g12,  4) * alpha * 
alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a2 *  (g12 * g12) * alpha * 
alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a2 * alpha * alpha * cos(theta) * 
sin(theta) - 0.4e1 * a2 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a2 * 
alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a2 * g22 *  g12 * alpha 
* alpha * cos(theta) * sin(theta) - 0.4e1 * a2 * g21 * g21 *  (g11 * 
g11) * alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a2 *  (g12 * 
g12) * alpha * alpha * cos(theta) * sin(theta) - 0.8e1 * a2 * g21 *  g11 
* alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a2 * g22 * g22 *  
(g12 * g12) * alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a2 * g21 * 
 g11 * alpha * alpha * cos(theta) * sin(theta) - 0.8e1 * a2 *  
(g11 * g11) *  (g12 * g12) * alpha * alpha * pow(cos(theta), 0.2e1) - 
0.8e1 * a2 *  (g11 * g11) * alpha * alpha * cos(theta) * sin(theta) - 
0.8e1 * a2 * g21 *  g11 *  (g12 * g12) * alpha * alpha * 
cos(theta) * sin(theta) - 0.8e1 * a2 * g21 *  g11 * g22 *  g12 * 
alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a2 *   
pow( g12,  3) * g22 * alpha * alpha * cos(theta) * sin(theta) - 
0.8e1 * a2 * g22 *  g12 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 
* a2 *   pow( g11,  3) * g21 * alpha * alpha * 
cos(theta) * sin(theta) - 0.8e1 * a2 *  (g11 * g11) * g22 *  g12 
* alpha * alpha * cos(theta) * sin(theta);

double B0 = -0.2e1 * beta * sin(theta) * (a2 * g21 * g11 * g22 * h12 * beta * beta * 
pow(cos(theta), 0.2e1) + a2 * pow(g12, 0.3e1) * h12 * beta * beta * 
pow(sin(theta), 0.2e1) + a2 * g21 * h11 * g22 * g12 * beta * beta * 
pow(cos(theta), 0.2e1) + a2 * g11 * h11 * g12 * g12 * beta * beta * 
pow(sin(theta), 0.2e1) - a2 * g11 * h11 * beta * beta * sin(theta) * cos(theta) 
- a4 * a4 * h11 * g11 + a2 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) + 
a2 * g22 * g22 * h12 * g12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g22 * 
h12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g12 * h12 * beta * beta * 
pow(sin(theta), 0.2e1) - a2 * g11 * h11 * g22 * g12 * beta * beta * sin(theta) * 
cos(theta) - a2 * g12 * h12 * beta * beta * sin(theta) * cos(theta) - a2 * g21 * 
h11 * g12 * g12 * beta * beta * sin(theta) * cos(theta) - a2 * g21 * h11 * beta 
* beta * sin(theta) * cos(theta) - a2 * g11 * g11 * g22 * h12 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * 
sin(theta) * cos(theta) - a2 * g22 * h12 * beta * beta * sin(theta) * cos(theta) 
+ a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g21 * 
g21 * h11 * g11 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g21 * h11 * beta * 
beta * pow(cos(theta), 0.2e1) - a2 * g21 * g11 * g12 * h12 * beta * beta * 
sin(theta) * cos(theta) + a2 * g11 * g11 * g12 * h12 * beta * beta * 
pow(sin(theta), 0.2e1) - a4 * a4 * h12 * g12);

double B1 = -0.2e1 * beta * sin(theta) * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * 
sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g22 * h12 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g22 * h12 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sin(theta) * beta 
* cos(theta) - 0.2e1 * a2 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.2e1 * a2 * g12 * h12 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * 
g22 * h12 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g11 * h11 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * h11 * g21 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * g21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * h11 * g22 * g12 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * h12 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * 
sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g21 * h11 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g21 * h11 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a2 * g21 * h11 * alpha * sin(theta) * beta * cos(theta) 
- 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sin(theta) * beta * cos(theta) + 
0.4e1 * a2 * g11 * h11 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * 
g21 * g11 * g22 * h12 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a2 * 
g11 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * g21 * h11 
* g11 * alpha * sin(theta) * beta * cos(theta)) - 0.4e1 * alpha * cos(theta) * 
(a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * 
pow(g12, 0.3e1) * h12 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g21 * h11 * 
g22 * g12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g11 * h11 * g12 * g12 * 
beta * beta * pow(sin(theta), 0.2e1) - a2 * g11 * h11 * beta * beta * sin(theta) 
* cos(theta) - a4 * a4 * h11 * g11 + a2 * g11 * h11 * beta * beta * 
pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * g12 * beta * beta * 
pow(cos(theta), 0.2e1) + a2 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) + 
a2 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - a2 * g11 * h11 * g22 * 
g12 * beta * beta * sin(theta) * cos(theta) - a2 * g12 * h12 * beta * beta * 
sin(theta) * cos(theta) - a2 * g21 * h11 * g12 * g12 * beta * beta * sin(theta) 
* cos(theta) - a2 * g21 * h11 * beta * beta * sin(theta) * cos(theta) - a2 * g11 
* g11 * g22 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12 * 
g12 * h12 * g22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * g11 
* h11 * g21 * beta * beta * sin(theta) * cos(theta) - a2 * g22 * h12 * beta * 
beta * sin(theta) * cos(theta) + a2 * pow(g11, 0.3e1) * h11 * beta * beta * 
pow(sin(theta), 0.2e1) + a2 * g21 * g21 * h11 * g11 * beta * beta * 
pow(cos(theta), 0.2e1) + a2 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) - 
a2 * g21 * g11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) + a2 * g11 * 
g11 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a4 * h12 * g12);

double B2 = -0.2e1 * beta * sin(theta) * (0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * 
alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * 
alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g11 * h11 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * g22 * h12 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * alpha * alpha * cos(theta) * 
sin(theta) + 0.8e1 * a2 * g11 * g11 * h11 * g21 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a2 * g11 * h11 * g22 * g12 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * g12 * g12 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a2 * g22 * h12 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g22 * h12 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a2 * g21 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cos(theta) * sin(theta) + 
0.8e1 * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) - 
0.2e1 * a4 * a4 * h11 * g11 - 0.2e1 * a4 * a4 * h12 * g12 - 0.2e1 * a2 * g22 * 
g22 * h12 * g12 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g12 * h12 
* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a2 * g22 * h12 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g21 * g21 * h11 * g11 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g22 * h12 * beta * beta * pow(cos(theta), 
0.2e1) - 0.2e1 * a2 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * 
a2 * g11 * h11 * g12 * g12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * 
g11 * g11 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 * g11 
* g11 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a2 * g11 * 
h11 * g22 * g12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12 * h12 
* beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * h11 * beta * beta * pow(sin(theta), 
0.2e1) - 0.2e1 * a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sin(theta), 
0.2e1) + 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cos(theta), 
0.2e1) + 0.2e1 * a2 * g21 * h11 * beta * beta * sin(theta) * cos(theta) - 0.2e1 
* a2 * pow(g12, 0.3e1) * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 
* g21 * h11 * g22 * g12 * beta * beta * pow(cos(theta), 0.2e1)) - 0.4e1 * alpha 
* cos(theta) * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sin(theta) * beta 
* cos(theta) - 0.2e1 * a2 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.2e1 * a2 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * 
g22 * g22 * h12 * g12 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 * 
g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g12 * h12 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g22 * h12 * alpha * 
sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g11 * h11 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * h12 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * sin(theta) * 
beta * cos(theta) - 0.2e1 * a2 * g21 * h11 * alpha * pow(cos(theta), 0.2e1) * 
beta + 0.2e1 * a2 * g21 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * 
a2 * g21 * h11 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g21 * g11 
* g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * 
g12 * h12 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g11 * g22 * h12 * 
alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g11 * h11 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * 
sin(theta) * beta * cos(theta)) + 0.2e1 * beta * sin(theta) * (a2 * g21 * g11 * 
g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * pow(g12, 0.3e1) * h12 * 
beta * beta * pow(sin(theta), 0.2e1) + a2 * g21 * h11 * g22 * g12 * beta * beta 
* pow(cos(theta), 0.2e1) + a2 * g11 * h11 * g12 * g12 * beta * beta * 
pow(sin(theta), 0.2e1) - a2 * g11 * h11 * beta * beta * sin(theta) * cos(theta) 
- a4 * a4 * h11 * g11 + a2 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) + 
a2 * g22 * g22 * h12 * g12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g22 * 
h12 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g12 * h12 * beta * beta * 
pow(sin(theta), 0.2e1) - a2 * g11 * h11 * g22 * g12 * beta * beta * sin(theta) * 
cos(theta) - a2 * g12 * h12 * beta * beta * sin(theta) * cos(theta) - a2 * g21 * 
h11 * g12 * g12 * beta * beta * sin(theta) * cos(theta) - a2 * g21 * h11 * beta 
* beta * sin(theta) * cos(theta) - a2 * g11 * g11 * g22 * h12 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * 
sin(theta) * cos(theta) - a2 * g22 * h12 * beta * beta * sin(theta) * cos(theta) 
+ a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g21 * 
g21 * h11 * g11 * beta * beta * pow(cos(theta), 0.2e1) + a2 * g21 * h11 * beta * 
beta * pow(cos(theta), 0.2e1) - a2 * g21 * g11 * g12 * h12 * beta * beta * 
sin(theta) * cos(theta) + a2 * g11 * g11 * g12 * h12 * beta * beta * 
pow(sin(theta), 0.2e1) - a4 * a4 * h12 * g12);

double B3 = -0.2e1 * beta * sin(theta) * (-0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * 
sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g22 * h12 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g22 * h12 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * sin(theta) * beta 
* cos(theta) + 0.2e1 * a2 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta - 
0.2e1 * a2 * g12 * h12 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * 
g22 * h12 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g11 * h11 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * g21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * h11 * g21 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g21 * h11 * g22 * g12 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g12 * h12 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * 
sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g21 * h11 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a2 * g21 * h11 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.4e1 * a2 * g21 * h11 * alpha * sin(theta) * beta * cos(theta) 
+ 0.2e1 * a2 * g21 * g11 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta - 
0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * sin(theta) * beta * cos(theta) - 
0.4e1 * a2 * g11 * h11 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * 
g21 * g11 * g22 * h12 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 * 
g11 * h11 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g21 * g21 * h11 
* g11 * alpha * sin(theta) * beta * cos(theta)) - 0.4e1 * alpha * cos(theta) * 
(0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a2 * g11 * g11 * g12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.4e1 * a2 * g11 * h11 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 * 
g11 * g11 * g22 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * 
g11 * h11 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1 * a2 * g11 * g11 * 
h11 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * pow(g12, 
0.3e1) * h12 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * 
alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * 
alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * g22 * g12 * 
alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * alpha * alpha 
* cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * g22 * g12 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g21 * h11 * g12 * g12 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a2 * g22 * h12 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a2 * g22 * g22 * h12 * g12 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g22 * h12 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g21 * g11 * g22 * h12 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a2 * g21 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a2 * g21 * g11 * g12 * h12 * alpha * alpha * cos(theta) * sin(theta) + 
0.8e1 * a2 * g12 * g12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) - 
0.2e1 * a4 * a4 * h11 * g11 - 0.2e1 * a4 * a4 * h12 * g12 - 0.2e1 * a2 * g22 * 
g22 * h12 * g12 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g12 * h12 
* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a2 * g22 * h12 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g21 * g21 * h11 * g11 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g22 * h12 * beta * beta * pow(cos(theta), 
0.2e1) - 0.2e1 * a2 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * 
a2 * g11 * h11 * g12 * g12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * 
g11 * g11 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 * g11 
* g11 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a2 * g11 * 
h11 * g22 * g12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12 * h12 
* beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a2 * g11 * h11 * beta * beta * pow(sin(theta), 
0.2e1) - 0.2e1 * a2 * pow(g11, 0.3e1) * h11 * beta * beta * pow(sin(theta), 
0.2e1) + 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a2 * g21 * g11 * g22 * h12 * beta * beta * pow(cos(theta), 
0.2e1) + 0.2e1 * a2 * g21 * h11 * beta * beta * sin(theta) * cos(theta) - 0.2e1 
* a2 * pow(g12, 0.3e1) * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 
* g21 * h11 * g22 * g12 * beta * beta * pow(cos(theta), 0.2e1)) + 0.2e1 * beta * 
sin(theta) * (0.2e1 * a2 * g11 * g11 * g22 * h12 * beta * pow(sin(theta), 0.2e1) 
* alpha - 0.2e1 * a2 * g11 * g11 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * 
beta - 0.2e1 * a2 * g21 * h11 * g12 * g12 * alpha * pow(cos(theta), 0.2e1) * 
beta + 0.2e1 * a2 * g21 * g11 * g12 * h12 * beta * pow(sin(theta), 0.2e1) * 
alpha - 0.4e1 * a2 * g12 * g12 * h12 * g22 * alpha * pow(cos(theta), 0.2e1) * 
beta + 0.4e1 * a2 * g12 * g12 * h12 * g22 * beta * pow(sin(theta), 0.2e1) * 
alpha + 0.2e1 * a2 * g21 * h11 * g12 * g12 * beta * pow(sin(theta), 0.2e1) * 
alpha - 0.2e1 * a2 * g11 * h11 * g22 * g12 * alpha * pow(cos(theta), 0.2e1) * 
beta + 0.2e1 * a2 * g11 * h11 * g22 * g12 * beta * pow(sin(theta), 0.2e1) * 
alpha + 0.4e1 * a2 * g11 * h11 * g12 * g12 * alpha * sin(theta) * beta * 
cos(theta) - 0.2e1 * a2 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.2e1 * a2 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * 
g22 * g22 * h12 * g12 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 * 
g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a2 * g12 * h12 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g22 * h12 * alpha * 
sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g11 * h11 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * g11 * h11 * g21 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * g21 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * h11 * g22 * g12 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * h12 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a2 * pow(g12, 0.3e1) * h12 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a2 * pow(g11, 0.3e1) * h11 * alpha * sin(theta) * 
beta * cos(theta) - 0.2e1 * a2 * g21 * h11 * alpha * pow(cos(theta), 0.2e1) * 
beta + 0.2e1 * a2 * g21 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * 
a2 * g21 * h11 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a2 * g21 * g11 
* g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * 
g12 * h12 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g11 * g22 * h12 * 
alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a2 * g11 * h11 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g21 * g21 * h11 * g11 * alpha * 
sin(theta) * beta * cos(theta));

double C0 = -beta * beta * pow(sin(theta), 0.2e1) * (-a4 * a4 * h12 * h12 + 0.2e1 * a2 
* g21 * h11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * 
g11 * h11 * h11 * g21 * beta * beta * sin(theta) * cos(theta) - a4 * a4 * h11 * 
h11 - 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sin(theta) * cos(theta) 
+ a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 
* g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) + a2 * g12 * g12 
* h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * h12 
* beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * g12 * h12 * 
beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta 
* beta * sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 * beta * beta * 
pow(sin(theta), 0.2e1));

double C1 = -beta * beta * pow(sin(theta), 0.2e1) * (0.8e1 * a2 * g11 * h11 * g12 * 
h12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g22 * g22 * h12 * 
h12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g21 * h11 * 
h11 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * h11 * g12 * 
h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g12 * h12 * h12 * g22 
* alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g21 * h11 * g12 * h12 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * h11 * h11 * g21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * h11 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * g22 * h12 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * h12 * h12 * g22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.8e1 * a2 * g21 * h11 * g22 * h12 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * g12 * h12 * h12 * 
alpha * sin(theta) * beta * cos(theta)) - 0.4e1 * beta * sin(theta) * alpha * 
cos(theta) * (-a4 * a4 * h12 * h12 + 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * 
beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta 
* sin(theta) * cos(theta) - a4 * a4 * h11 * h11 - 0.2e1 * a2 * g12 * h12 * h12 * 
g22 * beta * beta * sin(theta) * cos(theta) + a2 * g21 * g21 * h11 * h11 * beta 
* beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * 
beta * sin(theta) * cos(theta) + a2 * g12 * g12 * h12 * h12 * beta * beta * 
pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * h12 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * 
sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 * beta * beta * 
pow(sin(theta), 0.2e1));

double C2 = -beta * beta * pow(sin(theta), 0.2e1) * (-0.4e1 * a2 * g11 * h11 * g12 * 
h12 * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * h12 * g22 
* beta * beta * sin(theta) * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * 
beta * beta * sin(theta) * cos(theta) - 0.2e1 * a4 * a4 * h12 * h12 + 0.4e1 * a2 
* g21 * g21 * h11 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2 * 
g21 * h11 * g12 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * 
g12 * g12 * h12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * 
g22 * g22 * h12 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g21 
* h11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g11 * 
h11 * g22 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * 
h11 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.8e1 * a2 * g11 * h11 
* g12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 * 
g22 * h12 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * g11 * 
h11 * h11 * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * h11 
* h11 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a4 * a4 * h11 * h11 + 
0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.8e1 * a2 * g12 * h12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) + 
0.8e1 * a2 * g11 * h11 * h11 * g21 * alpha * alpha * cos(theta) * sin(theta) + 
0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) - 
0.2e1 * a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - 
0.2e1 * a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1)) - 
0.4e1 * beta * sin(theta) * alpha * cos(theta) * (0.8e1 * a2 * g11 * h11 * g12 * 
h12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g22 * g22 * h12 * 
h12 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g21 * h11 * 
h11 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * h11 * g12 * 
h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g12 * h12 * h12 * g22 
* alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g21 * h11 * g12 * h12 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * h11 * h11 * g21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * h11 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * g22 * h12 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * h12 * h12 * g22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.8e1 * a2 * g21 * h11 * g22 * h12 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * g12 * h12 * h12 * 
alpha * sin(theta) * beta * cos(theta)) - (-0.2e1 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.4e1 * alpha * alpha * pow(cos(theta), 0.2e1)) * (-a4 
* a4 * h12 * h12 + 0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * 
sin(theta) * cos(theta) - a4 * a4 * h11 * h11 - 0.2e1 * a2 * g12 * h12 * h12 * 
g22 * beta * beta * sin(theta) * cos(theta) + a2 * g21 * g21 * h11 * h11 * beta 
* beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * 
beta * sin(theta) * cos(theta) + a2 * g12 * g12 * h12 * h12 * beta * beta * 
pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * h12 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta * beta * 
sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 * beta * beta * 
pow(sin(theta), 0.2e1));

double C3 = -beta * beta * pow(sin(theta), 0.2e1) * (0.4e1 * a2 * g21 * h11 * g12 * h12 
* alpha * pow(cos(theta), 0.2e1) * beta - 0.8e1 * a2 * g11 * h11 * g12 * h12 * 
alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a2 * g21 * h11 * g22 * h12 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g12 * h12 * h12 * g22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g12 * g12 * h12 * h12 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g12 * h12 * h12 * g22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g21 * g21 * h11 * h11 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g11 * h11 * h11 * g21 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g11 * h11 * g22 * h12 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g21 * h11 * g12 * h12 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * g11 * h11 * h11 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g22 * g22 * h12 * h12 * 
alpha * sin(theta) * beta * cos(theta)) - 0.4e1 * beta * sin(theta) * alpha * 
cos(theta) * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * 
sin(theta) * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a4 * a4 * h12 * h12 + 0.4e1 * a2 * g21 * g21 * 
h11 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 * 
g12 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g12 * g12 * 
h12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g22 * g22 * 
h12 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g21 * h11 * g22 
* h12 * beta * beta * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g11 * h11 * g22 * 
h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * g22 * 
h12 * beta * beta * sin(theta) * cos(theta) + 0.8e1 * a2 * g11 * h11 * g12 * h12 
* alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 * g22 * h12 * 
alpha * alpha * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * g11 * h11 * h11 * 
beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * h11 * h11 * 
alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a4 * a4 * h11 * h11 + 0.4e1 * 
a2 * g22 * g22 * h12 * h12 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2 
* g12 * h12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1 * a2 * 
g11 * h11 * h11 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * 
g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12 
* g12 * h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g21 * 
g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1)) - (-0.2e1 * beta * beta 
* pow(sin(theta), 0.2e1) + 0.4e1 * alpha * alpha * pow(cos(theta), 0.2e1)) * 
(0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sin(theta) * beta * cos(theta) - 
0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sin(theta) * beta * cos(theta) - 
0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sin(theta) * beta * cos(theta) - 
0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta - 
0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sin(theta), 0.2e1) * alpha - 
0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cos(theta), 0.2e1) * beta - 
0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sin(theta) * beta * cos(theta) + 
0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha + 
0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sin(theta), 0.2e1) * alpha - 
0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sin(theta) * beta * cos(theta) + 
0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sin(theta), 0.2e1) * alpha + 
0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sin(theta) * beta * cos(theta)) + 
0.4e1 * beta * sin(theta) * alpha * cos(theta) * (-a4 * a4 * h12 * h12 + 0.2e1 * 
a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * 
g11 * h11 * h11 * g21 * beta * beta * sin(theta) * cos(theta) - a4 * a4 * h11 * 
h11 - 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sin(theta) * cos(theta) 
+ a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a2 
* g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) + a2 * g12 * g12 
* h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g22 * g22 * h12 * h12 
* beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 * g12 * h12 * 
beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * g22 * h12 * beta 
* beta * sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 * beta * beta * 
pow(sin(theta), 0.2e1));

double C4 = -0.2e1 * beta * beta * pow(sin(theta), 0.2e1) * (-a4 * a4 * h12 * h12 + 
0.2e1 * a2 * g21 * h11 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 
0.2e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * sin(theta) * cos(theta) - a4 
* a4 * h11 * h11 - 0.2e1 * a2 * g12 * h12 * h12 * g22 * beta * beta * sin(theta) 
* cos(theta) + a2 * g21 * g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1) 
- 0.2e1 * a2 * g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) + 
a2 * g12 * g12 * h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a2 * g22 * 
g22 * h12 * h12 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a2 * g11 * h11 
* g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * h11 * 
g22 * h12 * beta * beta * sin(theta) * cos(theta) + a2 * g11 * g11 * h11 * h11 * 
beta * beta * pow(sin(theta), 0.2e1)) - 0.4e1 * beta * sin(theta) * alpha * 
cos(theta) * (0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * sin(theta) * beta * 
cos(theta) + 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * sin(theta) * beta * 
cos(theta) - 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * pow(sin(theta), 0.2e1) 
* alpha - 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * sin(theta) * beta * 
cos(theta) + 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * pow(cos(theta), 0.2e1) 
* beta + 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * sin(theta) * beta * 
cos(theta) - 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * pow(sin(theta), 0.2e1) 
* alpha + 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * 
beta - 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * pow(sin(theta), 0.2e1) * 
alpha - 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * sin(theta) * beta * 
cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * pow(cos(theta), 0.2e1) 
* beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * 
alpha + 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * sin(theta) * beta * 
cos(theta)) - (-0.2e1 * beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * alpha * 
alpha * pow(cos(theta), 0.2e1)) * (-0.4e1 * a2 * g11 * h11 * g12 * h12 * beta * 
beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * beta 
* sin(theta) * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a4 * a4 * h12 * h12 + 0.4e1 * a2 * g21 * g21 * 
h11 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 * 
g12 * h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g12 * g12 * 
h12 * h12 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a2 * g22 * g22 * 
h12 * h12 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a2 * g21 * h11 * g22 
* h12 * beta * beta * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g11 * h11 * g22 * 
h12 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * g11 * h11 * g22 * 
h12 * beta * beta * sin(theta) * cos(theta) + 0.8e1 * a2 * g11 * h11 * g12 * h12 
* alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a2 * g21 * h11 * g22 * h12 * 
alpha * alpha * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g11 * g11 * h11 * h11 * 
beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a2 * g11 * g11 * h11 * h11 * 
alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a4 * a4 * h11 * h11 + 0.4e1 * 
a2 * g22 * g22 * h12 * h12 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a2 
* g12 * h12 * h12 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1 * a2 * 
g11 * h11 * h11 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a2 * 
g21 * h11 * g12 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a2 * g12 
* g12 * h12 * h12 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a2 * g21 * 
g21 * h11 * h11 * beta * beta * pow(cos(theta), 0.2e1)) + 0.4e1 * beta * 
sin(theta) * alpha * cos(theta) * (0.8e1 * a2 * g11 * h11 * g12 * h12 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g22 * g22 * h12 * h12 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * g21 * h11 * h11 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a2 * g21 * h11 * g12 * h12 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g12 * h12 * h12 * g22 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g21 * h11 * g12 * h12 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a2 * g11 * h11 * h11 * g21 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a2 * g11 * h11 * g22 * h12 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a2 * g11 * g11 * h11 * h11 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * g22 * h12 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * h12 * h12 * g22 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.8e1 * a2 * g21 * h11 * g22 * h12 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a2 * g11 * h11 * h11 * g21 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a2 * g12 * g12 * h12 * h12 * alpha * 
sin(theta) * beta * cos(theta));

double E0 = 0.2e1 * a3 * g21 * g21 * g12 * g22 * beta * beta * cos(theta) * 
sin(theta) + 0.2e1 * a3 * g12 * g22 * beta * beta * cos(theta) * sin(theta) + 
0.2e1 * a3 * g12 * pow(g22, 0.3e1) * beta * beta * cos(theta) * sin(theta) - a3 
* beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g12 * g22 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * g21 * beta * beta * pow(sin(theta), 
0.2e1) - a3 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a3 * g11 * g21 * 
beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g22 * g22 * beta * beta * 
cos(theta) * sin(theta) - a3 * pow(g21, 0.4e1) * beta * beta * pow(cos(theta), 
0.2e1) + 0.2e1 * a3 * g21 * g21 * beta * beta * cos(theta) * sin(theta) - a3 * 
g12 * g12 * g22 * g22 * beta * beta * pow(sin(theta), 0.2e1) - a3 * g11 * g11 * 
g21 * g21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g21 * g21 * g22 
* g22 * beta * beta * pow(cos(theta), 0.2e1) + a6 * a6 * g21 * g21 + a6 * a6 * 
g22 * g22 + 0.2e1 * a3 * g11 * pow(g21, 0.3e1) * beta * beta * cos(theta) * 
sin(theta) + 0.2e1 * a3 * g11 * g21 * g22 * g22 * beta * beta * cos(theta) * 
sin(theta) - 0.2e1 * a3 * g11 * g21 * g12 * g22 * beta * beta * pow(sin(theta), 
0.2e1) + a6 * a6 - 0.2e1 * a3 * g21 * g21 * beta * beta * pow(cos(theta), 0.2e1) 
- 0.2e1 * a3 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) - a3 * pow(g22, 
0.4e1) * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a3 * beta * beta * 
cos(theta) * sin(theta);

double E1 = -0.4e1 * a3 * g11 * g11 * g21 * g21 * alpha * sin(theta) * beta * 
cos(theta) + 0.8e1 * a3 * g21 * g21 * g22 * g22 * alpha * sin(theta) * beta * 
cos(theta) - 0.4e1 * a3 * g12 * g12 * g22 * g22 * alpha * sin(theta) * beta * 
cos(theta) + 0.4e1 * a3 * g22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha - 
0.4e1 * a3 * g22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta - 0.8e1 * a3 * 
g11 * g21 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g12 * g22 * 
alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * g22 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g21 * g21 * g12 * g22 * beta * 
pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g12 * pow(g22, 0.3e1) * alpha * 
pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * pow(g22, 0.3e1) * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.8e1 * a3 * g21 * g21 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a3 * pow(g21, 0.4e1) * alpha * sin(theta) * beta * 
cos(theta) - 0.4e1 * a3 * g11 * g21 * alpha * pow(sin(theta), 0.2e1) * beta + 
0.4e1 * a3 * g11 * g21 * beta * pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 * 
g12 * g22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * 
g12 * g22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * g21 * g22 
* g22 * beta * pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 * g11 * g21 * g12 * 
g22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g11 * pow(g21, 
0.3e1) * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * pow(g21, 
0.3e1) * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g21 * g21 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.8e1 * a3 * g22 * g22 * alpha * sin(theta) * 
beta * cos(theta) - 0.4e1 * a3 * alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * 
a3 * g21 * g21 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * pow(g22, 
0.4e1) * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g11 * g21 * g22 * 
g22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * beta * 
pow(cos(theta), 0.2e1) * alpha;

double E2 = -0.4e1 * a3 * g21 * g21 * g12 * g22 * beta * beta * cos(theta) * sin(theta) 
- 0.4e1 * a3 * g12 * g22 * beta * beta * cos(theta) * sin(theta) - 0.4e1 * a3 * 
g12 * pow(g22, 0.3e1) * beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * 
beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g12 * g22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 * g21 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a3 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a3 * g11 * 
g21 * beta * beta * cos(theta) * sin(theta) - 0.4e1 * a3 * g22 * g22 * beta * 
beta * cos(theta) * sin(theta) + 0.2e1 * a3 * pow(g21, 0.4e1) * beta * beta * 
pow(cos(theta), 0.2e1) - 0.4e1 * a3 * g21 * g21 * beta * beta * cos(theta) * 
sin(theta) + 0.2e1 * a3 * g12 * g12 * g22 * g22 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a3 * g11 * g11 * g21 * g21 * beta * beta * pow(sin(theta), 
0.2e1) + 0.4e1 * a3 * g21 * g21 * g22 * g22 * beta * beta * pow(cos(theta), 
0.2e1) + 0.2e1 * a6 * a6 * g21 * g21 + 0.2e1 * a6 * a6 * g22 * g22 - 0.4e1 * a3 
* g11 * pow(g21, 0.3e1) * beta * beta * cos(theta) * sin(theta) - 0.4e1 * a3 * 
g11 * g21 * g22 * g22 * beta * beta * cos(theta) * sin(theta) + 0.4e1 * a3 * g11 
* g21 * g12 * g22 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a6 * a6 - 
0.4e1 * a3 * g11 * g11 * g21 * g21 * alpha * alpha * pow(cos(theta), 0.2e1) - 
0.8e1 * a3 * g12 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a3 * 
g11 * g21 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a3 * g12 * g22 * 
alpha * alpha * sin(theta) * cos(theta) - 0.8e1 * a3 * g12 * pow(g22, 0.3e1) * 
alpha * alpha * sin(theta) * cos(theta) - 0.8e1 * a3 * g21 * g21 * g22 * g22 * 
alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a3 * g22 * g22 * alpha * alpha 
* sin(theta) * cos(theta) - 0.8e1 * a3 * g21 * g21 * g12 * g22 * alpha * alpha * 
sin(theta) * cos(theta) - 0.8e1 * a3 * g11 * pow(g21, 0.3e1) * alpha * alpha * 
sin(theta) * cos(theta) - 0.4e1 * a3 * alpha * alpha * pow(sin(theta), 0.2e1) - 
0.8e1 * a3 * g11 * g21 * g12 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) - 
0.4e1 * a3 * g12 * g12 * g22 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) - 
0.8e1 * a3 * g21 * g21 * alpha * alpha * sin(theta) * cos(theta) - 0.8e1 * a3 * 
g11 * g21 * alpha * alpha * sin(theta) * cos(theta) - 0.4e1 * a3 * alpha * alpha 
* pow(cos(theta), 0.2e1) - 0.8e1 * a3 * g11 * g21 * g22 * g22 * alpha * alpha * 
sin(theta) * cos(theta) + 0.4e1 * a3 * g21 * g21 * beta * beta * pow(cos(theta), 
0.2e1) + 0.4e1 * a3 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * 
a3 * pow(g22, 0.4e1) * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a3 * beta 
* beta * cos(theta) * sin(theta) - 0.4e1 * a3 * pow(g21, 0.4e1) * alpha * alpha 
* pow(sin(theta), 0.2e1) - 0.8e1 * a3 * g21 * g21 * alpha * alpha * 
pow(sin(theta), 0.2e1) - 0.8e1 * a3 * alpha * alpha * sin(theta) * cos(theta) - 
0.4e1 * a3 * pow(g22, 0.4e1) * alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * 
a3 * g22 * g22 * alpha * alpha * pow(sin(theta), 0.2e1);

double F0 = -0.2e1 * beta * cos(theta) * (-a6 * a6 * h22 * g22 - a6 * a6 * h21 * g21 - 
a3 * g11 * h21 * g22 * g22 * beta * beta * cos(theta) * sin(theta) + a3 * g11 * 
h21 * g12 * g22 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 
* g21 * g21 * beta * beta * cos(theta) * sin(theta) + a3 * pow(g22, 0.3e1) * h22 
* beta * beta * pow(cos(theta), 0.2e1) + a3 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) - a3 * g21 * h21 * g12 * g22 * beta * beta * cos(theta) * 
sin(theta) + a3 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g12 * 
g12 * h22 * g22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g21 * h21 * g22 * 
g22 * beta * beta * pow(cos(theta), 0.2e1) - a3 * g21 * g21 * g12 * h22 * beta * 
beta * cos(theta) * sin(theta) - a3 * g11 * h21 * beta * beta * cos(theta) * 
sin(theta) + a3 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - a3 * g11 * 
g21 * g22 * h22 * beta * beta * cos(theta) * sin(theta) - a3 * g22 * h22 * beta 
* beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g12 * h22 * g22 * g22 * beta * 
beta * cos(theta) * sin(theta) + a3 * g21 * h21 * beta * beta * pow(cos(theta), 
0.2e1) + a3 * g21 * g21 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - a3 
* g21 * h21 * beta * beta * cos(theta) * sin(theta) - a3 * g12 * h22 * beta * 
beta * cos(theta) * sin(theta) + a3 * g11 * g11 * h21 * g21 * beta * beta * 
pow(sin(theta), 0.2e1) + a3 * pow(g21, 0.3e1) * h21 * beta * beta * 
pow(cos(theta), 0.2e1) + a3 * g11 * g21 * g12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1));

double F1 = -0.2e1 * beta * cos(theta) * (-0.4e1 * a3 * g21 * h21 * alpha * sin(theta) 
* beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sin(theta) * 
beta * cos(theta) - 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha * sin(theta) * 
beta * cos(theta) + 0.2e1 * a3 * g21 * h21 * alpha * pow(sin(theta), 0.2e1) * 
beta - 0.2e1 * a3 * g21 * h21 * beta * pow(cos(theta), 0.2e1) * alpha + 0.2e1 * 
a3 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * g12 
* h22 * g22 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * h22 * 
g22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * h22 * g22 
* g22 * beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g21 * h21 * g22 * 
g22 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g21 * g21 * g12 * 
h22 * alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g21 * g21 * g12 * h22 
* beta * pow(cos(theta), 0.2e1) * alpha - 0.2e1 * a3 * g11 * h21 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * h21 * alpha * sin(theta) * 
beta * cos(theta) - 0.4e1 * a3 * g22 * h22 * alpha * sin(theta) * beta * 
cos(theta) + 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * pow(sin(theta), 0.2e1) 
* beta - 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * pow(cos(theta), 0.2e1) * 
alpha + 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sin(theta) * beta * 
cos(theta) + 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * pow(sin(theta), 0.2e1) 
* beta - 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * pow(cos(theta), 0.2e1) * 
alpha + 0.4e1 * a3 * g12 * h22 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 
* a3 * g11 * h21 * g12 * g22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * 
a3 * g11 * h21 * g21 * g21 * beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 
* pow(g22, 0.3e1) * h22 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a3 * 
g11 * h21 * g22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g11 
* h21 * alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g12 * h22 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * 
pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * 
pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * 
sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g22 * h22 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g22 * h22 * beta * pow(cos(theta), 
0.2e1) * alpha) + 0.4e1 * alpha * sin(theta) * (-a6 * a6 * h22 * g22 - a6 * a6 * 
h21 * g21 - a3 * g11 * h21 * g22 * g22 * beta * beta * cos(theta) * sin(theta) + 
a3 * g11 * h21 * g12 * g22 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * 
g11 * h21 * g21 * g21 * beta * beta * cos(theta) * sin(theta) + a3 * pow(g22, 
0.3e1) * h22 * beta * beta * pow(cos(theta), 0.2e1) + a3 * g22 * h22 * beta * 
beta * pow(cos(theta), 0.2e1) - a3 * g21 * h21 * g12 * g22 * beta * beta * 
cos(theta) * sin(theta) + a3 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) 
+ a3 * g12 * g12 * h22 * g22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g21 * 
h21 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) - a3 * g21 * g21 * g12 * 
h22 * beta * beta * cos(theta) * sin(theta) - a3 * g11 * h21 * beta * beta * 
cos(theta) * sin(theta) + a3 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) 
- a3 * g11 * g21 * g22 * h22 * beta * beta * cos(theta) * sin(theta) - a3 * g22 
* h22 * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g12 * h22 * g22 * 
g22 * beta * beta * cos(theta) * sin(theta) + a3 * g21 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) + a3 * g21 * g21 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) - a3 * g21 * h21 * beta * beta * cos(theta) * sin(theta) 
- a3 * g12 * h22 * beta * beta * cos(theta) * sin(theta) + a3 * g11 * g11 * h21 
* g21 * beta * beta * pow(sin(theta), 0.2e1) + a3 * pow(g21, 0.3e1) * h21 * beta 
* beta * pow(cos(theta), 0.2e1) + a3 * g11 * g21 * g12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1));

double F2 = -0.2e1 * beta * cos(theta) * (- (2 * a6 * a6 * h22 * g22) -  (2 * a6 * a6 * 
h21 * g21) + 0.2e1 * a3 * g11 *  h21 *  (g22 * g22) * beta * beta * cos(theta) * 
sin(theta) - 0.2e1 * a3 * g11 *  h21 * g12 *  g22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 *  h21 *  (g21 * g21) * beta * beta * 
cos(theta) * sin(theta) - 0.2e1 * a3 *   pow( g22,  3) *  h22 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 *  g22 *  h22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a3 *  g21 *  h21 * g12 *  g22 * beta * beta * 
cos(theta) * sin(theta) - 0.2e1 * a3 * g12 *  h22 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g12 * g12 *  h22 *  g22 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a3 *  g21 *  h21 *  (g22 * g22) * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a3 *  (g21 * g21) * g12 *  h22 * beta * beta * 
cos(theta) * sin(theta) + 0.2e1 * a3 * g11 *  h21 * beta * beta * cos(theta) * 
sin(theta) - 0.2e1 * a3 * g11 *  h21 * beta * beta * pow(sin(theta), 0.2e1) + 
0.2e1 * a3 * g11 *  g21 *  g22 *  h22 * beta * beta * cos(theta) * sin(theta) + 
0.2e1 * a3 *  g22 *  h22 * beta * beta * cos(theta) * sin(theta) + 0.4e1 * a3 * 
g12 *  h22 *  (g22 * g22) * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * 
g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 *  (g21 * g21) * 
g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a3 *  g21 *  h21 * 
beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g12 *  h22 * beta * beta * 
cos(theta) * sin(theta) - 0.2e1 * a3 * g11 * g11 *  h21 *  g21 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a3 *   pow( g21,  3) *  h21 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 *  g21 * g12 *  h22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 *  g21 *  h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 *  (g21 * g21) *  g22 *  h22 * alpha * alpha 
* pow(sin(theta), 0.2e1) + 0.4e1 * a3 *   pow( g21,  3) *  h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 *  g21 *  h21 * alpha * alpha * sin(theta) * 
cos(theta) + 0.4e1 * a3 * g12 *  h22 * alpha * alpha * sin(theta) * cos(theta) + 
0.4e1 * a3 * g12 * g12 *  h22 *  g22 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.4e1 * a3 *  (g21 * g21) * g12 *  h22 * alpha * alpha * sin(theta) * cos(theta) 
+ 0.4e1 * a3 * g11 *  h21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 
* g11 *  h21 * alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 *  g21 * 
h21 * g12 *  g22 * alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 * g11 * 
g21 * g12 *  h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * 
h21 * g12 *  g22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * 
h21 *  (g22 * g22) * alpha * alpha * sin(theta) * cos(theta) + 0.8e1 * a3 * g11 
*  h21 *  (g21 * g21) * alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 * 
g22 *  h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 *  g21 *  h21 * 
(g22 * g22) * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 *   pow( g22, 
3) *  h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 * g11 * 
h21 *  g21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 *  h22 * 
alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 *  g21 *  g22 *  h22 * 
alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 *  g22 *  h22 * alpha * 
alpha * sin(theta) * cos(theta) + 0.8e1 * a3 * g12 *  h22 *  (g22 * g22) * alpha 
* alpha * sin(theta) * cos(theta)) + 0.4e1 * alpha * sin(theta) * (-0.4e1 * a3 * 
g21 *  h21 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 *  (g21 * g21) 
*  g22 *  h22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 *   pow( 
g21,  3) *  h21 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 *  g21 * 
h21 * alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 *  g21 *  h21 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g12 *  h22 * alpha * 
pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * g12 *  h22 *  g22 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 *  h22 *  (g22 * g22) * alpha 
* pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 *  h22 *  (g22 * g22) * beta 
* pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 *  g21 *  h21 *  (g22 * g22) * 
alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 *  (g21 * g21) * g12 *  h22 
* alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 *  (g21 * g21) * g12 *  h22 
* beta * pow(cos(theta), 0.2e1) * alpha - 0.2e1 * a3 * g11 *  h21 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 *  h21 * alpha * sin(theta) * 
beta * cos(theta) - 0.4e1 * a3 *  g22 *  h22 * alpha * sin(theta) * beta * 
cos(theta) + 0.2e1 * a3 *  g21 *  h21 * g12 *  g22 * alpha * pow(sin(theta), 
0.2e1) * beta - 0.2e1 * a3 *  g21 *  h21 * g12 *  g22 * beta * pow(cos(theta), 
0.2e1) * alpha + 0.4e1 * a3 * g11 *  g21 * g12 *  h22 * alpha * sin(theta) * 
beta * cos(theta) + 0.2e1 * a3 * g11 *  g21 *  g22 *  h22 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g11 *  g21 *  g22 *  h22 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 *  h22 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a3 * g11 *  h21 * g12 *  g22 * alpha * sin(theta) * 
beta * cos(theta) - 0.4e1 * a3 * g11 *  h21 *  (g21 * g21) * beta * 
pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 *   pow( g22,  3) *  h22 * alpha * 
sin(theta) * beta * cos(theta) - 0.2e1 * a3 * g11 *  h21 *  (g22 * g22) * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g11 *  h21 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g12 *  h22 * beta * pow(cos(theta), 
0.2e1) * alpha + 0.4e1 * a3 * g11 *  h21 *  (g21 * g21) * alpha * 
pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g11 *  h21 *  (g22 * g22) * alpha * 
pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * g11 *  h21 *  g21 * alpha * 
sin(theta) * beta * cos(theta) + 0.2e1 * a3 *  g22 *  h22 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 *  g22 *  h22 * beta * 
pow(cos(theta), 0.2e1) * alpha) + 0.2e1 * beta * cos(theta) * (- (a6 * a6 * h22 
* g22) -  (a6 * a6 * h21 * g21) - a3 * g11 *  h21 *  (g22 * g22) * beta * beta * 
cos(theta) * sin(theta) + a3 * g11 *  h21 * g12 *  g22 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 *  h21 *  (g21 * g21) * beta * beta * 
cos(theta) * sin(theta) + a3 *   pow( g22,  3) *  h22 * beta * beta * 
pow(cos(theta), 0.2e1) + a3 *  g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1) 
- a3 *  g21 *  h21 * g12 *  g22 * beta * beta * cos(theta) * sin(theta) + a3 * 
g12 *  h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g12 * g12 *  h22 *  g22 
* beta * beta * pow(sin(theta), 0.2e1) + a3 *  g21 *  h21 *  (g22 * g22) * beta 
* beta * pow(cos(theta), 0.2e1) - a3 *  (g21 * g21) * g12 *  h22 * beta * beta * 
cos(theta) * sin(theta) - a3 * g11 *  h21 * beta * beta * cos(theta) * 
sin(theta) + a3 * g11 *  h21 * beta * beta * pow(sin(theta), 0.2e1) - a3 * g11 * 
g21 *  g22 *  h22 * beta * beta * cos(theta) * sin(theta) - a3 *  g22 *  h22 * 
beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g12 *  h22 *  (g22 * g22) * 
beta * beta * cos(theta) * sin(theta) + a3 *  g21 *  h21 * beta * beta * 
pow(cos(theta), 0.2e1) + a3 *  (g21 * g21) *  g22 *  h22 * beta * beta * 
pow(cos(theta), 0.2e1) - a3 *  g21 *  h21 * beta * beta * cos(theta) * 
sin(theta) - a3 * g12 *  h22 * beta * beta * cos(theta) * sin(theta) + a3 * g11 
* g11 *  h21 *  g21 * beta * beta * pow(sin(theta), 0.2e1) + a3 *   pow( g21, 
3) *  h21 * beta * beta * pow(cos(theta), 0.2e1) + a3 * g11 *  g21 * g12 *  h22 
* beta * beta * pow(sin(theta), 0.2e1));

double F3 = -0.2e1 * beta * cos(theta) * (0.4e1 * a3 * g21 * h21 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha * sin(theta) * 
beta * cos(theta) - 0.2e1 * a3 * g21 * h21 * alpha * pow(sin(theta), 0.2e1) * 
beta + 0.2e1 * a3 * g21 * h21 * beta * pow(cos(theta), 0.2e1) * alpha - 0.2e1 * 
a3 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * g12 
* h22 * g22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g12 * h22 * 
g22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g12 * h22 * g22 
* g22 * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g21 * h21 * g22 * 
g22 * alpha * sin(theta) * beta * cos(theta) - 0.2e1 * a3 * g21 * g21 * g12 * 
h22 * alpha * pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g21 * g21 * g12 * h22 
* beta * pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g11 * h21 * beta * 
pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * alpha * sin(theta) * 
beta * cos(theta) + 0.4e1 * a3 * g22 * h22 * alpha * sin(theta) * beta * 
cos(theta) - 0.2e1 * a3 * g21 * h21 * g12 * g22 * alpha * pow(sin(theta), 0.2e1) 
* beta + 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * pow(cos(theta), 0.2e1) * 
alpha - 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * sin(theta) * beta * 
cos(theta) - 0.2e1 * a3 * g11 * g21 * g22 * h22 * alpha * pow(sin(theta), 0.2e1) 
* beta + 0.2e1 * a3 * g11 * g21 * g22 * h22 * beta * pow(cos(theta), 0.2e1) * 
alpha - 0.4e1 * a3 * g12 * h22 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 
* a3 * g11 * h21 * g12 * g22 * alpha * sin(theta) * beta * cos(theta) + 0.4e1 * 
a3 * g11 * h21 * g21 * g21 * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 
* pow(g22, 0.3e1) * h22 * alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 * 
g11 * h21 * g22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha - 0.2e1 * a3 * g11 
* h21 * alpha * pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g12 * h22 * beta * 
pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * 
sin(theta) * beta * cos(theta) - 0.2e1 * a3 * g22 * h22 * alpha * 
pow(sin(theta), 0.2e1) * beta + 0.2e1 * a3 * g22 * h22 * beta * pow(cos(theta), 
0.2e1) * alpha) + 0.4e1 * alpha * sin(theta) * (-0.2e1 * a6 * a6 * h22 * g22 - 
0.2e1 * a6 * a6 * h21 * g21 + 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * beta * 
cos(theta) * sin(theta) - 0.2e1 * a3 * g11 * h21 * g12 * g22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * beta * 
cos(theta) * sin(theta) - 0.2e1 * a3 * pow(g22, 0.3e1) * h22 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g22 * h22 * beta * beta * pow(cos(theta), 
0.2e1) + 0.2e1 * a3 * g21 * h21 * g12 * g22 * beta * beta * cos(theta) * 
sin(theta) - 0.2e1 * a3 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) - 
0.2e1 * a3 * g12 * g12 * h22 * g22 * beta * beta * pow(sin(theta), 0.2e1) - 
0.2e1 * a3 * g21 * h21 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) + 
0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * beta * cos(theta) * sin(theta) + 
0.2e1 * a3 * g11 * h21 * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * 
g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a3 * g11 * g21 * g22 
* h22 * beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g22 * h22 * beta * 
beta * cos(theta) * sin(theta) + 0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * 
beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g21 * g21 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a3 * g21 * h21 * beta * beta * cos(theta) * 
sin(theta) + 0.2e1 * a3 * g12 * h22 * beta * beta * cos(theta) * sin(theta) - 
0.2e1 * a3 * g11 * g11 * h21 * g21 * beta * beta * pow(sin(theta), 0.2e1) - 
0.2e1 * a3 * pow(g21, 0.3e1) * h21 * beta * beta * pow(cos(theta), 0.2e1) - 
0.2e1 * a3 * g11 * g21 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) + 
0.4e1 * a3 * g21 * h21 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 * 
g21 * g21 * g22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 * 
pow(g21, 0.3e1) * h21 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a3 * 
g21 * h21 * alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 * g12 * h22 * 
alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 * g12 * g12 * h22 * g22 * 
alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g21 * g21 * g12 * h22 * 
alpha * alpha * sin(theta) * cos(theta) + 0.4e1 * a3 * g11 * h21 * alpha * alpha 
* pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * alpha * alpha * sin(theta) * 
cos(theta) + 0.4e1 * a3 * g21 * h21 * g12 * g22 * alpha * alpha * sin(theta) * 
cos(theta) + 0.4e1 * a3 * g11 * g21 * g12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * g22 * g22 * alpha * alpha * 
sin(theta) * cos(theta) + 0.8e1 * a3 * g11 * h21 * g21 * g21 * alpha * alpha * 
sin(theta) * cos(theta) + 0.4e1 * a3 * g22 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 * pow(g22, 0.3e1) * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * g21 * g22 * h22 * alpha * alpha * 
sin(theta) * cos(theta) + 0.4e1 * a3 * g22 * h22 * alpha * alpha * sin(theta) * 
cos(theta) + 0.8e1 * a3 * g12 * h22 * g22 * g22 * alpha * alpha * sin(theta) * 
cos(theta)) + 0.2e1 * beta * cos(theta) * (-0.4e1 * a3 * g21 * h21 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * g22 * h22 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a3 * pow(g21, 0.3e1) * h21 * alpha * 
sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g21 * h21 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g21 * h21 * beta * pow(cos(theta), 
0.2e1) * alpha + 0.2e1 * a3 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta 
+ 0.4e1 * a3 * g12 * g12 * h22 * g22 * alpha * sin(theta) * beta * cos(theta) + 
0.4e1 * a3 * g12 * h22 * g22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta - 
0.4e1 * a3 * g12 * h22 * g22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha - 
0.4e1 * a3 * g21 * h21 * g22 * g22 * alpha * sin(theta) * beta * cos(theta) + 
0.2e1 * a3 * g21 * g21 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta - 
0.2e1 * a3 * g21 * g21 * g12 * h22 * beta * pow(cos(theta), 0.2e1) * alpha - 
0.2e1 * a3 * g11 * h21 * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * 
g11 * h21 * alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g22 * h22 * 
alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g21 * h21 * g12 * g22 * 
alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g21 * h21 * g12 * g22 * 
beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * g21 * g12 * h22 * 
alpha * sin(theta) * beta * cos(theta) + 0.2e1 * a3 * g11 * g21 * g22 * h22 * 
alpha * pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g11 * g21 * g22 * h22 * 
beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 * h22 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * g12 * g22 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g11 * h21 * g21 * g21 * beta * 
pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * pow(g22, 0.3e1) * h22 * alpha * 
sin(theta) * beta * cos(theta) - 0.2e1 * a3 * g11 * h21 * g22 * g22 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.2e1 * a3 * g11 * h21 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.2e1 * a3 * g12 * h22 * beta * pow(cos(theta), 
0.2e1) * alpha + 0.4e1 * a3 * g11 * h21 * g21 * g21 * alpha * pow(sin(theta), 
0.2e1) * beta + 0.2e1 * a3 * g11 * h21 * g22 * g22 * alpha * pow(sin(theta), 
0.2e1) * beta + 0.4e1 * a3 * g11 * g11 * h21 * g21 * alpha * sin(theta) * beta * 
cos(theta) + 0.2e1 * a3 * g22 * h22 * alpha * pow(sin(theta), 0.2e1) * beta - 
0.2e1 * a3 * g22 * h22 * beta * pow(cos(theta), 0.2e1) * alpha);



double G0 = -beta * beta * pow(cos(theta), 0.2e1) * (-a6 * a6 * h21 * h21 - a6 * a6 * 
h22 * h22 + a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sin(theta), 0.2e1) + 
0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * 
g11 * g11 * h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 
* h21 * h21 * g21 * beta * beta * cos(theta) * sin(theta) + a3 * g21 * g21 * h21 
* h21 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * h22 * h22 * 
g22 * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * g12 * h22 
* beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g21 * h21 * g22 * h22 * 
beta * beta * pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 * h22 * beta * beta 
* pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * 
cos(theta) * sin(theta));

G1= -beta * beta * pow(cos(theta), 0.2e1) * (-0.4e1 * a3 * g21 * h21 * g12 * h22 
* beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * h21 * g21 * 
beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * g11 * h21 * h21 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * h21 * h21 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * g12 * h22 * h22 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * g22 * h22 * 
alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * h22 * h22 * g22 * 
beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 * h22 * h22 * g22 * 
alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g11 * h21 * g22 * h22 * 
beta * pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 * g21 * h21 * g22 * h22 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g22 * g22 * h22 * h22 * 
alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a3 * g11 * h21 * g12 * h22 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g21 * h21 * g12 * h22 * 
alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * h21 * h21 * g21 * 
alpha * pow(sin(theta), 0.2e1) * beta) + 0.4e1 * beta * cos(theta) * alpha * 
sin(theta) * (-a6 * a6 * h21 * h21 - a6 * a6 * h22 * h22 + a3 * g12 * g12 * h22 
* h22 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a3 * g11 * h21 * g12 * 
h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g11 * g11 * h21 * h21 * beta * 
beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta 
* cos(theta) * sin(theta) + a3 * g21 * g21 * h21 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * 
cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * 
cos(theta) * sin(theta) + 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * 
cos(theta) * sin(theta));

double G2 = -beta * beta * pow(cos(theta), 0.2e1) * (-0.4e1 * a3 * g11 * h21 * g12 * 
h22 * beta * beta * pow(sin(theta), 0.2e1) + 0.8e1 * a3 * g11 * h21 * g12 * h22 
* alpha * alpha * pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g22 * g22 * h22 * h22 * 
beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta 
* beta * cos(theta) * sin(theta) - 0.2e1 * a6 * a6 * h21 * h21 + 0.8e1 * a3 * 
g21 * h21 * g12 * h22 * alpha * alpha * sin(theta) * cos(theta) - 0.2e1 * a3 * 
g11 * g11 * h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a6 * a6 * 
h22 * h22 + 0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sin(theta) * 
cos(theta) + 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cos(theta) * 
sin(theta) + 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) - 0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * 
sin(theta) * cos(theta) - 0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * 
cos(theta) * sin(theta) + 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * 
cos(theta) * sin(theta) + 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * 
sin(theta) * cos(theta)) + 0.4e1 * beta * cos(theta) * alpha * sin(theta) * 
(-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cos(theta), 0.2e1) * alpha - 
0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cos(theta), 0.2e1) * alpha + 
0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sin(theta) * beta * cos(theta) - 
0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sin(theta) * beta * cos(theta) + 
0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sin(theta) * beta * cos(theta) + 
0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sin(theta), 0.2e1) * beta - 
0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cos(theta), 0.2e1) * alpha + 
0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sin(theta), 0.2e1) * beta - 
0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cos(theta), 0.2e1) * alpha - 
0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sin(theta) * beta * cos(theta) - 
0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sin(theta) * beta * cos(theta) + 
0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sin(theta) * beta * cos(theta) + 
0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * beta + 
0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sin(theta), 0.2e1) * beta) - 
(-0.2e1 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * alpha * alpha * 
pow(sin(theta), 0.2e1)) * (-a6 * a6 * h21 * h21 - a6 * a6 * h22 * h22 + a3 * g12 
* g12 * h22 * h22 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a3 * g11 * 
h21 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) + a3 * g11 * g11 * h21 * 
h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * h21 * g21 
* beta * beta * cos(theta) * sin(theta) + a3 * g21 * g21 * h21 * h21 * beta * 
beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta 
* cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * 
cos(theta) * sin(theta) + 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * 
cos(theta) * sin(theta));

double G3 = -beta * beta * pow(cos(theta), 0.2e1) * (0.4e1 * a3 * g22 * g22 * h22 * h22 
* alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g21 * g21 * h21 * h21 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * h21 * g21 * 
beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g12 * h22 * h22 * g22 * 
alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * h21 * g22 * h22 * 
beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * g11 * h21 * h21 * 
alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a3 * g21 * h21 * g22 * h22 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g12 * g12 * h22 * h22 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * h22 * h22 * g22 * 
beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * g22 * h22 * 
alpha * pow(sin(theta), 0.2e1) * beta - 0.8e1 * a3 * g11 * h21 * g12 * h22 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g11 * h21 * h21 * g21 * 
alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g21 * h21 * g12 * h22 * 
alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g21 * h21 * g12 * h22 * 
beta * pow(cos(theta), 0.2e1) * alpha) + 0.4e1 * beta * cos(theta) * alpha * 
sin(theta) * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * 
cos(theta) * sin(theta) - 0.2e1 * a6 * a6 * h21 * h21 + 0.8e1 * a3 * g21 * h21 * 
g12 * h22 * alpha * alpha * sin(theta) * cos(theta) - 0.2e1 * a3 * g11 * g11 * 
h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a6 * a6 * h22 * h22 + 
0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sin(theta) * cos(theta) + 
0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cos(theta) * sin(theta) + 
0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * pow(sin(theta), 0.2e1) - 
0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 
0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sin(theta), 0.2e1) + 
0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sin(theta) * cos(theta) - 
0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + 
0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cos(theta) * sin(theta) + 
0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cos(theta) * sin(theta) + 
0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sin(theta) * cos(theta)) - 
(-0.2e1 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * alpha * alpha * 
pow(sin(theta), 0.2e1)) * (-0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * 
pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * 
pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * 
pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * 
pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * 
sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * 
sin(theta) * beta * cos(theta) + 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * 
sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * 
pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * 
pow(sin(theta), 0.2e1) * beta) - 0.4e1 * beta * cos(theta) * alpha * sin(theta) 
* (-a6 * a6 * h21 * h21 - a6 * a6 * h22 * h22 + a3 * g12 * g12 * h22 * h22 * 
beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta 
* beta * pow(sin(theta), 0.2e1) + a3 * g11 * g11 * h21 * h21 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * 
cos(theta) * sin(theta) + a3 * g21 * g21 * h21 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * 
cos(theta) * sin(theta) - 0.2e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * 
cos(theta) * sin(theta) + 0.2e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * 
cos(theta) * sin(theta));

double G4 = -0.2e1 * beta * beta * pow(cos(theta), 0.2e1) * (-a6 * a6 * h21 * h21 - a6 
* a6 * h22 * h22 + a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * pow(sin(theta), 
0.2e1) + a3 * g11 * g11 * h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 
0.2e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cos(theta) * sin(theta) + a3 
* g21 * g21 * h21 * h21 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 * 
g12 * h22 * h22 * g22 * beta * beta * cos(theta) * sin(theta) - 0.2e1 * a3 * g21 
* h21 * g12 * h22 * beta * beta * cos(theta) * sin(theta) + 0.2e1 * a3 * g21 * 
h21 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) + a3 * g22 * g22 * h22 * 
h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g11 * h21 * g22 * h22 
* beta * beta * cos(theta) * sin(theta)) + 0.4e1 * beta * cos(theta) * alpha * 
sin(theta) * (0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * sin(theta) * beta * 
cos(theta) + 0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * sin(theta) * beta * 
cos(theta) + 0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * pow(cos(theta), 0.2e1) 
* alpha - 0.4e1 * a3 * g12 * h22 * h22 * g22 * alpha * pow(sin(theta), 0.2e1) * 
beta + 0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * pow(cos(theta), 0.2e1) * 
alpha - 0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * sin(theta) * beta * 
cos(theta) + 0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * sin(theta) * beta * 
cos(theta) - 0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * sin(theta) * beta * 
cos(theta) + 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * pow(cos(theta), 0.2e1) 
* alpha - 0.4e1 * a3 * g11 * h21 * g22 * h22 * alpha * pow(sin(theta), 0.2e1) * 
beta - 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * sin(theta) * beta * 
cos(theta) - 0.4e1 * a3 * g11 * h21 * h21 * g21 * alpha * pow(sin(theta), 0.2e1) 
* beta - 0.4e1 * a3 * g21 * h21 * g12 * h22 * alpha * pow(sin(theta), 0.2e1) * 
beta + 0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * pow(cos(theta), 0.2e1) * 
alpha) - (-0.2e1 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * alpha * alpha 
* pow(sin(theta), 0.2e1)) * (-0.4e1 * a3 * g11 * h21 * g12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.8e1 * a3 * g11 * h21 * g12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) - 0.2e1 * a3 * g22 * g22 * h22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.4e1 * a3 * g12 * h22 * h22 * g22 * beta * beta * 
cos(theta) * sin(theta) - 0.2e1 * a6 * a6 * h21 * h21 + 0.8e1 * a3 * g21 * h21 * 
g12 * h22 * alpha * alpha * sin(theta) * cos(theta) - 0.2e1 * a3 * g11 * g11 * 
h21 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a6 * a6 * h22 * h22 + 
0.8e1 * a3 * g11 * h21 * h21 * g21 * alpha * alpha * sin(theta) * cos(theta) + 
0.4e1 * a3 * g11 * h21 * h21 * g21 * beta * beta * cos(theta) * sin(theta) + 
0.4e1 * a3 * g22 * g22 * h22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a3 * g12 * g12 * h22 * h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.4e1 * a3 * g21 * g21 * h21 * h21 * alpha * alpha * pow(sin(theta), 0.2e1) - 
0.4e1 * a3 * g21 * h21 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 
0.2e1 * a3 * g12 * g12 * h22 * h22 * beta * beta * pow(sin(theta), 0.2e1) + 
0.8e1 * a3 * g12 * h22 * h22 * g22 * alpha * alpha * sin(theta) * cos(theta) - 
0.2e1 * a3 * g21 * g21 * h21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + 
0.4e1 * a3 * g11 * h21 * g22 * h22 * beta * beta * cos(theta) * sin(theta) + 
0.8e1 * a3 * g21 * h21 * g22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a3 * g21 * h21 * g12 * h22 * beta * beta * cos(theta) * sin(theta) + 
0.4e1 * a3 * g11 * g11 * h21 * h21 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.8e1 * a3 * g11 * h21 * g22 * h22 * alpha * alpha * sin(theta) * cos(theta)) - 
0.4e1 * beta * cos(theta) * alpha * sin(theta) * (-0.4e1 * a3 * g21 * h21 * g12 
* h22 * beta * pow(cos(theta), 0.2e1) * alpha - 0.4e1 * a3 * g11 * h21 * h21 * 
g21 * beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g11 * g11 * h21 * h21 
* alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g21 * g21 * h21 * h21 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g12 * g12 * h22 * h22 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g11 * h21 * g22 * h22 * 
alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g12 * h22 * h22 * g22 * 
beta * pow(cos(theta), 0.2e1) * alpha + 0.4e1 * a3 * g12 * h22 * h22 * g22 * 
alpha * pow(sin(theta), 0.2e1) * beta - 0.4e1 * a3 * g11 * h21 * g22 * h22 * 
beta * pow(cos(theta), 0.2e1) * alpha - 0.8e1 * a3 * g21 * h21 * g22 * h22 * 
alpha * sin(theta) * beta * cos(theta) - 0.4e1 * a3 * g22 * g22 * h22 * h22 * 
alpha * sin(theta) * beta * cos(theta) + 0.8e1 * a3 * g11 * h21 * g12 * h22 * 
alpha * sin(theta) * beta * cos(theta) + 0.4e1 * a3 * g21 * h21 * g12 * h22 * 
alpha * pow(sin(theta), 0.2e1) * beta + 0.4e1 * a3 * g11 * h21 * h21 * g21 * 
alpha * pow(sin(theta), 0.2e1) * beta);

double H0 = -beta * beta * sin(theta) * cos(theta) * (a5 * g11 * g11 * h11 * h21 * beta 
* beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 * a6 * h12 * h22 - a5 
* g22 * h12 * g11 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g22 * h12 
* g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g22 * g22 * h12 * h22 
* beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * h11 * g21 * h21 * 
beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * h22 * beta * beta 
* pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 * beta * beta * sin(theta) 
* cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta * pow(cos(theta), 0.2e1) 
- a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta) * cos(theta) + a5 * g12 
* h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g21 * h11 * g22 
* h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g12 * h12 * g22 * 
h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12 * g21 * h21 * beta 
* beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1));

double H1 = -beta * beta * sin(theta) * cos(theta) * (0.4e1 * a5 * g11 * g11 * h11 * 
h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * 
h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 * 
h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g11 * h21 
* alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * h12 * g21 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * h11 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * h11 * g22 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g22 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * beta * pow(sin(theta), 0.2e1) 
* alpha + 0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta) * (a5 * g11 * g11 * h11 
* h21 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 * a6 * 
h12 * h22 - a5 * g22 * h12 * g11 * h21 * beta * beta * sin(theta) * cos(theta) + 
a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g22 * 
g22 * h12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * h11 
* g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * h22 
* beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 * beta * 
beta * sin(theta) * cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) - a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta) * 
cos(theta) + a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) + 
a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * 
g12 * h12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12 * 
g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22 * 
beta * beta * pow(sin(theta), 0.2e1));

double H2 = -beta * beta * sin(theta) * cos(theta) * (-0.2e1 * a5 * g11 * g11 * h11 * 
h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a4 * a6 * h11 * h21 - 0.2e1 
* a4 * a6 * h12 * h22 + 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha * 
cos(theta) * sin(theta)) - (-0.2e1 * beta * pow(sin(theta), 0.2e1) * alpha + 
0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta) * (0.4e1 * a5 * g11 * g11 * h11 * 
h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * 
h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 * 
h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g11 * h21 
* alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * h12 * g21 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * h11 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * h11 * g22 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g22 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * beta * beta * sin(theta) * 
cos(theta) - 0.4e1 * alpha * alpha * cos(theta) * sin(theta)) * (a5 * g11 * g11 
* h11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 * 
a6 * h12 * h22 - a5 * g22 * h12 * g11 * h21 * beta * beta * sin(theta) * 
cos(theta) + a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + 
a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * 
g11 * h11 * g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * 
g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 * 
beta * beta * sin(theta) * cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta 
* pow(cos(theta), 0.2e1) - a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta) 
* cos(theta) + a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) 
+ a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 
* g12 * h12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12 
* g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22 
* beta * beta * pow(sin(theta), 0.2e1));

double H3 = -beta * beta * sin(theta) * cos(theta) * (0.4e1 * a5 * g21 * g21 * h11 * 
h21 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g11 * h11 * g12 * 
h22 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 * 
h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * h12 * g11 * h21 
* alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h11 * g12 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g21 * h11 * g22 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g22 * h12 * g21 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * g22 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g11 * g11 * h11 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 * g12 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * g12 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * beta * pow(sin(theta), 0.2e1) 
* alpha + 0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta) * (-0.2e1 * a5 * g11 * 
g11 * h11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a4 * a6 * h11 * 
h21 - 0.2e1 * a4 * a6 * h12 * h22 + 0.2e1 * a5 * g22 * h12 * g11 * h21 * beta * 
beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * h21 * beta * 
beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta * beta 
* pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha * 
cos(theta) * sin(theta)) - (-0.2e1 * beta * beta * sin(theta) * cos(theta) - 
0.4e1 * alpha * alpha * cos(theta) * sin(theta)) * (0.4e1 * a5 * g11 * g11 * h11 
* h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * 
h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 * 
h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g11 * h21 
* alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * h12 * g21 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * h11 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * h11 * g22 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g22 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.2e1 * beta * pow(sin(theta), 0.2e1) * alpha) * (a5 * g11 * g11 
* h11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 * 
a6 * h12 * h22 - a5 * g22 * h12 * g11 * h21 * beta * beta * sin(theta) * 
cos(theta) + a5 * g22 * h12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + 
a5 * g22 * g22 * h12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * 
g11 * h11 * g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * 
g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 * 
beta * beta * sin(theta) * cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta 
* pow(cos(theta), 0.2e1) - a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta) 
* cos(theta) + a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) 
+ a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 
* g12 * h12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12 
* g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22 
* beta * beta * pow(sin(theta), 0.2e1));

double H4 = -0.2e1 * beta * beta * sin(theta) * cos(theta) * (a5 * g11 * g11 * h11 * 
h21 * beta * beta * pow(sin(theta), 0.2e1) - a4 * a6 * h11 * h21 - a4 * a6 * h12 
* h22 - a5 * g22 * h12 * g11 * h21 * beta * beta * sin(theta) * cos(theta) + a5 
* g22 * h12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g22 * g22 
* h12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * h11 * 
g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * h22 * 
beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * h22 * beta * beta 
* sin(theta) * cos(theta) + a5 * g21 * g21 * h11 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) - a5 * g21 * h11 * g12 * h22 * beta * beta * sin(theta) * 
cos(theta) + a5 * g12 * h12 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) + 
a5 * g21 * h11 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * 
g12 * h12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h12 * 
g21 * h21 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * h22 * 
beta * beta * pow(sin(theta), 0.2e1)) - (-0.2e1 * beta * pow(sin(theta), 0.2e1) 
* alpha + 0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta) * (0.4e1 * a5 * g21 * 
g21 * h11 * h21 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g11 * 
h11 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * 
h12 * g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * h12 
* g11 * h21 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * 
g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 
* h21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g22 * 
h22 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g11 * h11 * g22 * h22 
* alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h11 * g12 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g21 * h11 * g22 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g22 * h12 * g21 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * g22 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g11 * g11 * h11 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 * g12 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * g12 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta) - (-0.2e1 * beta * beta * sin(theta) * 
cos(theta) - 0.4e1 * alpha * alpha * cos(theta) * sin(theta)) * (-0.2e1 * a5 * 
g11 * g11 * h11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a4 * a6 * 
h11 * h21 - 0.2e1 * a4 * a6 * h12 * h22 + 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * h21 * 
beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 * h12 * h22 * beta 
* beta * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g21 * h21 * beta * 
beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * h22 * beta * 
beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * beta * beta 
* sin(theta) * cos(theta) - 0.2e1 * a5 * g21 * g21 * h11 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * h22 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * h21 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a5 * g12 * h12 * g21 * h21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * g22 * h22 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * g21 * h11 * h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.8e1 * a5 * g11 * h11 * g21 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.8e1 * a5 * g12 * h12 * g22 * h22 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * g22 * h12 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * g11 * h11 * h21 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g11 * h21 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * h11 * g12 * h22 * alpha * alpha * 
cos(theta) * sin(theta)) - (-0.2e1 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.2e1 * beta * pow(sin(theta), 0.2e1) * alpha) * (0.4e1 * a5 * g11 * g11 * h11 * 
h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * 
h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 * 
h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g12 * h12 * g11 * h21 
* alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g22 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * h12 * g21 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * h11 * h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * h11 * g22 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * g21 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g22 * h12 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * h12 * g22 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta);

double J0 = -beta * cos(theta) * (-a4 * a6 * g12 * h22 - a4 * a6 * g11 * h21 + a5 * g12 
* h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * g12 * g21 * h21 * beta 
* beta * pow(cos(theta), 0.2e1) + a5 * g12 * g12 * g11 * h21 * beta * beta * 
pow(sin(theta), 0.2e1) - a5 * g12 * g12 * g21 * h21 * beta * beta * sin(theta) * 
cos(theta) - a5 * g22 * g12 * g11 * h21 * beta * beta * sin(theta) * cos(theta) 
- a5 * g22 * h22 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h21 * beta 
* beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g11 * g11 * g21 * h21 * beta * 
beta * sin(theta) * cos(theta) - a5 * g12 * h22 * beta * beta * sin(theta) * 
cos(theta) - a5 * g11 * h21 * beta * beta * sin(theta) * cos(theta) - a5 * g21 * 
g11 * g12 * h22 * beta * beta * sin(theta) * cos(theta) + a5 * pow(g12, 0.3e1) * 
h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * g22 * g12 * h22 * beta * 
beta * pow(cos(theta), 0.2e1) - a5 * g21 * h21 * beta * beta * sin(theta) * 
cos(theta) + a5 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 * 
pow(g11, 0.3e1) * h21 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g11 * g11 * 
g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * g11 * g22 * h22 * 
beta * beta * sin(theta) * cos(theta) + a5 * g21 * g21 * g11 * h21 * beta * beta 
* pow(cos(theta), 0.2e1) + a5 * g21 * g11 * g22 * h22 * beta * beta * 
pow(cos(theta), 0.2e1) + a5 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) - 
0.2e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sin(theta) * cos(theta));

double J1 = -beta * cos(theta) * (0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * 
cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * 
cos(theta) * beta * sin(theta) + 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha * 
cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g12 * h22 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a5 * g21 * h21 * alpha * cos(theta) * beta * sin(theta) 
+ 0.4e1 * a5 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * 
pow(g11, 0.3e1) * h21 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * 
g12 * h22 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * g12 * g21 
* h21 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * g12 * g11 * 
h21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 * g12 * g11 * h21 
* beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * h21 * alpha * 
cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h21 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g22 * h22 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cos(theta) * beta * 
sin(theta) - 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) 
* beta + 0.2e1 * a5 * g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 
* a5 * g22 * g22 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * 
a5 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h21 
* beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g11 * g22 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g22 * h22 * alpha * 
cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g11 * h21 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cos(theta) * beta 
* sin(theta)) + 0.2e1 * alpha * sin(theta) * (-a4 * a6 * g12 * h22 - a4 * a6 * 
g11 * h21 + a5 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * 
g12 * g21 * h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g12 * g12 * g11 * 
h21 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g12 * g12 * g21 * h21 * beta * 
beta * sin(theta) * cos(theta) - a5 * g22 * g12 * g11 * h21 * beta * beta * 
sin(theta) * cos(theta) - a5 * g22 * h22 * beta * beta * sin(theta) * cos(theta) 
+ a5 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g11 * g11 
* g21 * h21 * beta * beta * sin(theta) * cos(theta) - a5 * g12 * h22 * beta * 
beta * sin(theta) * cos(theta) - a5 * g11 * h21 * beta * beta * sin(theta) * 
cos(theta) - a5 * g21 * g11 * g12 * h22 * beta * beta * sin(theta) * cos(theta) 
+ a5 * pow(g12, 0.3e1) * h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * 
g22 * g12 * h22 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g21 * h21 * beta * 
beta * sin(theta) * cos(theta) + a5 * g21 * h21 * beta * beta * pow(cos(theta), 
0.2e1) + a5 * pow(g11, 0.3e1) * h21 * beta * beta * pow(sin(theta), 0.2e1) + a5 
* g11 * g11 * g12 * h22 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * g11 
* g22 * h22 * beta * beta * sin(theta) * cos(theta) + a5 * g21 * g21 * g11 * h21 
* beta * beta * pow(cos(theta), 0.2e1) + a5 * g21 * g11 * g22 * h22 * beta * 
beta * pow(cos(theta), 0.2e1) + a5 * g22 * h22 * beta * beta * pow(cos(theta), 
0.2e1) - 0.2e1 * a5 * g12 * g12 * g22 * h22 * beta * beta * sin(theta) * 
cos(theta));

double J2 = -beta * cos(theta) * (- (2 * a4 * a6 * g12 * h22) -  (2 * a4 * a6 * g11 * 
h21) - 0.2e1 * a5 *  g12 *  h22 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * 
a5 * g22 *  g12 * g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 
*  (g12 * g12) *  g11 *  h21 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 
*  (g12 * g12) * g21 *  h21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 
* g22 *  g12 *  g11 *  h21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 
* g22 *  h22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 *  g11 *  h21 
* beta * beta * pow(sin(theta), 0.2e1) + 0.4e1 * a5 *  (g11 * g11) * g21 *  h21 
* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 *  g12 *  h22 * beta * beta 
* sin(theta) * cos(theta) + 0.2e1 * a5 *  g11 *  h21 * beta * beta * sin(theta) 
* cos(theta) + 0.2e1 * a5 * g21 *  g11 *  g12 *  h22 * beta * beta * sin(theta) 
* cos(theta) - 0.2e1 * a5 *   pow( g12,  3) *  h22 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 *  g12 *  h22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 *  h21 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) - 
0.2e1 * a5 *   pow( g11,  3) *  h21 * beta * beta * pow(sin(theta), 0.2e1) - 
0.2e1 * a5 *  (g11 * g11) *  g12 *  h22 * beta * beta * pow(sin(theta), 0.2e1) + 
0.2e1 * a5 *  (g11 * g11) * g22 *  h22 * beta * beta * sin(theta) * cos(theta) - 
0.2e1 * a5 * g21 * g21 *  g11 *  h21 * beta * beta * pow(cos(theta), 0.2e1) - 
0.2e1 * a5 * g21 *  g11 * g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1) - 
0.2e1 * a5 * g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * 
(g12 * g12) * g22 *  h22 * beta * beta * sin(theta) * cos(theta) + 0.4e1 * a5 * 
(g11 * g11) *  g12 *  h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 
*  (g11 * g11) * g22 *  h22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * 
a5 * g21 *  g11 *  g12 *  h22 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1 
* a5 *  (g11 * g11) * g21 *  h21 * alpha * alpha * cos(theta) * sin(theta) + 
0.4e1 * a5 *   pow( g12,  3) *  h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.4e1 * a5 * g22 * g22 *  g12 *  h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a5 *  g11 *  h21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * 
g21 *  h21 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 *  g12 *  h22 * 
alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 *   pow( g11,  3) *  h21 * 
alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 *  g12 *  h22 * alpha * 
alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 *  g12 * g21 *  h21 * alpha * 
alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 *  g12 *  g11 *  h21 * alpha * 
alpha * cos(theta) * sin(theta) + 0.4e1 * a5 *  g11 *  h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g21 *  g11 * g22 *  h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 *  h22 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a5 *  (g12 * g12) * g21 *  h21 * alpha * alpha * cos(theta) 
* sin(theta) + 0.4e1 * a5 * g21 *  h21 * alpha * alpha * cos(theta) * sin(theta) 
+ 0.4e1 * a5 *  (g12 * g12) *  g11 *  h21 * alpha * alpha * pow(cos(theta), 
0.2e1) + 0.4e1 * a5 * g22 *  h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.8e1 * a5 *  (g12 * g12) * g22 *  h22 * alpha * alpha * cos(theta) * sin(theta) 
+ 0.4e1 * a5 * g21 * g21 *  g11 *  h21 * alpha * alpha * pow(sin(theta), 0.2e1)) 
+ 0.2e1 * alpha * sin(theta) * (0.4e1 * a5 *  (g11 * g11) *  g12 *  h22 * alpha 
* cos(theta) * beta * sin(theta) - 0.2e1 * a5 *  (g11 * g11) * g22 *  h22 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 *  g11 *  g12 *  h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 *  g11 *  g12 *  h22 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 *  (g11 * g11) * g21 *  h21 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 *  (g11 * g11) * g21 *  h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 *  g12 * g21 *  h21 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 *   pow( g12,  3) *  h22 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 *  g12 *  h22 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 *  h21 * alpha * cos(theta) * 
beta * sin(theta) + 0.4e1 * a5 *  g12 *  h22 * alpha * cos(theta) * beta * 
sin(theta) + 0.4e1 * a5 *   pow( g11,  3) *  h21 * alpha * cos(theta) * beta * 
sin(theta) - 0.2e1 * a5 *  g12 *  h22 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.2e1 * a5 *  (g12 * g12) * g21 *  h21 * beta * pow(sin(theta), 0.2e1) * alpha - 
0.2e1 * a5 * g22 *  g12 *  g11 *  h21 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.2e1 * a5 * g22 *  g12 *  g11 *  h21 * beta * pow(sin(theta), 0.2e1) * alpha + 
0.4e1 * a5 *  g11 *  h21 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * 
g11 *  h21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g22 *  h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 *  (g12 * g12) *  g11 *  h21 
* alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 *  (g12 * g12) * g21 * 
h21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 *  h22 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g22 * g22 *  g12 *  h22 * alpha * 
cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g21 *  h21 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 *  h21 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a5 * g21 *  g11 * g22 *  h22 * alpha * cos(theta) * 
beta * sin(theta) - 0.4e1 * a5 * g22 *  h22 * alpha * cos(theta) * beta * 
sin(theta) - 0.4e1 * a5 *  (g12 * g12) * g22 *  h22 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.4e1 * a5 *  (g12 * g12) * g22 *  h22 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.2e1 * a5 *  g11 *  h21 * beta * pow(sin(theta), 0.2e1) * 
alpha + 0.2e1 * a5 *  (g11 * g11) * g22 *  h22 * beta * pow(sin(theta), 0.2e1) * 
alpha - 0.4e1 * a5 * g21 * g21 *  g11 *  h21 * alpha * cos(theta) * beta * 
sin(theta)) + beta * cos(theta) * (- (a4 * a6 * g12 * h22) -  (a4 * a6 * g11 * 
h21) + a5 *  g12 *  h22 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 *  g12 
* g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 *  (g12 * g12) *  g11 * 
h21 * beta * beta * pow(sin(theta), 0.2e1) - a5 *  (g12 * g12) * g21 *  h21 * 
beta * beta * sin(theta) * cos(theta) - a5 * g22 *  g12 *  g11 *  h21 * beta * 
beta * sin(theta) * cos(theta) - a5 * g22 *  h22 * beta * beta * sin(theta) * 
cos(theta) + a5 *  g11 *  h21 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * 
a5 *  (g11 * g11) * g21 *  h21 * beta * beta * sin(theta) * cos(theta) - a5 * 
g12 *  h22 * beta * beta * sin(theta) * cos(theta) - a5 *  g11 *  h21 * beta * 
beta * sin(theta) * cos(theta) - a5 * g21 *  g11 *  g12 *  h22 * beta * beta * 
sin(theta) * cos(theta) + a5 *   pow( g12,  3) *  h22 * beta * beta * 
pow(sin(theta), 0.2e1) + a5 * g22 * g22 *  g12 *  h22 * beta * beta * 
pow(cos(theta), 0.2e1) - a5 * g21 *  h21 * beta * beta * sin(theta) * cos(theta) 
+ a5 * g21 *  h21 * beta * beta * pow(cos(theta), 0.2e1) + a5 *   pow( g11,  3) 
*  h21 * beta * beta * pow(sin(theta), 0.2e1) + a5 *  (g11 * g11) *  g12 *  h22 
* beta * beta * pow(sin(theta), 0.2e1) - a5 *  (g11 * g11) * g22 *  h22 * beta * 
beta * sin(theta) * cos(theta) + a5 * g21 * g21 *  g11 *  h21 * beta * beta * 
pow(cos(theta), 0.2e1) + a5 * g21 *  g11 * g22 *  h22 * beta * beta * 
pow(cos(theta), 0.2e1) + a5 * g22 *  h22 * beta * beta * pow(cos(theta), 0.2e1) 
- 0.2e1 * a5 *  (g12 * g12) * g22 *  h22 * beta * beta * sin(theta) * 
cos(theta));

double J3 = -beta * cos(theta) * (-0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * 
cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * 
cos(theta) * beta * sin(theta) - 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha * 
cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g12 * h22 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.4e1 * a5 * g21 * h21 * alpha * cos(theta) * beta * sin(theta) 
- 0.4e1 * a5 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * 
pow(g11, 0.3e1) * h21 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * 
g12 * h22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g12 * g12 * g21 
* h21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g22 * g12 * g11 * 
h21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g22 * g12 * g11 * h21 
* beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g11 * h21 * alpha * 
cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h21 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 * h22 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cos(theta) * beta * 
sin(theta) + 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) 
* beta - 0.2e1 * a5 * g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 
* a5 * g22 * g22 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * 
a5 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h21 
* beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g21 * g11 * g22 * h22 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * h22 * alpha * 
cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g11 * h21 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cos(theta) * beta 
* sin(theta)) + 0.2e1 * alpha * sin(theta) * (-0.2e1 * a4 * a6 * g12 * h22 - 
0.2e1 * a4 * a6 * g11 * h21 - 0.2e1 * a5 * g12 * h22 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * g12 * g21 * h21 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g12 * g12 * g11 * h21 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g12 * g12 * g21 * h21 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a5 * g22 * g12 * g11 * h21 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a5 * g22 * h22 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g11 * h21 * beta * beta * pow(sin(theta), 0.2e1) + 
0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * beta * sin(theta) * cos(theta) + 
0.2e1 * a5 * g12 * h22 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * 
g11 * h21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * g21 * g11 * g12 
* h22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * pow(g12, 0.3e1) * 
h22 * beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * g22 * g12 * h22 
* beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g21 * h21 * beta * beta * pow(cos(theta), 
0.2e1) - 0.2e1 * a5 * pow(g11, 0.3e1) * h21 * beta * beta * pow(sin(theta), 
0.2e1) - 0.2e1 * a5 * g11 * g11 * g12 * h22 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g21 * g21 * g11 * h21 * beta * beta * pow(cos(theta), 
0.2e1) - 0.2e1 * a5 * g21 * g11 * g22 * h22 * beta * beta * pow(cos(theta), 
0.2e1) - 0.2e1 * a5 * g22 * h22 * beta * beta * pow(cos(theta), 0.2e1) + 0.4e1 * 
a5 * g12 * g12 * g22 * h22 * beta * beta * sin(theta) * cos(theta) + 0.4e1 * a5 
* g11 * g11 * g12 * h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * 
g11 * g11 * g22 * h22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * 
g21 * g11 * g12 * h22 * alpha * alpha * cos(theta) * sin(theta) + 0.8e1 * a5 * 
g11 * g11 * g21 * h21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * 
pow(g12, 0.3e1) * h22 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * 
g22 * g22 * g12 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * 
g11 * h21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g21 * h21 * 
alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h22 * alpha * alpha 
* pow(cos(theta), 0.2e1) + 0.4e1 * a5 * pow(g11, 0.3e1) * h21 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h22 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * g12 * g11 * h21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h21 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a5 * g21 * g11 * g22 * h22 * alpha * alpha * 
pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h22 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a5 * g12 * g12 * g21 * h21 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a5 * g21 * h21 * alpha * alpha * cos(theta) * sin(theta) + 
0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.4e1 * a5 * g22 * h22 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.8e1 * a5 * 
g12 * g12 * g22 * h22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * 
g21 * g21 * g11 * h21 * alpha * alpha * pow(sin(theta), 0.2e1)) + beta * 
cos(theta) * (0.4e1 * a5 * g11 * g11 * g12 * h22 * alpha * cos(theta) * beta * 
sin(theta) - 0.2e1 * a5 * g11 * g11 * g22 * h22 * alpha * pow(cos(theta), 0.2e1) 
* beta - 0.2e1 * a5 * g21 * g11 * g12 * h22 * alpha * pow(cos(theta), 0.2e1) * 
beta + 0.2e1 * a5 * g21 * g11 * g12 * h22 * beta * pow(sin(theta), 0.2e1) * 
alpha - 0.4e1 * a5 * g11 * g11 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) * 
beta + 0.4e1 * a5 * g11 * g11 * g21 * h21 * beta * pow(sin(theta), 0.2e1) * 
alpha - 0.4e1 * a5 * g22 * g12 * g21 * h21 * alpha * cos(theta) * beta * 
sin(theta) + 0.4e1 * a5 * pow(g12, 0.3e1) * h22 * alpha * cos(theta) * beta * 
sin(theta) + 0.2e1 * a5 * g12 * h22 * beta * pow(sin(theta), 0.2e1) * alpha - 
0.4e1 * a5 * g21 * h21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * 
g12 * h22 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * pow(g11, 
0.3e1) * h21 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g12 * h22 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * g12 * g21 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * g12 * g11 * h21 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 * g12 * g11 * h21 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * h21 * alpha * 
cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h21 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g22 * h22 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.4e1 * a5 * g12 * g12 * g11 * h21 * alpha * cos(theta) * beta * 
sin(theta) - 0.2e1 * a5 * g12 * g12 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) 
* beta + 0.2e1 * a5 * g22 * h22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 
* a5 * g22 * g22 * g12 * h22 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * 
a5 * g21 * h21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h21 
* beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g21 * g11 * g22 * h22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g22 * h22 * alpha * 
cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g12 * g12 * g22 * h22 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * g12 * g22 * h22 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g11 * h21 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.2e1 * a5 * g11 * g11 * g22 * h22 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a5 * g21 * g21 * g11 * h21 * alpha * cos(theta) * beta 
* sin(theta));

double K0 = -beta * sin(theta) * (a5 * pow(g21, 0.3e1) * h11 * beta * beta * 
pow(cos(theta), 0.2e1) - a5 * g12 * h12 * beta * beta * sin(theta) * cos(theta) 
+ a5 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g21 * h11 * beta 
* beta * pow(cos(theta), 0.2e1) - a5 * g22 * h12 * g11 * g21 * beta * beta * 
sin(theta) * cos(theta) - a5 * g11 * h11 * beta * beta * sin(theta) * cos(theta) 
+ a5 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * h12 * beta 
* beta * pow(cos(theta), 0.2e1) + a5 * pow(g22, 0.3e1) * h12 * beta * beta * 
pow(cos(theta), 0.2e1) + a5 * g11 * g11 * h11 * g21 * beta * beta * 
pow(sin(theta), 0.2e1) - a5 * g21 * h11 * g12 * g22 * beta * beta * sin(theta) * 
cos(theta) + a5 * g22 * h12 * g21 * g21 * beta * beta * pow(cos(theta), 0.2e1) + 
a5 * g21 * h11 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g21 * 
h11 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * h12 * g11 * g21 * beta 
* beta * pow(sin(theta), 0.2e1) - a5 * g12 * h12 * g21 * g21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * 
sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * g22 * beta * beta * 
pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * g22 * beta * beta * sin(theta) * 
cos(theta) + a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sin(theta), 0.2e1) - 
a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 
* g22 * g22 * beta * beta * sin(theta) * cos(theta) - a4 * a6 * h11 * g21 - a4 * 
a6 * h12 * g22);

double K1 = -beta * sin(theta) * (0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12 * h12 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cos(theta) * beta * 
sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cos(theta) * beta * 
sin(theta) - 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cos(theta) * beta * 
sin(theta) + 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cos(theta) * beta * 
sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cos(theta) * beta * 
sin(theta) + 0.4e1 * a5 * g12 * h12 * alpha * cos(theta) * beta * sin(theta) - 
0.4e1 * a5 * g21 * h11 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * 
g21 * h11 * g12 * g22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 
* h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * h12 * alpha * 
cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * h11 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * cos(theta) * beta 
* sin(theta) - 0.2e1 * a5 * g11 * h11 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha + 
0.4e1 * a5 * g11 * h11 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * 
g11 * h11 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g11 
* h11 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g11 * 
h11 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g11 * h11 
* g22 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * 
g12 * g22 * alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g12 
* g22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g22 * h12 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha * 
cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 * alpha * pow(cos(theta), 0.2e1) * beta 
+ 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha) - 
0.2e1 * alpha * cos(theta) * (a5 * pow(g21, 0.3e1) * h11 * beta * beta * 
pow(cos(theta), 0.2e1) - a5 * g12 * h12 * beta * beta * sin(theta) * cos(theta) 
+ a5 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g21 * h11 * beta 
* beta * pow(cos(theta), 0.2e1) - a5 * g22 * h12 * g11 * g21 * beta * beta * 
sin(theta) * cos(theta) - a5 * g11 * h11 * beta * beta * sin(theta) * cos(theta) 
+ a5 * g12 * h12 * beta * beta * pow(sin(theta), 0.2e1) + a5 * g22 * h12 * beta 
* beta * pow(cos(theta), 0.2e1) + a5 * pow(g22, 0.3e1) * h12 * beta * beta * 
pow(cos(theta), 0.2e1) + a5 * g11 * g11 * h11 * g21 * beta * beta * 
pow(sin(theta), 0.2e1) - a5 * g21 * h11 * g12 * g22 * beta * beta * sin(theta) * 
cos(theta) + a5 * g22 * h12 * g21 * g21 * beta * beta * pow(cos(theta), 0.2e1) + 
a5 * g21 * h11 * g22 * g22 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g21 * 
h11 * beta * beta * sin(theta) * cos(theta) + a5 * g12 * h12 * g11 * g21 * beta 
* beta * pow(sin(theta), 0.2e1) - a5 * g12 * h12 * g21 * g21 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * 
sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * g22 * beta * beta * 
pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * g22 * beta * beta * sin(theta) * 
cos(theta) + a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sin(theta), 0.2e1) - 
a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 
* g22 * g22 * beta * beta * sin(theta) * cos(theta) - a4 * a6 * h11 * g21 - a4 * 
a6 * h12 * g22);

double K2 = -beta * sin(theta) * (-0.2e1 * a5 * pow(g21, 0.3e1) * h11 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g12 * h12 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) - 
0.2e1 * a5 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g22 
* h12 * g11 * g21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * g11 * 
h11 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * beta * 
beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * h12 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a5 * pow(g22, 0.3e1) * h12 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * g11 * h11 * g21 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * g21 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * g22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * beta * sin(theta) * 
cos(theta) + 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * g22 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.4e1 
* a5 * g12 * h12 * g22 * g22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * 
a4 * a6 * h11 * g21 - 0.2e1 * a4 * a6 * h12 * g22 + 0.4e1 * a5 * pow(g22, 0.3e1) 
* h12 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * pow(g21, 0.3e1) * 
h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * alpha * 
alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * g21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a5 * g21 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.4e1 * a5 * g11 * h11 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * 
g12 * h12 * g11 * g21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a5 * 
g11 * h11 * g21 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * 
g11 * h11 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * 
g22 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * 
g12 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * 
alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * g21 * 
alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * g21 * 
alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 * 
alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a5 * g12 * h12 * g22 * g22 * 
alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * alpha * alpha 
* cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * alpha * 
pow(sin(theta), 0.2e1)) - 0.2e1 * alpha * cos(theta) * (0.2e1 * a5 * g22 * h12 * 
g11 * g21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12 * h12 * 
beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * g21 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * h12 * g21 * g21 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * g22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g21 * h11 * g22 * g22 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g11 * g21 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 * 
alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * alpha * 
cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g21 * h11 * alpha * cos(theta) * 
beta * sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta 
- 0.4e1 * a5 * g22 * h12 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * 
g11 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * pow(g21, 0.3e1) 
* h11 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * alpha 
* pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * h11 * alpha * cos(theta) * 
beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cos(theta) * beta * 
sin(theta) + 0.2e1 * a5 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha - 
0.2e1 * a5 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * 
g12 * h12 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g12 
* h12 * g22 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * 
pow(g22, 0.3e1) * h12 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * 
g21 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * h12 * g22 * g22 * 
beta * pow(sin(theta), 0.2e1) * alpha) + beta * sin(theta) * (a5 * pow(g21, 
0.3e1) * h11 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g12 * h12 * beta * 
beta * sin(theta) * cos(theta) + a5 * g11 * h11 * beta * beta * pow(sin(theta), 
0.2e1) + a5 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g22 * h12 
* g11 * g21 * beta * beta * sin(theta) * cos(theta) - a5 * g11 * h11 * beta * 
beta * sin(theta) * cos(theta) + a5 * g12 * h12 * beta * beta * pow(sin(theta), 
0.2e1) + a5 * g22 * h12 * beta * beta * pow(cos(theta), 0.2e1) + a5 * pow(g22, 
0.3e1) * h12 * beta * beta * pow(cos(theta), 0.2e1) + a5 * g11 * g11 * h11 * g21 
* beta * beta * pow(sin(theta), 0.2e1) - a5 * g21 * h11 * g12 * g22 * beta * 
beta * sin(theta) * cos(theta) + a5 * g22 * h12 * g21 * g21 * beta * beta * 
pow(cos(theta), 0.2e1) + a5 * g21 * h11 * g22 * g22 * beta * beta * 
pow(cos(theta), 0.2e1) - a5 * g21 * h11 * beta * beta * sin(theta) * cos(theta) 
+ a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sin(theta), 0.2e1) - a5 * g12 * 
h12 * g21 * g21 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g11 * h11 
* g21 * g21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * h11 * g12 * g22 
* beta * beta * pow(sin(theta), 0.2e1) - a5 * g11 * h11 * g22 * g22 * beta * 
beta * sin(theta) * cos(theta) + a5 * g12 * g12 * h12 * g22 * beta * beta * 
pow(sin(theta), 0.2e1) - a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta) 
- 0.2e1 * a5 * g12 * h12 * g22 * g22 * beta * beta * sin(theta) * cos(theta) - 
a4 * a6 * h11 * g21 - a4 * a6 * h12 * g22);

double K3 = -beta * sin(theta) * (-0.2e1 * a5 * g22 * h12 * g11 * g21 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g12 * h12 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * cos(theta) * beta * 
sin(theta) - 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * cos(theta) * beta * 
sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * cos(theta) * beta * 
sin(theta) - 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * cos(theta) * beta * 
sin(theta) - 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * cos(theta) * beta * 
sin(theta) - 0.4e1 * a5 * g12 * h12 * alpha * cos(theta) * beta * sin(theta) + 
0.4e1 * a5 * g21 * h11 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * 
g21 * h11 * g12 * g22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12 
* h12 * alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g22 * h12 * alpha * 
cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.4e1 * a5 * pow(g21, 0.3e1) * h11 * alpha * cos(theta) * beta 
* sin(theta) + 0.2e1 * a5 * g11 * h11 * alpha * pow(cos(theta), 0.2e1) * beta - 
0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha - 
0.4e1 * a5 * g11 * h11 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * 
g11 * h11 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 
* h11 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g11 * 
h11 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g11 * h11 
* g22 * g22 * alpha * pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g21 * h11 * 
g12 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g11 * h11 * g12 
* g22 * alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g22 * h12 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g22 * h12 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * g21 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * h12 * g22 * g22 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * pow(g22, 0.3e1) * h12 * alpha * 
cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g21 * h11 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.2e1 * a5 * g21 * h11 * alpha * pow(cos(theta), 0.2e1) * beta 
- 0.4e1 * a5 * g12 * h12 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha) - 
0.2e1 * alpha * cos(theta) * (-0.2e1 * a5 * pow(g21, 0.3e1) * h11 * beta * beta 
* pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g12 * h12 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g11 * h11 * beta * beta * pow(sin(theta), 0.2e1) - 
0.2e1 * a5 * g21 * h11 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g22 
* h12 * g11 * g21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * g11 * 
h11 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 * g12 * h12 * beta * 
beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 * g22 * h12 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a5 * pow(g22, 0.3e1) * h12 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g11 * g11 * h11 * g21 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * beta * 
sin(theta) * cos(theta) - 0.2e1 * a5 * g22 * h12 * g21 * g21 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a5 * g21 * h11 * g22 * g22 * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a5 * g21 * h11 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g12 * h12 * g11 * g21 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * beta * sin(theta) * 
cos(theta) + 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g11 * h11 * g12 * g22 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * beta * sin(theta) * 
cos(theta) - 0.2e1 * a5 * g12 * g12 * h12 * g22 * beta * beta * pow(sin(theta), 
0.2e1) + 0.2e1 * a5 * g22 * h12 * beta * beta * sin(theta) * cos(theta) + 0.4e1 
* a5 * g12 * h12 * g22 * g22 * beta * beta * sin(theta) * cos(theta) - 0.2e1 * 
a4 * a6 * h11 * g21 - 0.2e1 * a4 * a6 * h12 * g22 + 0.4e1 * a5 * pow(g22, 0.3e1) 
* h12 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * pow(g21, 0.3e1) * 
h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * alpha * 
alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g22 * h12 * g11 * g21 * alpha * alpha * 
cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * alpha * 
pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * alpha * alpha * cos(theta) * 
sin(theta) + 0.4e1 * a5 * g21 * h11 * alpha * alpha * pow(sin(theta), 0.2e1) + 
0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * alpha * pow(cos(theta), 0.2e1) + 
0.4e1 * a5 * g11 * h11 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * 
g12 * h12 * g11 * g21 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a5 * 
g11 * h11 * g21 * g21 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * 
g11 * h11 * alpha * alpha * pow(cos(theta), 0.2e1) + 0.4e1 * a5 * g11 * h11 * 
g22 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * 
g12 * g22 * alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * 
alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g22 * h12 * g21 * g21 * 
alpha * alpha * pow(sin(theta), 0.2e1) + 0.4e1 * a5 * g12 * h12 * g21 * g21 * 
alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 * 
alpha * alpha * pow(cos(theta), 0.2e1) + 0.8e1 * a5 * g12 * h12 * g22 * g22 * 
alpha * alpha * cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * alpha * alpha 
* cos(theta) * sin(theta) + 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * alpha * 
pow(sin(theta), 0.2e1)) + beta * sin(theta) * (0.2e1 * a5 * g22 * h12 * g11 * 
g21 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12 * h12 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g22 * h12 * g11 * g21 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * h12 * g21 * g21 * alpha * 
cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * g12 * h12 * g22 * alpha * 
cos(theta) * beta * sin(theta) - 0.4e1 * a5 * g21 * h11 * g22 * g22 * alpha * 
cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * g11 * g21 * alpha * 
cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g11 * g11 * h11 * g21 * alpha * 
cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g12 * h12 * alpha * cos(theta) * 
beta * sin(theta) - 0.4e1 * a5 * g21 * h11 * alpha * cos(theta) * beta * 
sin(theta) + 0.2e1 * a5 * g21 * h11 * g12 * g22 * beta * pow(sin(theta), 0.2e1) 
* alpha - 0.2e1 * a5 * g12 * h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 
* a5 * g22 * h12 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * 
h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * pow(g21, 0.3e1) * h11 
* alpha * cos(theta) * beta * sin(theta) - 0.2e1 * a5 * g11 * h11 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g12 * h12 * g21 * g21 * beta * 
pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g11 * h11 * alpha * cos(theta) * 
beta * sin(theta) + 0.4e1 * a5 * g11 * h11 * g21 * g21 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.4e1 * a5 * g11 * h11 * g21 * g21 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.2e1 * a5 * g11 * h11 * g22 * g22 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.2e1 * a5 * g11 * h11 * g22 * g22 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.2e1 * a5 * g21 * h11 * g12 * g22 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.4e1 * a5 * g11 * h11 * g12 * g22 * alpha * cos(theta) * beta * 
sin(theta) + 0.2e1 * a5 * g22 * h12 * beta * pow(sin(theta), 0.2e1) * alpha - 
0.2e1 * a5 * g22 * h12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * 
g12 * h12 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g12 
* h12 * g22 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * 
pow(g22, 0.3e1) * h12 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * 
g21 * h11 * beta * pow(sin(theta), 0.2e1) * alpha - 0.2e1 * a5 * g21 * h11 * 
alpha * pow(cos(theta), 0.2e1) * beta + 0.4e1 * a5 * g12 * h12 * g22 * g22 * 
beta * pow(sin(theta), 0.2e1) * alpha);

double L0 = a4 * a6 * g11 * g21 + a4 * a6 * g12 * g22 + a4 * a6 - a5 * g22 * g22 * beta 
* beta * pow(cos(theta), 0.2e1) - a5 * g21 * g21 * beta * beta * pow(cos(theta), 
0.2e1) + a5 * g11 * g11 * g22 * g22 * beta * beta * sin(theta) * cos(theta) - a5 
* beta * beta * pow(sin(theta), 0.2e1) - a5 * beta * beta * pow(cos(theta), 
0.2e1) - a5 * pow(g12, 0.3e1) * g22 * beta * beta * pow(sin(theta), 0.2e1) - a5 
* g12 * g12 * g11 * g21 * beta * beta * pow(sin(theta), 0.2e1) - a5 * pow(g21, 
0.3e1) * g11 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g22 * g12 * beta * 
beta * pow(cos(theta), 0.2e1) - a5 * g21 * g11 * beta * beta * pow(sin(theta), 
0.2e1) - a5 * pow(g22, 0.3e1) * g12 * beta * beta * pow(cos(theta), 0.2e1) + 
0.2e1 * a5 * g21 * g11 * beta * beta * sin(theta) * cos(theta) - a5 * g22 * g12 
* beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g21 * g11 * g12 * g22 * 
beta * beta * sin(theta) * cos(theta) - a5 * g21 * g11 * beta * beta * 
pow(cos(theta), 0.2e1) + a5 * g12 * g12 * beta * beta * sin(theta) * cos(theta) 
- a5 * g12 * g12 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 * g11 * g11 
* g21 * g21 * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * g22 * g12 * 
beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * beta * beta * sin(theta) * 
cos(theta) - a5 * g22 * g12 * g21 * g21 * beta * beta * pow(cos(theta), 0.2e1) + 
0.2e1 * a5 * g12 * g12 * g22 * g22 * beta * beta * sin(theta) * cos(theta) + a5 
* g12 * g12 * g21 * g21 * beta * beta * sin(theta) * cos(theta) + a5 * g11 * g11 
* beta * beta * sin(theta) * cos(theta) - a5 * g11 * g11 * g12 * g22 * beta * 
beta * pow(sin(theta), 0.2e1) - a5 * pow(g11, 0.3e1) * g21 * beta * beta * 
pow(sin(theta), 0.2e1) + a5 * g21 * g21 * beta * beta * sin(theta) * cos(theta) 
+ a5 * g22 * g22 * beta * beta * sin(theta) * cos(theta) - a5 * g21 * g11 * g22 
* g22 * beta * beta * pow(cos(theta), 0.2e1) - a5 * g11 * g11 * beta * beta * 
pow(sin(theta), 0.2e1);

double L1 = 0.4e1 * a5 * pow(g22, 0.3e1) * g12 * alpha * cos(theta) * beta * sin(theta) 
+ 0.4e1 * a5 * pow(g21, 0.3e1) * g11 * alpha * cos(theta) * beta * sin(theta) + 
0.4e1 * a5 * g21 * g21 * alpha * cos(theta) * beta * sin(theta) + 0.4e1 * a5 * 
g21 * g11 * g12 * g22 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g21 
* g11 * g12 * g22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g21 * 
g11 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g21 * g11 * beta * 
pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * g12 * g11 * g21 * alpha * 
cos(theta) * beta * sin(theta) - 0.4e1 * a5 * pow(g12, 0.3e1) * g22 * alpha * 
cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g11 * g11 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g11 * g11 * beta * pow(sin(theta), 
0.2e1) * alpha + 0.4e1 * a5 * g12 * g12 * g22 * g22 * alpha * pow(cos(theta), 
0.2e1) * beta + 0.4e1 * a5 * alpha * pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 
* beta * pow(sin(theta), 0.2e1) * alpha - 0.4e1 * a5 * g12 * g12 * alpha * 
cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * g12 * g21 * g21 * alpha * 
cos(theta) * beta * sin(theta) + 0.4e1 * a5 * g22 * g12 * alpha * 
pow(cos(theta), 0.2e1) * beta - 0.4e1 * a5 * g22 * g12 * beta * pow(sin(theta), 
0.2e1) * alpha - 0.2e1 * a5 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha 
- 0.4e1 * a5 * g11 * g11 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha + 
0.4e1 * a5 * g11 * g11 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta + 
0.4e1 * a5 * g22 * g22 * alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * 
g12 * g12 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g12 * g12 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g11 * g11 * g22 * g22 * 
alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g11 * g11 * g22 * g22 * 
beta * pow(sin(theta), 0.2e1) * alpha + 0.4e1 * a5 * g21 * g11 * g22 * g22 * 
alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * pow(g11, 0.3e1) * g21 * 
alpha * cos(theta) * beta * sin(theta) + 0.2e1 * a5 * g21 * g21 * alpha * 
pow(cos(theta), 0.2e1) * beta + 0.2e1 * a5 * g22 * g22 * alpha * pow(cos(theta), 
0.2e1) * beta - 0.2e1 * a5 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha - 
0.4e1 * a5 * g11 * g11 * g12 * g22 * alpha * cos(theta) * beta * sin(theta) - 
0.4e1 * a5 * g11 * g11 * alpha * cos(theta) * beta * sin(theta) - 0.4e1 * a5 * 
g12 * g12 * g22 * g22 * beta * pow(sin(theta), 0.2e1) * alpha + 0.2e1 * a5 * g12 
* g12 * g21 * g21 * alpha * pow(cos(theta), 0.2e1) * beta - 0.2e1 * a5 * g12 * 
g12 * g21 * g21 * beta * pow(sin(theta), 0.2e1) * alpha;

double L2 =  (2 * a4 * a6 * g11 * g21) +  (2 * a4 * a6 * g12 * g22) +  (2 * a4 * a6) + 
0.2e1 * a5 *  (g22 * g22) * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 * 
(g21 * g21) * beta * beta * pow(cos(theta), 0.2e1) - 0.2e1 * a5 *  (g11 * g11) * 
(g22 * g22) * beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 * beta * beta * 
pow(sin(theta), 0.2e1) + 0.2e1 * a5 * beta * beta * pow(cos(theta), 0.2e1) + 
0.2e1 * a5 *   pow( g12,  3) *  g22 * beta * beta * pow(sin(theta), 0.2e1) + 
0.2e1 * a5 *  (g12 * g12) *  g11 *  g21 * beta * beta * pow(sin(theta), 0.2e1) + 
0.2e1 * a5 *   pow( g21,  3) *  g11 * beta * beta * pow(cos(theta), 0.2e1) + 
0.2e1 * a5 *  g22 *  g12 * beta * beta * pow(cos(theta), 0.2e1) + 0.2e1 * a5 * 
g21 *  g11 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 *   pow( g22,  3) 
*  g12 * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  g21 *  g11 * beta 
* beta * sin(theta) * cos(theta) + 0.2e1 * a5 *  g22 *  g12 * beta * beta * 
pow(sin(theta), 0.2e1) - 0.4e1 * a5 *  g21 *  g11 *  g12 *  g22 * beta * beta * 
sin(theta) * cos(theta) + 0.2e1 * a5 *  g21 *  g11 * beta * beta * 
pow(cos(theta), 0.2e1) - 0.2e1 * a5 *  (g12 * g12) * beta * beta * sin(theta) * 
cos(theta) + 0.2e1 * a5 *  (g12 * g12) * beta * beta * pow(sin(theta), 0.2e1) - 
0.4e1 * a5 *  (g11 * g11) *  (g21 * g21) * beta * beta * sin(theta) * cos(theta) 
- 0.4e1 * a5 *  g22 *  g12 * beta * beta * sin(theta) * cos(theta) - 0.4e1 * a5 
* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 *  g22 *  g12 *  (g21 * 
g21) * beta * beta * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g12 * g12) *  (g22 
* g22) * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 *  (g12 * g12) * 
(g21 * g21) * beta * beta * sin(theta) * cos(theta) - 0.2e1 * a5 *  (g11 * g11) 
* beta * beta * sin(theta) * cos(theta) + 0.2e1 * a5 *  (g11 * g11) *  g12 * 
g22 * beta * beta * pow(sin(theta), 0.2e1) + 0.2e1 * a5 *   pow( g11,  3) *  g21 
* beta * beta * pow(sin(theta), 0.2e1) - 0.2e1 * a5 *  (g21 * g21) * beta * beta 
* sin(theta) * cos(theta) - 0.2e1 * a5 *  (g22 * g22) * beta * beta * sin(theta) 
* cos(theta) + 0.2e1 * a5 *  g21 *  g11 *  (g22 * g22) * beta * beta * 
pow(cos(theta), 0.2e1) + 0.2e1 * a5 *  (g11 * g11) * beta * beta * 
pow(sin(theta), 0.2e1) - 0.4e1 * a5 * alpha * alpha * pow(cos(theta), 0.2e1) - 
0.4e1 * a5 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a5 *   pow( g22, 
3) *  g12 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a5 *  (g21 * g21) * 
alpha * alpha * pow(sin(theta), 0.2e1) - 0.8e1 * a5 *  (g11 * g11) *  (g21 * 
g21) * alpha * alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  g22 *  g12 * 
alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g11 * g11) * alpha * 
alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g22 * g22) * alpha * alpha * 
pow(sin(theta), 0.2e1) - 0.8e1 * a5 * alpha * alpha * cos(theta) * sin(theta) - 
0.4e1 * a5 *   pow( g21,  3) *  g11 * alpha * alpha * pow(sin(theta), 0.2e1) - 
0.4e1 * a5 *  g22 *  g12 * alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a5 * 
(g12 * g12) *  g11 *  g21 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 
*  g21 *  g11 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *   pow( 
g12,  3) *  g22 * alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g11 * 
g11) * alpha * alpha * cos(theta) * sin(theta) - 0.8e1 * a5 *  (g12 * g12) * 
(g22 * g22) * alpha * alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  (g12 * 
g12) * alpha * alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a5 *  g21 *  g11 *  g12 
*  g22 * alpha * alpha * cos(theta) * sin(theta) - 0.8e1 * a5 *  g22 *  g12 * 
alpha * alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  (g21 * g21) * alpha * 
alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  g21 *  g11 * alpha * alpha * 
pow(sin(theta), 0.2e1) - 0.4e1 * a5 *  (g11 * g11) *  (g22 * g22) * alpha * 
alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  (g12 * g12) *  (g21 * g21) * 
alpha * alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *   pow( g11,  3) *  g21 * 
alpha * alpha * pow(cos(theta), 0.2e1) - 0.4e1 * a5 *  (g22 * g22) * alpha * 
alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  (g11 * g11) *  g12 *  g22 * 
alpha * alpha * pow(cos(theta), 0.2e1) - 0.8e1 * a5 *  g21 *  g11 * alpha * 
alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  g22 *  g12 *  (g21 * g21) * 
alpha * alpha * pow(sin(theta), 0.2e1) - 0.4e1 * a5 *  (g12 * g12) * alpha * 
alpha * cos(theta) * sin(theta) - 0.4e1 * a5 *  g21 *  g11 *  (g22 * g22) * 
alpha * alpha * pow(sin(theta), 0.2e1);
 
//% t-dependent stuff:

//%rf_pose_from_point_tangents_2_fn_t(t);

//% Not used - just to check magnitute of poly.
double int_coef[] = {A0, A1, A2, B0, B1, B2, B3, C0, C1, C2, C3, C4, E0, E1, E2, F0, F1, F2, F3, G0, G1, G2, G3, G4, H0, H1, H2, H3, H4, J0, J1, J2, J3, K0, K1, K2, K3, L0, L1, L2};
common::norm(int_coef);

