#ifndef GLOBAL_VARIABLE_CONSTANTS_HPP
#define GLOBAL_VARIABLE_CONSTANTS_HPP


#include<map>
#include<vector>
#include<cmath>
#include<gsl/gsl_const_cgs.h>
#include<iostream>
#include<boost/math/interpolators/pchip.hpp>

using boost::math::interpolators::pchip;
using namespace std;



//----------------------------------------------------------------------------------------------
//                                            Global Variables 
//----------------------------------------------------------------------------------------------



//=======parameterization method(method of joint low density table and high density)=======


bool interp_only = 0;///< interpolation only.
bool piecewise = 0;///< piecewise polytropic expansion, conflict with causal.
bool piecewise_p = 0;///< piecewise polytropic expansion with pressure at different densities.
bool spectral = 0;///< spectral expansion and interpolation.
bool spectral_causal = 0;///< note that Gamma will have different meaning, see Lindblom_18.
bool phase_trans = 0;///< spectral eos causal, phase transition(https://arxiv.org/pdf/1811.10929.pdf, Model2).
bool quark_star = 0;///< quark star model(Alford_2005_ApJ_629_969).
bool adapt_piecewise = 0;///< adaptive rho border, with parameters: delta_rho*n+speed_of_sound*n.
bool cons_cs = 0;///< constant speed of sound.
bool piece_spec_phtr_css = 0;//piecewise+spectral+phase transition+constant speed of sound


//=======Control and storage related=======


//verbose level
bool verbose = 0;///< (command: -v).
bool vverbose = 0;///< (command: -vv).
bool vvverbose = 0;///< (command: -vvv).
//tov mass allowed
double minm_tov = 0.0;///<minimum limit of TOV mass
double maxm_tov = 0.0;///<maximum limit of TOV mass
//constant integration step
bool consid_const_inter_step = false;///< constant integration step
double sg_const_step = 0.;///< single step in constant integration step
//G deviate from standard value
double deviation_G = 0.;///< G'=G*(1+deviation_G)

//important global volatile variables
typedef vector< double > state_type;
map<double, state_type>::iterator coeftb_idx;///< index iterator-get the corresponding interval of h to use corresponding interpolation coefficients.
map<double, state_type> coeftb;///< coefficient of p, e and c, taking h as key.
double mrl_result[3], eos_props[6];
state_type e_borders;
state_type eos_params;
//to store the eos table
state_type eos_table_h;
state_type eos_table_p;
state_type eos_table_e;
state_type eos_table_rho;
vector< pchip<state_type> > eos_table_function_h_base;
vector< pchip<state_type> > eos_table_function_rho_base;
vector< pchip<state_type> > eos_table_function_e_base;

//author suggested parameters
bool check_causal = 1;
double minimum_allowed_h = 0.0316;
double suggested_start_point = minimum_allowed_h+0.1;
//array<double, 4> eos_params = {1.007, 0.404, -0.3709, 0.0696};///< H4(Lindblom-10 TableII).



//----------------------------------------------------------------------------------------------
//                                         Global Constants 
//----------------------------------------------------------------------------------------------


//=======unit transformation related constants=======


//meta unit
double deviation_G_plus1 = deviation_G+1.;
double G = GSL_CONST_CGS_GRAVITATIONAL_CONSTANT*deviation_G_plus1;
const double C = GSL_CONST_CGS_SPEED_OF_LIGHT;
const double Hb = GSL_CONST_CGS_PLANCKS_CONSTANT_HBAR;
const double eV = GSL_CONST_CGS_ELECTRON_VOLT;
const double Ms = GSL_CONST_CGS_SOLAR_MASS;
//natural to cgs
const double MeV_to_ifm = eV/(Hb*C*1.e7);//original: 1e6*eV/(Hb*C*1e13)
const double MeV3_to_ifm3 = pow(MeV_to_ifm, 3);
const double MeVEifm3_to_ergEicm3 = 1.e45*eV;//original: 1e6*eV/(1e-13)^3
const double MeVEifm3_to_gEicm3 = 1.e45*eV/pow(C, 2);//original: MeVEifm3_to_ergEicm3/C^2
//cgs to cactus
double gEicm3_to_dmls = pow(G, 3)*pow(Ms, 2)/pow(C, 6);//original: 1/Ms*(G*Ms/C^2)^3
double ergEicm3_to_dmls = pow(G, 3)*pow(Ms, 2)/pow(C, 8);//original: gEicm3_to_dmls/C^2
//natural to cactus
double MeV4_to_dmls = MeV3_to_ifm3*MeVEifm3_to_ergEicm3*ergEicm3_to_dmls;


//=======integrate related constants=======


//important parameters for integration control
const double sg_step = 1.e-10;///< single step.
const double h_surface = 1e-10;///< surface of the star.
const double h_ig=1e-7;///< where to start integrate TOV equation.
const double h_0 = 2.161113031691454245e-02;///< border of two eos approximation methods. default: 0.0216
//constants (transform factor with CGS unit)
// const double rho_trans = 6.17714470405638e+17;//from Rezzolla unit.py 
// const double p_trans = 5.55174079257738e+38;//from Rezzolla unit.py
// const double r_trans = 1.4766250385;//from Rezzolla unit.py
double length_trans = G*Ms/pow(C, 2);
double rho_trans = pow(C, 6)/(pow(G, 3)*pow(Ms, 2));
double p_trans = rho_trans*pow(C, 2);
double r_trans = length_trans/1.e5;
//core and border
double e_0 = 9.075720079516068750e+13/rho_trans;
double p_0 = 2.979148040306152863e+32/p_trans;
double rho_0 = (e_0+p_0)/exp(h_0);//8.914126511164630e+13/rho_trans;
double rho_sat = 2.7e+14/rho_trans;
double max_possible_e = 1000*rho_sat;
state_type rho_borders{1.*rho_sat, 1.85*rho_sat, 3.7*rho_sat, 7.4*rho_sat};
double eos_table_max_h = 1.5;
char init_interpolation_function_type[3] = {'1', '1', '1'};



#if false
int main(){
    std::cout.precision(10);
    std::cout<<rho_trans<<"\t\t"<<p_trans<<"\t\t"<<r_trans<<std::endl;
    return 0;
}
#endif

#endif
