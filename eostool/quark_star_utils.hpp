#ifndef QUARK_STAR_UTILS_HPP
#define QUARK_STAR_UTILS_HPP


#include<cmath>
#include<vector>
#include<iostream>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_roots.h>
#include"global_variable_constants.hpp"

using namespace std;
typedef vector< double > state_type;


const double m_qk = 100.; // mass of strange quark in MeV
const double m_br = 930.; // rest baryon mass on stellar surface in MeV

struct quark_params{
	double a4, Beff, egap;
};

inline double Omega_CFL(double mu, double a4, double Beff, double egap){return 3.*pow(m_qk,2.)*pow(mu,2.)/4./M_PI/M_PI- \
	3.*pow(mu,4.)/4./M_PI/M_PI-(1.-12.*log(m_qk/2./mu))/32./M_PI/M_PI*pow(m_qk,4.)+(1.-a4)*3.*pow(mu,4.)/4./M_PI/M_PI-3.*pow(egap,2.)*pow(mu,2.)/M_PI/M_PI+Beff;}
inline double n_CFL(double mu, double a4, double Beff, double egap){return -(3.*pow(m_qk,2.)*mu/2./M_PI/M_PI-3.*pow(mu,3.)/M_PI/M_PI \
	-3.*pow(m_qk,4.)/8./M_PI/M_PI/mu+(1.-a4)*3.*pow(mu,3.)/M_PI/M_PI-6.*pow(egap,2.)*mu/M_PI/M_PI)/3.;}
inline double ph_qks_org(double h, state_type par_eos){return -Omega_CFL(exp(h)*m_br/3., par_eos[0], par_eos[1], par_eos[2]);}
inline double eh_qks_org(double h, double p, state_type par_eos){return -p+3.*n_CFL(exp(h)*m_br/3., par_eos[0], par_eos[1], par_eos[2])*(exp(h)*m_br/3.);}
inline double ph_qks(double h, state_type par_eos){return ph_qks_org(h, par_eos)*MeV4_to_dmls;}
inline double eh_qks(double h, double p, state_type par_eos){return eh_qks_org(h, p/MeV4_to_dmls, par_eos)*MeV4_to_dmls;}

double gammah_qks(double h, double p, double e, state_type par_eos){
//As soon as the unit of p and e are the same, it will be OK,
//no matter it is MeV^4 or dimensionless, because they will cancel each other.
    double mu = exp(h)*m_br/3., n = n_CFL(mu, par_eos[0], par_eos[1], par_eos[2]);
    double mu_dot_n_mu =  -(3.*pow(m_qk,2.)*mu/2./M_PI/M_PI-9.*pow(mu,3.)/M_PI/M_PI \
    +3.*pow(m_qk,4.)/8./M_PI/M_PI/mu+(1.-par_eos[0])*9.*pow(mu,3.)/M_PI/M_PI-6.*pow(par_eos[2],2.)*mu/M_PI/M_PI)/3.;
    return (e+p)*n/(mu_dot_n_mu*p);
}

double minfunc_qks(double h, void *params){
	struct quark_params *p = static_cast<struct quark_params *>(params);
	return Omega_CFL(exp(h)*m_br/3., p->a4, p->Beff, p->egap);
}


double solve_surface_h(double lower_h, double upper_h, state_type par_eos, int max_iter=100){
    //initiate
    int status;
    int iter = 0;
    double guess = (lower_h+upper_h)/2., h_low, h_up, h_tol = 1e-8;
    const gsl_root_fsolver_type *ST;
    gsl_root_fsolver *solver;
    gsl_function Func;
    struct quark_params params = {par_eos[0], par_eos[1], par_eos[2]};
    Func.function = &minfunc_qks;
    Func.params = &params;
    ST = gsl_root_fsolver_brent;
    solver = gsl_root_fsolver_alloc(ST);
    //check root exist
    double f_lower = minfunc_qks(lower_h, &params), f_upper = minfunc_qks(upper_h, &params);
    if ((f_lower < 0.0 && f_upper < 0.0) || (f_lower > 0.0 && f_upper > 0.0)){
        cout<<"f_lower: "<<f_lower<<", f_upper: "<<f_upper<<endl;
        throw std::runtime_error("Do not have root in given region!(in function solve_surface_h)");
    }
    else{gsl_root_fsolver_set(solver, &Func, lower_h, upper_h);}
    //find root
    do{
    	iter++;
    	status = gsl_root_fsolver_iterate(solver);
    	guess = gsl_root_fsolver_root(solver);
    	h_low = gsl_root_fsolver_x_lower(solver);
    	h_up = gsl_root_fsolver_x_upper(solver);
    	status = gsl_root_test_interval(h_low, h_up, 0., h_tol);
        //cout<<h_low<<"\t"<<guess<<"\t"<<h_up<<"\t"<<minfunc_qks(guess, &params)<<endl;
    }
    while (status==GSL_CONTINUE && iter<max_iter);
    gsl_root_fsolver_free(solver);
    return guess;
}


void make_qks_eos_table(double h_begin, double h_end, int n_table, state_type par_eos, string unit="MeV"){
    double e, p, n;
    double h = h_begin;
    double h_step = (h_end-h_begin)/n_table;
    if (unit=="cgs"){cout<<"h"<<"\t\t"<<"n/fm-3"<<"\t\t"<<"e/(g*cm-3)"<<"\t\t"<<"p/(erg*cm-3)"<<endl;}
    else {cout<<"h"<<"\t\t"<<"n/fm-3"<<"\t\t"<<"e/(MeV*fm-3)"<<"\t\t"<<"p/(MeV*fm-3)"<<endl;}
    for (int i=0; i<n_table+1; i++){
        n = n_CFL(exp(h)*m_br/3., par_eos[0], par_eos[1], par_eos[2]);
        p = ph_qks_org(h, par_eos);
        e = eh_qks_org(h, p, par_eos);
        if (unit=="cgs"){cout<<h<<"\t\t"<<n*MeV3_to_ifm3<<"\t\t"<<e*MeV3_to_ifm3*MeVEifm3_to_gEicm3<<"\t\t"<<p*MeV3_to_ifm3*MeVEifm3_to_ergEicm3<<endl;}
        else {cout<<h<<"\t\t"<<n*MeV3_to_ifm3<<"\t\t"<<e*MeV3_to_ifm3<<"\t\t"<<p*MeV3_to_ifm3<<endl;}
        h += h_step;
    }
}


bool print_quark_star_surface(double fh_low, double fh_up, state_type par_eos, string unit="MeV"){
    double h_surf, e_surf, p_surf, n_surf;
    try{
        h_surf = solve_surface_h(fh_low, fh_up, par_eos);
        n_surf = n_CFL(exp(h_surf)*m_br/3., par_eos[0], par_eos[1], par_eos[2]);
        p_surf = ph_qks_org(h_surf, par_eos);
        e_surf = eh_qks_org(h_surf, p_surf, par_eos);
        if (unit=="cgs"){
            cout<<"h"<<"\t\t"<<"n/fm-3"<<"\t\t"<<"e/(g*cm-3)"<<"\t\t"<<"p/(erg*cm-3)"<<endl;
            cout<<h_surf<<"\t\t"<<n_surf*MeV3_to_ifm3<<"\t\t"<<e_surf*MeV3_to_ifm3*MeVEifm3_to_gEicm3<<"\t\t"<<p_surf*MeV3_to_ifm3*MeVEifm3_to_ergEicm3<<endl;
        }
        else {
            cout<<"h"<<"\t\t"<<"n/fm-3"<<"\t\t"<<"e/(MeV*fm-3)"<<"\t\t"<<"p/(MeV*fm-3)"<<endl;
            cout<<h_surf<<"\t"<<n_surf*MeV3_to_ifm3<<"\t\t"<<e_surf*MeV3_to_ifm3<<"\t\t"<<p_surf*MeV3_to_ifm3<<endl;
        }
    }
    catch(exception& e){
        cout<<e.what()<<endl;
        return false;
    }
	return true;
}


#endif
