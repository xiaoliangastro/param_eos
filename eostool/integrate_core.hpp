#ifndef INTEGRATE_CORE_HPP
#define INTEGRATE_CORE_HPP


/* ===================================================================================
 *	unit system: C = 1, G = 1, M_sun = 1
 *	using adaptive runge_kutta_fehlberg78 method implemented in boost to integrate
 *   Note: 1. x0 will change after calling integrate_adaptive, std::move will clear the table
 *         2. use abs path please
 * ===================================================================================
 */


#include<boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include<boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include<boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include<boost/numeric/odeint/integrate/integrate_const.hpp>
#include<exception>
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<array>
#include"global_variable_constants.hpp"
#include"quark_star_utils.hpp"
#include"data.hpp"


using namespace std;
using namespace boost::numeric::odeint;


//----------------------------------------------------------------------------------------------
//basic integrate controller initiate
//----------------------------------------------------------------------------------------------


struct push_back_state_and_time {
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;
    push_back_state_and_time( std::vector< state_type > &states, std::vector< double > &times ): m_states( states ), m_times( times ) { }
    void operator()( const state_type &x, double t ){
        if (m_states.size()%20==0){
            if (vvverbose) cout<<"        t: "<<t<<", x: ";
            for(int i=0; i<x.size();i++){
                if (vvverbose) cout<<x[i]<<"    ";
                if (x[i]>1e6 or std::isnan(x[i]) or m_times.size()>1e5) {
                    if (verbose) cout<<"Invalid value x["<<i<<"]="<<x[i]<<" encountered at t="<<t<<", with x size: "<<x.size()<<", iteration times: "<<m_times.size()<<endl;
                    throw exception();
                }
            }
            if (vvverbose) cout<<", x size: "<<x.size()<<endl;
        }
        m_states.push_back( x ); m_times.push_back( t );
    }
};

typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
double abs_err = 1.0e-18, rel_err = 1.0e-12, a_x = 1.0, a_dxdt = 1.0;
controlled_stepper_type controlled_stepper(default_error_checker< double, range_algebra, default_operations >( abs_err, rel_err, a_x, a_dxdt ) );//controlled_stepper_type
controlled_stepper_type controlled_stepper_cal_eos(default_error_checker< double, range_algebra, default_operations >( 1e-15, 1e-9, a_x, a_dxdt ) );//controlled_stepper_type
state_type integrate_func(state_type x0, void func(const state_type &, state_type &, double), double start_t, double end_t, bool reverse = true);
state_type int_whole_star(double hc);


//----------------------------------------------------------------------------------------------
//basic functions to calculate the eos properties 
//----------------------------------------------------------------------------------------------


inline double cal_lambda(double C, double Y){return 16*pow(1-2*C, 2.0)*(2+2*C*(Y-1)-Y)/(15*
            (4*pow(C, 3.0)*(13-11*Y+C*(3*Y-2)+2*pow(C, 2.0)*(1+Y))+3*pow(1-2*C, 2.0)*(2-Y+2*C*(Y-1))*log(1-2*C)+2*C*(6-3*Y+3*C*(5*Y-8))));}

inline double lambda_tilde(double m1, double m2, double l1, double l2){return 16.0/13.0*(((12.0*m2+m1)/pow((m1+m2), 5))*pow(m1, 4)*l1+((12.0*m1+m2)/pow((m1+m2), 5))*pow(m2, 4)*l2);}

double eh(double h){
    double co_p = coeftb_idx->second[0], co_e = coeftb_idx->second[1], co_c = coeftb_idx->second[2], co_h = coeftb_idx->first;
    if (coeftb_idx!= coeftb.begin()) return co_e*pow((co_e+co_p)/co_p*exp((co_c-1.0)/co_c*(h-co_h))-co_e/co_p, 1.0/(co_c-1.0));
    else return co_e*pow(co_e/co_p*(exp(2.0*h/5.0)-1.0), 1.5);
}

inline double ph(double h){return coeftb_idx->second[0]*pow(eh(h)/coeftb_idx->second[1], coeftb_idx->second[2]);}

inline double gammah_interpolation(double p, double e){return coeftb_idx->second[2]*(p/e+1.0); }

/** @brief Calculate gamma in the piece_spec_phtr_css model. */
double gammaeh_pspc(double h, double p, double e, state_type a){
//expected a: rhob1, gamma1, rhob2, h0, g0, g1, g2, rhob3, gamma_pt, cs
    double rho = (e+p)/exp(h), gamma, loggamma = 0.;
    if (rho<a[0]) gamma = a[1];
    else if (rho<a[2]) {
        for (int j = 0; j<3; j++){ loggamma += a[4+j]*pow(log(h/a[3]), j); }
        gamma = (e+p)/((exp(loggamma)+1.)*p);
    }
    else if (rho<a[7]) gamma = a[8];
    else gamma = (e+p)*pow(a[9], 2)/p;
    return gamma;
}

/** @brief Calculate gamma in the spectral expansion model. */
double gammah_spectral(double h, double p, double e, state_type a, double h_init=h_0){
    double loggamma = 0.0, ret;
    for(int j = 0; j<a.size(); j++) loggamma += a[j]*pow(log(h/h_init), j);
    ret = exp(loggamma);
    if (spectral_causal) ret = (e+p)/((ret+1.)*p);
    return ret;
}

/** @brief Calculate gamma in the gamma piecewise model. */
double gammae_piecewise(double e, state_type a){
    auto upper = std::lower_bound(e_borders.begin()+1, e_borders.end(), e);
    int idx = std::distance(e_borders.begin(), upper)-1;
    if (upper==e_borders.end()) {idx -= 1;}
    return a[idx];
}

/** @brief Calculate gamma in the pressure piecewise model. */
double gammae_piecewise_p(double e, state_type a){
    auto upper = std::lower_bound(e_borders.begin()+1, e_borders.end(), e);
    int idx = std::distance(e_borders.begin(), upper)-1;
    if (upper==e_borders.end()) {idx -= 1;}
    if(idx==0) return log(a[0]/p_0)/log(rho_borders[0]/rho_0);
    else return log(a[idx]/a[idx-1])/log(rho_borders[idx]/rho_borders[idx-1]);
}

/** @brief Calculate gamma in the phase transition model. */
double gammah_phase_trans(double p, double e, state_type a){
    auto upper = std::lower_bound(e_borders.begin()+1, e_borders.end(), e);
    int idx = std::distance(e_borders.begin(), upper)-1;
    if(idx<=2) return a[idx];
    else if(idx==3) return e/p+1.;
    else return (e/p+1.)*(1./3.);
}

/** @brief Calculate gamma in the constant speed of sound model.
  * @attention Define v^2=dp/drho, and p=K*rho^{gamma}.
*/
double gammae_cons_cs(double p, double e, state_type a){
    return (e+p)*pow(a[0], 2)/p;
}

/** @brief Choose a proper function to calculate the adiabatic index.
  * @param h Enthalpy at where to calculate the gamma.
  * @param p Pressure at where to calculate the gamma.
  * @param e Energy density at where to calculate the gamma.
  * @param a State vector to store parameterization parameters.
*/
double gamma_which(double h, double p, double e, state_type a){
    if (piecewise) return gammae_piecewise(e, a);
    else if (cons_cs) return gammae_cons_cs(p, e, a);
    else if (piece_spec_phtr_css) return gammaeh_pspc(h, p, e, a);
    else if (piecewise_p or adapt_piecewise) return gammae_piecewise_p(e, a);
    else if (spectral or spectral_causal) return gammah_spectral(h, p, e, a);
    else if (phase_trans) return gammah_phase_trans(p, e, a);
    else { cout<<"Unknown method of jointing low density eos table and high density one."<<endl; exit(0);}
}


//----------------------------------------------------------------------------------------------
//ODE group to be solved, x = (p, e, m, r, y)
//----------------------------------------------------------------------------------------------


inline double dpdh(double p, double e){return e+p;}

inline double dedh(double p, double e, double Gamma){return pow((e+p), 2.0)/(p*Gamma);}

inline double dmdh(double p, double e, double m, double r){return -(4*M_PI*e*pow(r, 3.0)*(r-2*m))/(m+4*M_PI*pow(r, 3.0)*p);}

inline double drdh(double p, double e, double m, double r){return -(r*(r-2*m))/(m+4*M_PI*pow(r, 3.0)*p);}

double dydh(double p, double e, double m, double r, double y, double Gamma){
    double dydh1 = ((r-2*m)*(y+1)*y)/(m+4*M_PI*pow(r, 3.0)*p)+y+((m-4*M_PI*pow(r, 3.0)*e)*y+4*M_PI*pow(r, 3.0)*(5*e+9*p)-6*r)/(m+4*M_PI*pow(r, 3.0)*p);
    double dydh2 = (4*M_PI*pow(r, 3.0)*pow(e+p, 2))/((m+4*M_PI*pow(r, 3.0)*p)*p*Gamma)-4*(m+4*M_PI*pow(r, 3.0)*p)/(r-2*m); 
    return (dydh1+dydh2);
}

//inline double dbmdh(double p, double rho, double m, double r){return -(4*M_PI*rho*pow(r, 3.0)*(r-2*m))/(m+4*M_PI*pow(r, 3.0)*p);}
//inline double dIdh(double p, double rho, double m, double r){return -(8./3.*M_PI*rho*pow(r, 5.0)*(r-2*m))/(m+4*M_PI*pow(r, 3.0)*p);}


//----------------------------------------------------------------------------------------------
//second order utilities to integrate
//----------------------------------------------------------------------------------------------


void cal_eos(const state_type &x, state_type &dxdt, double h){
    double gamma = gamma_which(h, x[0], x[1], eos_params);
    dxdt[0] = dpdh(x[0], x[1]);
    dxdt[1] = dedh(x[0], x[1], gamma);
    dxdt[2] = x[2]/(gamma*x[0])*dxdt[0];
    //cout<<h<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<dxdt[0]<<" "<<dxdt[1]<<" "<<dxdt[2]<<" "<<gamma<<endl;
}

void inner_to_outer(const state_type &x, state_type &dxdt, double h){
    double p_h = eos_table_function_h_base[0](h), e_h = eos_table_function_h_base[1](h);
    double gamma = gamma_which(h, p_h, e_h, eos_params);
    dxdt[0] = dmdh(p_h, e_h, x[0]*deviation_G_plus1, x[1]);
    dxdt[1] = drdh(p_h, e_h, x[0]*deviation_G_plus1, x[1]);
    dxdt[2] = dydh(p_h, e_h, x[0]*deviation_G_plus1, x[1], x[2], gamma);
}

void inner_to_outer_interpolation(const state_type &x, state_type &dxdt, double h){
    double p_h, e_h, gamma;
    if (quark_star){
        p_h = ph_qks(h, eos_params); e_h = eh_qks(h, p_h, eos_params);
        gamma = gammah_qks(h, p_h, e_h, eos_params);
    }
    else{
        coeftb_idx = --coeftb.upper_bound(h);
        p_h = ph(h); e_h = eh(h);
        gamma = gammah_interpolation(p_h, e_h);
    }
    dxdt[0] = dmdh(p_h, e_h, x[0]*deviation_G_plus1, x[1]);
    dxdt[1] = drdh(p_h, e_h, x[0]*deviation_G_plus1, x[1]);
    dxdt[2] = dydh(p_h, e_h, x[0]*deviation_G_plus1, x[1], x[2], gamma);
    //cout<<h<<" "<<x[0]<<" "<<x[1]*r_trans<<endl;
}

void outer_to_inner(const state_type &x, state_type &dxdt, double h){
    dxdt[0] = dpdh(x[0], x[1]);
    dxdt[1] = dedh(x[0], x[1], gamma_which(h, x[0], x[1], eos_params));
}


//----------------------------------------------------------------------------------------------
//third order utilities to integrate the whole star
//----------------------------------------------------------------------------------------------


/** @brief Integrated portable function to integrate ODE. */
state_type integrate_func(state_type x0, void func(const state_type &, state_type &, double), double start_t, double end_t, bool reverse){
    vector<state_type> x_o;//x for output
    state_type h_o;//h for output
    double dt;
    //lower bound must smaller than upper bound, so special care should be taken
    try{
        if (consid_const_inter_step){
            if (reverse) dt = -sg_const_step; else dt = sg_const_step;
            integrate_const(error_stepper_type(), func, x0, start_t, end_t, dt, push_back_state_and_time(x_o, h_o));
        }
        else{
            if (reverse) dt = -sg_step; else dt = sg_step;
            integrate_adaptive(controlled_stepper, func, x0, start_t, end_t, dt, push_back_state_and_time(x_o, h_o));
        }
    }
    catch (exception &){
        throw;
    }
    //debug
    if(vverbose){
        cout<<"iter steps:"<<h_o.size()<<endl<<endl;
        cout<<"h_i"<<"\t\t"<<"p/m"<<"\t\t"<<"e/r"<<"\t\t"<<"y"<<endl;
        for(int j = 0; j<h_o.size(); j++){
            cout<<h_o[j]<<"\t\t"; 
            for (int i = 0; i<x_o[j].size(); i++) cout<<x_o[j][i]<<"\t\t";
            cout<<endl;
        }
    }
    return x_o.back();
}

/** @brief Integrate the most central core part of the compact star to avoid singular point. */
state_type initiate_core(double hc, double hig=h_ig){
    double ec, pc, Gamma_c;
    state_type init_x(3);
    if (interp_only){coeftb_idx = --coeftb.upper_bound(hc); ec = eh(hc); pc = ph(hc); Gamma_c = gammah_interpolation(pc, ec);}
    else if (quark_star){pc = ph_qks(hc, eos_params); ec = eh_qks(hc, pc, eos_params); Gamma_c = gammah_qks(hc, pc, ec, eos_params);}
    else {
        try{
            pc = eos_table_function_h_base[0](hc), ec = eos_table_function_h_base[1](hc);
            Gamma_c = gamma_which(hc, pc, ec, eos_params);
        }
        catch (exception & except){
            if (verbose) cout<<"Invalid value encountered in initiate_core: "<<except.what()<<endl;
            throw;
        }
    }
    double r1 = pow(3.0/(2*M_PI*(ec+3*pc)), 0.5);
    double r3 = -r1/(4*(ec+3*pc))*(ec-3*pc-3*pow(ec+pc, 2.0)/(5*pc*Gamma_c));
    double m3 = 4*M_PI/3.0*ec*pow(r1, 3.0);
    double m5 = 4*M_PI*pow(r1, 3.0)*(r3*ec/r1-pow(ec+pc, 2.0)/(5*pc*Gamma_c));
    double y2 = -6*(ec/3+11*pc+pow(ec+pc, 2.0)/(pc*Gamma_c))/(7*(ec+3*pc));
    init_x[0] = m3*pow(hig, 1.5)+m5*pow(hig, 2.5);
    init_x[1] = r1*pow(hig, 0.5)+r3*pow(hig, 1.5);
    init_x[2] = 2+y2*(hig);
    if (vverbose) cout<<"M_init:"<<init_x[0]<<"    R_init:"<<init_x[1]<<endl;
    return init_x;
}

/** @brief Integrate the whole structure of the star.
  * @details Give h_c, get m, r and y.
  * @note Radius r in this function should multiply half Schwarzschild radius of the sun.
  * @param hc enthalpy at the center of the compact star.
*/
state_type int_whole_star(double hc){
    state_type  x_result(3);
    double hc_lowerbound = h_0+0.01;
    if ((not quark_star) and (hc<hc_lowerbound)){
        cout<<"Invalid value encountered in int_whole_star, hc="<<hc<<" < "<<hc_lowerbound<<" not valid"<<endl;
        throw exception();
    }
    try{
        state_type init_x = initiate_core(hc);
        if (interp_only){x_result = integrate_func(init_x, inner_to_outer_interpolation, hc-h_ig, h_surface);}
        else if (quark_star){
            double h_surf_qks = eos_params[3];
            double p_surf_qks = ph_qks(h_surf_qks, eos_params);
            if (p_surf_qks<1e-300){ //avoid divided by zero in calculating lambda
                eos_params[3] += 1e-8;
                h_surf_qks = eos_params[3];
                p_surf_qks = ph_qks(h_surf_qks, eos_params);
            }
            double e_surf_qks = eh_qks(h_surf_qks, p_surf_qks, eos_params);
            x_result = integrate_func(init_x, inner_to_outer_interpolation, hc-h_ig, h_surf_qks);
            double average_e = x_result[0]/(4.*M_PI*pow(x_result[1],3)/3.);
            x_result[2] -= 3*e_surf_qks/average_e;
        }
        else {
            state_type x0(init_x), x_inter_result(3);
            x_inter_result = integrate_func(x0, inner_to_outer, hc-h_ig, h_0);
            x_result = integrate_func(x_inter_result, inner_to_outer_interpolation, h_0, h_surface);
        }
    }
    catch (exception &){
        if (verbose) cout<<"Invalid value encountered in int_whole_star("<<hc<<"), set results to 0"<<endl;
        x_result[0]=0;x_result[1]=0;x_result[2]=0;
        throw;
    }
    return x_result;
}


#endif