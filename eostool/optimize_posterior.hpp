#ifndef OPTIMIZE_POSTERIOR_HPP
#define OPTIMIZE_POSTERIOR_HPP


//----------------------------------------------------------------------------------------------
//optimization related functions
//----------------------------------------------------------------------------------------------
#include<gsl/gsl_min.h>
#include<gsl/gsl_errno.h>
#include<iostream>
#include<fstream>
#include"integrate_core.hpp"

using namespace std;
map<double, state_type> h_mrl;//store state of minimization
map<double, double> toolm_h;//tool table of h, taking m as key
typedef map<double, double>::iterator itmap;

double minfunc_opt_post(double h, void *params){
//function to minimize
    if (vvverbose) cout<<"evaluated hc:    "<<h<<endl;
    double *m_aim = static_cast<double*>(params);
    double test_m;
    state_type mrl;
    try {
        state_type result = int_whole_star(h);
        double lambda = cal_lambda(*m_aim/result[1], result[2]);
        test_m = result[0];
        mrl.push_back(test_m); mrl.push_back(result[1]*r_trans); mrl.push_back(lambda);
        h_mrl[h] = mrl;
    }
    catch (exception &){
        cout<<"Invalid value encountered in minfunc_opt_post, set zero to h_mrl table at h="<<h<<endl;
        mrl.push_back(0.), mrl.push_back(0.), mrl.push_back(0.);
        h_mrl[h] = mrl;
        return pow((*m_aim-0.), 2);
    }
    return pow((*m_aim-test_m), 2);
}

double find_aim(double lower, double upper, double m_aim){
//controller of gsl minimization
    //initialize
    int status;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *solver;
    double guess = (lower+upper)/2.0;
    gsl_function F;
    F.function = &minfunc_opt_post;
    F.params = &m_aim;
    T = gsl_min_fminimizer_quad_golden;
    solver = gsl_min_fminimizer_alloc(T);
    gsl_min_fminimizer_set(solver, &F, guess, lower, upper);
    //solve
    do{
        status = gsl_min_fminimizer_iterate(solver);
        guess = gsl_min_fminimizer_x_minimum(solver);
        lower = gsl_min_fminimizer_x_lower(solver);
        upper = gsl_min_fminimizer_x_upper(solver);
        if(verbose) cout<<"lower:"<<lower<<", gusess:"<<guess<<", upper:"<<upper<<endl;
        status = gsl_min_test_interval(lower, upper, 1e-5, 1e-5);
    }while(status == GSL_CONTINUE);
    gsl_min_fminimizer_free(solver);
    return guess;
}

state_type optimize_single_point(double mass){
//return m, r, l, h
    double h;
    itmap id = toolm_h.lower_bound(mass);
    itmap id_l = id, id_u = id;
    id_l--; id_l--; id_u++;
    h_mrl.erase(h_mrl.begin(), h_mrl.end());
    //h = find_aim(toolm_h.begin()->second, toolm_h.end()->second, mass);//only if eos is strange!
    h = find_aim(id_l->second, id_u->second, mass);
    state_type mrlh = h_mrl[h];
    mrlh.push_back(h);
    if (verbose) cout<<"accepted values:   m:"<<mrlh[0]<<"\t\t r:"<<mrlh[1]<<"\t\t l:"<<mrlh[2]<<"\t\t h:"<<mrlh[3]<<endl;
    return mrlh;
}

state_type optimize_single_pair(double m1, double m2){
//return l1, l2, lt, m1, m2
    state_type mrlh1 = optimize_single_point(m1);
    state_type mrlh2 = optimize_single_point(m2);
    state_type l12t_m12;
    l12t_m12.push_back(mrlh1[2]); l12t_m12.push_back(mrlh2[2]);
    l12t_m12.push_back(lambda_tilde(m1, m2, mrlh1[2], mrlh2[2]));
    l12t_m12.push_back(mrlh1[0]); l12t_m12.push_back(mrlh2[0]);
    return l12t_m12;
}

void optimize_file(const string ifname, const string ofname){
//ifname contains mass pairs to be optimized, ofname accept output
    double mass1, mass2;
    state_type l12t_m12;
    int j = 0;
    //read posterior mass_tool
    ifstream fpi(ifname, ifstream::in);
    ofstream fpo(ofname, ofstream::out);
    while(fpi>>mass1>>mass2){
        j++;
        l12t_m12 = optimize_single_pair(mass1, mass2);
        fpo<<l12t_m12[2]<<"\t\t"<<l12t_m12[3]<<"\t\t"<<l12t_m12[4]<<"\t\t"<<l12t_m12[0]<<"\t\t"<<l12t_m12[1]<<endl;
        cout<<j<<"th point"<<endl;
    }
    fpi.close();
    fpo.close();
}


#endif
