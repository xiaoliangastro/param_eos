#ifndef PARAM_UTILS_HPP
#define PARAM_UTILS_HPP


#include"integrate_core.hpp"
#include<boost/math/interpolators/makima.hpp>
#include<cstdlib>

using boost::math::interpolators::makima;


//----------------------------------------------------------------------------------------------
//check eos & make tool table
//----------------------------------------------------------------------------------------------


//important global parameters

extern "C"{
    double* get_mrl(double hc);
    bool cool_eos(double *hc);
    bool check_mmax(double *hc, double *M_max, double h_start=suggested_start_point+0.08, bool cp_func_type_less=true, double jump_size=0.1, bool check_ok=true);
    bool check_mmax_pt_two_branch(double* h1, double* h2, double* h3, double* h_max, double* m1, double* m2, double* m3, double* m_max, double rho_tr, double drho);
    void make_eos_table(double h_start, double h_end, double dh, int precision, char *f_name, const char *unit="cgs", int get_type=1);
    double interp_pe_likelihood(double *e_org, double *p_org, int size_org);
    bool get_unknowns_from_knowns(double known1, double known2, double *unknown1, double *unknown2, double *h_max, int get_type);
    double* find_eos_properties(double known_aim, int find_type);
    bool find_closest(double m_aim, double h_i, double h_max, double *h_closest, double *lambda, bool use_user_start_point=false);
    bool find_closest_with_maxm_known(double known_aim, double h_i, double h_max, double *h_closest, double *unknown, int get_type, bool use_user_start_point=false);
}


/** @brief Python interface function for integrating the whole structure of the star.
  * @details Give h_c, get M/M_sun, R/km and L, which are stored in global variable double *mrl_result.
  * @attention Use this for python control only, because double* is convenient for python, and global variable is needed,  
          which may cause confuse problem if used in c++ main, so please use int_whole_star instead for c++, 
          because that do not need to return global variable.
  * @param hc enthalpy at the center of the compact star.
*/
double* get_mrl(double hc){
    state_type x_result(3);
    try {
        x_result = int_whole_star(hc);
        mrl_result[0] = x_result[0], mrl_result[1] = x_result[1]*r_trans;
        mrl_result[2] = cal_lambda(x_result[0]/x_result[1], x_result[2]);
    }
    catch (exception &){
        if (verbose) cout<<"Invalid value encountered in get_mrl("<<hc<<"), set results to 0"<<endl;
        mrl_result[0] = 0, mrl_result[1] = 0, mrl_result[2] = 0;
    }
    return mrl_result;
}

bool integrate_eos(double h_start, double h_end, double dh, state_type *h_tb, state_type *p_tb, state_type *e_tb, \
                   state_type *rho_tb, state_type *gamma_tb, state_type *v_sq_tb, int get_type){
    bool success = true;
    double h = h_end, p, e, gamma, rho, v_square;
    if (interp_only){
        for (int i = 0; i<floor((h_end-h_start)/dh); i++){
            coeftb_idx = --coeftb.upper_bound(h);
            p = ph(h), e = eh(h), rho = (p+e)/exp(h);
            gamma = coeftb_idx->second[2]*(p/e+1.0);
            v_square = (p*gamma)/(e+p);
            h_tb->push_back(h); p_tb->push_back(p); e_tb->push_back(e);
            rho_tb->push_back(rho); gamma_tb->push_back(gamma); v_sq_tb->push_back(v_square); 
            h -= dh;
        }
    }
    else if (quark_star){
        for (int i = 0; i<floor((h_end-h_start)/dh); i++){
            p = ph_qks(h, eos_params), e = eh_qks(h, p, eos_params);
            gamma = gammah_qks(h, p, e, eos_params);
            v_square = (p*gamma)/(e+p), rho = (p+e)/exp(h);
            h_tb->push_back(h); p_tb->push_back(p); e_tb->push_back(e);
            rho_tb->push_back(rho); gamma_tb->push_back(gamma); v_sq_tb->push_back(v_square);
            h -= dh;
        }       
    }
    else {
        vector <state_type> per_o;
        state_type per_c, h_o;
        state_type x0(3); x0[0] = p_0; x0[1] = e_0; x0[2] = rho_0;
        try{
            per_c = integrate_func(x0, cal_eos, h_0, h_end, false);
            for (int i = 0; i<floor((h_end-h_start)/dh); i++){
                if (h>h_0){
                    if (get_type==1){
                        p = eos_table_function_h_base[0](h), e = eos_table_function_h_base[1](h), rho = eos_table_function_h_base[2](h);
                    }
                    else{
                        integrate_adaptive(controlled_stepper, cal_eos, per_c, h, h-dh, -1e-6, push_back_state_and_time(per_o, h_o));
                        p = per_o.back()[0], e = per_o.back()[1], rho = per_o.back()[2];
                    }
                    gamma = gamma_which(h, p, e, eos_params);
                }
                else {
                    coeftb_idx = --coeftb.upper_bound(h);
                    p = ph(h), e = eh(h), rho = (p+e)/exp(h);
                    gamma = coeftb_idx->second[2]*(p/e+1.0);
                }
                v_square = (p*gamma)/(e+p);
                h_tb->push_back(h); p_tb->push_back(p); e_tb->push_back(e);
                rho_tb->push_back(rho); gamma_tb->push_back(gamma); v_sq_tb->push_back(v_square);
                h -= dh;
            }
        }
        catch (exception & except){
            cout<<"failed in integrate_eos, "<<except.what()<<endl;
            success = false;
        }
    }
    return success;
}

/** @brief Python interface function for finding eos properties from h(get_type=1), 
           rho(get_type=2) or e(get_type=3).
  * @params get_type must equal to init_interpolation_function_type-1.
*/
double* find_eos_properties(double known_aim, int find_type){
    // h, p, e, rho, gamma, vs
    double h, p, e, rho;
    #define fail_process_find_eos_prop() for (int i=0; i<6; i++) {eos_props[i] = 0.;} return eos_props;
    try{
        if (find_type==1){
            h = known_aim;
            p = eos_table_function_h_base[0](h);
            e = eos_table_function_h_base[1](h);
            rho = eos_table_function_h_base[2](h);
        }
        else if (find_type==2){
            e = known_aim;
            h = eos_table_function_e_base[0](e);
            p = eos_table_function_e_base[1](e);
            rho = eos_table_function_e_base[2](e);
        }
        else if (find_type==3){
            rho = known_aim;
            h = eos_table_function_rho_base[0](rho);
            p = eos_table_function_rho_base[1](rho);
            e = eos_table_function_rho_base[2](rho);
        }
        else {
            cout<<"Unknown find_type in find_eos_properties, return zero"<<endl;
            fail_process_find_eos_prop();
        }
    }
    catch (exception & except){
        cout<<"failed in find_eos_properties: "<<except.what()<<endl;
        fail_process_find_eos_prop();
    }
    double gamma = gamma_which(h, p, e, eos_params), v_square = (p*gamma)/(e+p);
    eos_props[0] = h, eos_props[1] = p, eos_props[2] = e;
    eos_props[3] = rho, eos_props[4] = gamma, eos_props[5] = v_square;
    return eos_props;
}

double interp_pe_likelihood(double *e_org, double *p_org, int size_org){
    double likelihood = 0., increase = 0.;
    try {
        for (int i=0; i<size_org; i++) {
            increase = pow(log(p_org[i]/eos_table_function_e_base[1](e_org[i])), 2.);
            if (not std::isnan(increase)) likelihood += increase;
        }
    }
    catch (exception & except) {
        cout<<"interpolate error: "<<except.what()<<endl;
        likelihood = +1.e100;
    }
    return likelihood;
}

double interp_pe_likelihood_precise(double *e_org, double *p_org, int size_org, double h_start, double h_end, double dh){
    double likelihood = 0., increase = 0.;
    state_type h_tb, p_tb, e_tb, rho_tb, gamma_tb, v_sq_tb;
    state_type e_buff, p_buff, eorg_buff;
    for (int j=0; j<size_org; j++) eorg_buff.push_back(e_org[j]);
    bool int_eos_success = integrate_eos(h_start, h_end, dh, &h_tb, &p_tb, &e_tb, &rho_tb, &gamma_tb, &v_sq_tb, 2);
    if (not int_eos_success) {
        cout<<"integrate eos error at interp_pe_likelihood."<<endl;
        return +1.e100;
    }
    state_type::iterator e_tb_ub = std::upper_bound(e_tb.begin(), e_tb.end(), e_org[0], [](double val, double element){ return val>element;});
    if (e_tb_ub != e_tb.end()) e_tb_ub += 1;
    int distl = std::distance(e_tb.begin(), e_tb_ub);
    std::reverse_copy(e_tb.begin(), e_tb_ub, std::back_inserter(e_buff));
    std::reverse_copy(p_tb.begin(), p_tb.begin()+distl, std::back_inserter(p_buff));
    state_type::iterator e_tb_lb = std::lower_bound(eorg_buff.begin(), eorg_buff.end(), e_buff[e_buff.size()-1]);
    int distu = std::distance(eorg_buff.begin(), e_tb_lb);
    try {
        auto f_interp = makima(std::move(e_buff), std::move(p_buff)); //'move' will clear the table
        for (int i=0; i<distu; i++) { 
            increase = pow(log(p_org[i]/f_interp(e_org[i])), 2.);
            if (not std::isnan(increase)) likelihood += increase;
        }
    }
    catch (exception & except) {
        cout<<"interpolate error: "<<except.what()<<endl;
        likelihood = +1.e100;
    }
    return likelihood;
}

void make_eos_table(double h_start, double h_end, double dh, int precision, char *f_name, const char *unit, int get_type){
//check whether e, p have been correctly interpolated, cgs: p, e, rho
    state_type h_tb, p_tb, e_tb, rho_tb, gamma_tb, v_sq_tb;
    if (verbose) cout<<"Output EoS table to file: "<<f_name<<endl;
    ofstream fp(f_name);
    fp.precision(precision);
    integrate_eos(h_start, h_end, dh, &h_tb, &p_tb, &e_tb, &rho_tb, &gamma_tb, &v_sq_tb, get_type);
    if (strcmp(unit, "cgs")==0) fp<<"h"<<"\t\t"<<"p/dyn*cm-2"<<"\t\t"<<"e/erg*cm-3"<<"\t\t"<<"(be careful)rho/g*cm-3"<<"\t\t"<<"gamma"<<"\t\t"<<"v_square"<<endl;
    else fp<<"h"<<"\t\t\t"<<"p"<<"\t\t\t"<<"e"<<"\t\t\t"<<"rho"<<"\t\t\t"<<"gamma"<<"\t\t\t"<<"v_square"<<endl;
    for (int i=0; i<h_tb.size(); i++){
        if (strcmp(unit, "cgs")==0) {fp<<h_tb[i]<<"\t\t"<<p_tb[i]*p_trans<<"\t\t"<<e_tb[i]*rho_trans<<"\t\t"<<rho_tb[i]*rho_trans<<"\t\t"<<gamma_tb[i]<<"\t\t"<<v_sq_tb[i]<<endl;}
        else fp<<h_tb[i]<<"\t\t"<<p_tb[i]<<"\t\t"<<e_tb[i]<<"\t\t"<<rho_tb[i]<<"\t\t"<<gamma_tb[i]<<"\t\t"<<v_sq_tb[i]<<endl;
    }
    fp.close();
}

void make_tool_table(double minh, double maxh, double dh){
//make m,h table to accelerate optimization
    double h = minh;
    state_type result_x(3);
    string tool_fname = "../../tool/eos_table/test_hmrl.txt";
    ofstream ofp(tool_fname);
    for (int i = 0; i<floor((maxh-minh)/dh); i++){
        try{
            result_x = int_whole_star(h);
            double tidal_deform = cal_lambda(result_x[0]/result_x[1], result_x[2]);
            ofp<<h<<"\t\t"<<result_x[0]<<"\t\t"<<result_x[1]*r_trans<<"\t\t"<<tidal_deform<<endl;
            //tidal_deform*pow(result_x[0]*r_trans, 5), tidal_deform, 3./2.*tidal_deform*pow(result_x[0]/result_x[1], 5)
        }
        catch (exception &){
            cout<<"Invalid value encountered in make_tool_table, h="<<h<<endl;
        }
        h += dh;
    }
    ofp.close();
}


//----------------------------------------------------------------------------------------------
//constraining eos related
//----------------------------------------------------------------------------------------------


bool cool_eos(double *hc){
//check whether gamma(h)<7
    bool valid = true;
    double h_lowest = suggested_start_point;
    double h = *hc, dh=0.001, gamma, p, e, v_square;
    vector <state_type> pe_o;
    state_type pe_c, h_o;
    try{
        for (int i = 0; i<floor((*hc-h_lowest)/dh); i++){
            if (interp_only){
                coeftb_idx = --coeftb.upper_bound(h);
                p = ph(h), e = eh(h);
                gamma = gammah_interpolation(p, e);
            }
            else if (quark_star){
                p = ph_qks(h, eos_params), e = eh_qks(h, p, eos_params);
                gamma = gammah_qks(h, p, e, eos_params);
            }
            else {
                p = eos_table_function_h_base[0](h), e = eos_table_function_h_base[1](h);
                gamma = gamma_which(h, p, e, eos_params);
            }
            v_square = (p*gamma)/(e+p);
            valid = (v_square<1.) and (gamma<10); // must be this way, double checked
            if (verbose) cout<<"h: "<<h<<"   p: "<<p<<"   e: "<<e<<"    gamma: "<<gamma<<"   v^2: "<<v_square<<endl;
            if (valid) break;
            h -= dh;
        }
    }
    catch (exception & except){
        cout<<"Invalid value encountered in cool eos: "<<except.what()<<endl;
        valid = false;
    }
    if (valid) *hc = h;
    return valid;
}

/**
* @brief Find non-rotating maximum supported mass of an EoS.
* @details
* The challenge is the "noise" or none monotonicity property, the solution is:
* 1. A bigger guard to check again, 2. Judge whether a step before is larger.
* error: ~0.005 M_sun
* @param hc  Pointer to the central enthalpy, return -1 if finding error.
* @param M_max Pointer to the maximum mass, return 0 if finding error.
* @retval bool whether finding process correctly worked.
*/
bool check_mmax(double *hc, double *M_max, double h_start, bool cp_func_type_less, double jump_size, bool check_ok){
    auto cp_func = [=](double cp1, double cp2){
        if (cp_func_type_less) return cp1<cp2;
        else return cp1>cp2;
    };
    bool jump_back=false, force_jump_back=false;
    state_type x_result(3), x_guard(3), x_check_again(3);
    double h_max=eos_table_max_h, h_guard=4e-3, h_check_again=0.01;
    double h = h_start;
    double local_minimum_allowed_h = h_start-1e-4;
    double M=0., M_guard=0., M_check_again=0., M_temp=0.;
    #define fail_process_check_mmax(message) \
        *hc = -1; *M_max=0.; \
        if (verbose) { \
            cout<<message<<", with h="<<h<<", h_guard="<<h+h_guard<<", h_check_again="<<h+h_check_again; \
            cout<<", M="<<M<<", M_guard="<<M_guard<<", M_check_again="<<M_check_again<<endl;\
        } \
        return false;
    //find m_max
    if(verbose) cout<<endl<<"find maximum mass, with h_min: "<<local_minimum_allowed_h<<", h_max: "<<h_max<<endl;
    while(h>local_minimum_allowed_h and h<h_max and jump_size>=5e-4){
        if (jump_back) {
            if (force_jump_back) h -= jump_size;
            h -= jump_size;
            jump_size /= 5.;
            h += jump_size;
        }
        else h += jump_size;
        try{
            x_result = int_whole_star(h);
            M = x_result[0];
            x_guard = int_whole_star(h+h_guard);
            M_guard = x_guard[0];
        }
        catch (exception &){
            if (verbose) cout<<"Invalid value encountered in check_mmax, with h="<<h<<", or h="<<h+h_guard<<endl;
            M=0; M_guard=0;
        }
        if (verbose) cout<<"h: "<<h<<"     M:"<<M<<"     M_guard:"<<M_guard<<"     "<<endl;
        jump_back = (M==0 or M_guard==0 or std::isnan(M_guard) or std::isnan(M));
        force_jump_back = (M==0 or M_guard==0);
        if ((not jump_back) and (cp_func(M_guard, M) or M_guard==M)){
            try {
                x_check_again = int_whole_star(h+h_check_again);
                M_check_again = x_check_again[0];
                if (verbose) cout<<"M_check_again: "<<M_check_again<<endl;
                if (cp_func(M_check_again, M)) {jump_back = true; if (cp_func(std::max({M, M_guard}, cp_func), M_temp)) force_jump_back = true;} 
            }
            catch (exception &) { if (verbose) cout<<"Invalid value encountered in check_mmax: h="<<h+0.01<<endl;}
        }
        M_temp = std::max({M, M_guard}, cp_func);
    }
    // take the extreme of M and h
    M = std::max({M, M_guard, M_check_again}, cp_func);
    if (M==M_guard) h = h+h_guard;
    else if (M==M_check_again) h = h+h_check_again;
    else {}// M = M, h = h
    if (verbose) cout<<"find h: "<<h<<",    with M:"<<M<<endl<<endl;
    if (std::isnan(M) or M==0 or h<=local_minimum_allowed_h) {fail_process_check_mmax("check_mmax: find max mass error(got M=(NaN, 0) or h too small)");}
    if (check_ok){
        if (M<minm_tov) {fail_process_check_mmax("check_mmax: find max mass error(M<min_m_tov)");}
        if (check_causal){
            if (verbose) cout<<endl<<"cooling: "<<endl;
            if (cool_eos(&h)){
                try {M = int_whole_star(h)[0];}
                catch (exception&){ fail_process_check_mmax("Invalid value encountered in integrate after cooling");}
                if (M>maxm_tov or M<minm_tov) {fail_process_check_mmax("check_mmax: find max mass error(maximum mass invalid after cooling)");}
                if (verbose) cout<<"Max mass after cooling: "<<M<<endl<<endl;
            }
            else{fail_process_check_mmax("check_mmax: find max mass error(cool_eos not success)");}
        }
        else{
            if (M>maxm_tov) {fail_process_check_mmax("check_mmax: find max mass error without cooling(M>max_m_tov)");}      
        }
    }
    *hc = h; *M_max=M;
    return true;
}

bool check_mmax_gd(double *hc, double *M_max, double h_start, bool cp_func_type_less, double jump_size, bool check_ok){
    auto cp_func = [=](double cp1, double cp2){
        if (cp_func_type_less) return cp1<cp2;
        else return cp1>cp2;
    };
    int step_i = 0, sum_sign = 0;
    bool jump_back=false, force_jump_back=false;
    state_type x_result(3), x_guard(3), x_check_again(3);
    double h_max=eos_table_max_h, h_guard=1e-3, h_check_again=0.01;
    double h = h_start;
    double local_minimum_allowed_h = h_start-1e-10;
    double M=0., M_guard=0., gradient=0.0;
    #define fail_process_check_mmax_gd(message) \
        *hc = -1; *M_max=0.; \
        if (verbose) { \
            cout<<message<<", with h="<<h<<", h_guard="<<h+h_guard; \
            cout<<", M="<<M<<", M_guard="<<M_guard<<endl;\
        } \
        return false;
    //find m_max
    if(verbose) cout<<endl<<"find maximum mass, with h_min: "<<local_minimum_allowed_h<<", h_max: "<<h_max<<endl;
    while(h>local_minimum_allowed_h and h<h_max){
        h += gradient*jump_size;//+rand()%10/1.e4; // avoid dead loop
        try{
            x_result = int_whole_star(h);
            M = x_result[0];
            x_guard = int_whole_star(h+h_guard);
            M_guard = x_guard[0];
            if (cp_func_type_less) gradient = (M_guard-M)/h_guard;
            else gradient = (M-M_guard)/h_guard;
        }
        catch (exception &){
            if (verbose) cout<<"Invalid value encountered in check_mmax, with h="<<h<<", or h="<<h+h_guard<<endl;
            M=0; M_guard=0;
        }
        step_i += 1;
        if (verbose) cout<<"step: "<<step_i<<"    h: "<<h<<"     M:"<<M<<"     M_guard:"<<M_guard<<"     gradient:"<<gradient<<"    jump_size:"<<jump_size<<"     dh:"<<gradient*jump_size<<endl;
        if (std::signbit(gradient)) sum_sign -= 1;
        else sum_sign += 1;
        if (abs(gradient*jump_size)<1e-5 or abs(M-M_guard)<1e-5 or step_i>200) break;
        if (step_i%8==0){
            if (sum_sign==0) jump_size /= 2.;
            else if (sum_sign==8) jump_size *= 5.;
            else {}
            sum_sign = 0.;
        }
    }
    // take the extreme of M and h
    M = std::max({M, M_guard}, cp_func);
    if (M==M_guard) h = h+h_guard;
    else {}// M = M, h = h
    if (verbose) cout<<"find h: "<<h<<",    with M:"<<M<<endl<<endl;
    if (std::isnan(M) or M==0 or h<=local_minimum_allowed_h) {fail_process_check_mmax_gd("check_mmax: find max mass error(got M=(NaN, 0) or h too small)");}
    if (check_ok){
        if (M<minm_tov) {fail_process_check_mmax_gd("check_mmax: find max mass error(M<min_m_tov)");}
        if (check_causal){
            if (verbose) cout<<endl<<"cooling: "<<endl;
            if (cool_eos(&h)){
                try {M = int_whole_star(h)[0];}
                catch (exception&){ fail_process_check_mmax_gd("Invalid value encountered in integrate after cooling");}
                if (M>maxm_tov or M<minm_tov) {fail_process_check_mmax_gd("check_mmax: find max mass error(maximum mass invalid after cooling)");}
                if (verbose) cout<<"Max mass after cooling: "<<M<<endl<<endl;
            }
            else{fail_process_check_mmax_gd("check_mmax: find max mass error(cool_eos not success)");}
        }
        else{
            if (M>maxm_tov) {fail_process_check_mmax_gd("check_mmax: find max mass error without cooling(M>max_m_tov)");}      
        }
    }
    *hc = h; *M_max=M;
    return true;
}

bool check_mmax_pt_two_branch(double* h1, double* h2, double* h3, double* h_max, double* m1, double* m2, double* m3, double* m_max, double rho_tr, double drho){
    double h_tr_start, h_tr_end, jump_size;
    bool ck_tbranch_l=false, ck_tbranch_u=false;
    *h1 = -1, *h2 = -1, *h3 = -1, *h_max = -1, *m1 = 0., *m2 = 0., *m3 = 0., *m_max=0.;
    #define fail_process_check_mmax_pt_two_branch(message) \
        *h1 = -1, *h2 = -1, *h3 = -1, *h_max = -1, *m1 = 0., *m2 = 0., *m3 = 0., *m_max=0.; \
        if (verbose) cout<<message<<endl; \
        return false;
    try {
        h_tr_start = eos_table_function_rho_base[0](rho_tr);
        h_tr_end = eos_table_function_rho_base[0](rho_tr+drho);
        jump_size = (h_tr_end-h_tr_start)/2.;
        if (jump_size<0.01) jump_size = 0.01;
        if (jump_size>0.1) jump_size = 0.1;
    }
    catch (exception&) {h_tr_start=0.; h_tr_end=0.; jump_size=0.01; }
    if (check_mmax_gd(h1, m1, 0.08, true, 0.01,  false)){
        if (verbose) cout<<"h_tr_start="<<h_tr_start<<", h_tr_end="<<h_tr_end<<", h_find="<<*h1<<", delta_h="<<h_tr_start-*h1<<endl;
        if (abs(h_tr_start-*h1)<0.03 or abs(h_tr_end-*h1)<0.03){
            ck_tbranch_l = check_mmax_gd(h2, m2, *h1+0.01, false, 0.02, false);
            if (ck_tbranch_l){
                cout<<"two branches!"<<endl;
                ck_tbranch_u = check_mmax_gd(h3, m3, *h2+0.01, true, 0.2, false);
                if (ck_tbranch_u) {*h_max = *h3; *m_max = *m3;}
                else {fail_process_check_mmax_pt_two_branch("check second branch error");}
            }
        }
        else {
            *h_max = *h1; *m_max = *m1;
        }
        if ((not ck_tbranch_l) or (ck_tbranch_u and *m1>*m3)) {
            *h_max = *h1; *m_max = *m1;
        }
        if (*m_max<minm_tov) {fail_process_check_mmax_pt_two_branch("check_mmax: find max mass error(M<min_m_tov)");}
        if (check_causal){
            if (verbose) cout<<endl<<"cooling: "<<endl;
            if (cool_eos(h_max)){
                try {*m_max = int_whole_star(*h_max)[0];}
                catch (exception&){ fail_process_check_mmax_pt_two_branch("Invalid value encountered in integrate after cooling");}
                if (*m_max>maxm_tov or *m_max<minm_tov) {fail_process_check_mmax_pt_two_branch("check_mmax: find max mass error(maximum mass invalid after cooling)");}
                if (verbose) cout<<"Max mass after cooling: "<<*m_max<<endl<<endl;
            }
            else{fail_process_check_mmax_pt_two_branch("check_mmax: find max mass error(cool_eos not success)");}
        }
        else{
            if (*m_max>maxm_tov) {fail_process_check_mmax_pt_two_branch("check_mmax: find max mass error without cooling(M>max_m_tov)");}      
        }
        return true;
    }
    *h1 = -1, *h2 = -1, *h3 = -1, *h_max = -1, *m1 = 0., *m2 = 0., *m3 = 0., *m_max=0.;
    return false;
}

bool find_closest(double m_aim, double h_i, double h_max, double *h_closest, double *lambda, bool use_user_start_point){
//find nearest mass, return lambda;
    bool jump_back=false, force_jump_back=false;
    state_type x_result(3), x_guard(3), x_check_again(3);
    double h_guard=4e-3;
    double h, h_init, h_save, jump_size = 0.1;
    double M=0., M_guard=0., M_check_again=0., M_temp=0.;
    if (use_user_start_point) h_init = h_i;
    else h_init = suggested_start_point;
    if(verbose){
        cout<<endl<<"find closest aim (known maximum mass): "<<m_aim<<", user specified start point: "<<h_i<<", suggested start point: "<<suggested_start_point;
        cout<<", minimum_allowed_h: "<<minimum_allowed_h<<", use user start point: "<<use_user_start_point<<", h start: "<<h_init<<", h end: "<<h_max<<endl;
    }
    h = h_init, h_save = h_init;
    while(h>minimum_allowed_h and h<h_max and jump_size>=1e-3){
        try{
            x_result = int_whole_star(h);
            M = x_result[0];
            x_guard = int_whole_star(h+h_guard);
            M_guard = x_guard[0];
        }
        catch (exception &){
            cout<<"Invalid value encountered in find nearest mass: h="<<h<<endl;
            M=0; M_guard=0;
        }
        if (verbose) cout<<h<<"    "<<M_guard<<"     "<<M<<endl;
        jump_back = (M==0 or M_guard==0 or std::isnan(M_guard) or std::isnan(M) or max(M, M_guard)>m_aim);
        force_jump_back = false;
        if ((not jump_back) and M_guard<=M){
            try {
                x_check_again = int_whole_star(h+0.01);
                M_check_again = x_check_again[0];
                if (verbose) cout<<"M_check_again: "<<M_check_again<<endl;
                if (M_check_again<M) {jump_back = true; if (max(M, M_guard)<M_temp) force_jump_back=true;} 
            }
            catch (exception &) { if (verbose) cout<<"Invalid value encountered in find maximum mass: h="<<h+0.01<<endl;}
        }
        h_save = h;
        if (jump_back) {
            if (force_jump_back) h -= jump_size;
            h -= jump_size;
            jump_size /= 5.;
            h += jump_size;
        }
        else h += jump_size;
        M_temp = max(M, M_guard);
    }
    if (verbose) cout<<"find h: "<<h_save<<",    with M:"<<M<<endl<<endl;
    if (std::isnan(M) or h_save<=minimum_allowed_h or force_jump_back) {*h_closest = -1; *lambda=0; return false;}
    else {
        *h_closest = h_save;
        if (x_result[0]>0. and x_result[1]>0.) *lambda = cal_lambda(x_result[0]/x_result[1], x_result[2]);
        else *lambda = 0.;
        return true;
    }
}

/**
* @brief Find a global_a with another global_b.
* @param known_aim Global_b.
* @param h_i Initial guess.
* @param h_max Known maximum enthalpy.
* @param h_closest Pointer to the closest enthalpy, return -1 if finding error.
* @param unknown Pointer to the global_a to be found, return 0 if finding error.
* @param get_type 1. mass to lambda; 2. mass to radius; 3. lambda to mass.
* @param use_user_start_point Whether to use user specified initial guess h_i or use system default value.
* @retval bool whether finding process correctly worked.
*/
bool find_closest_with_maxm_known(double known_aim, double h_i, double h_max, double *h_closest, double *unknown, int get_type, bool use_user_start_point){
    int iter = 0;
    bool jump_back = false;
    state_type x_result(3);
    double h_guard=4e-3;
    double h, h_init, h_save, jump_size=0.1, test_known=0.;
    if (use_user_start_point) h_init = h_i;
    else h_init = suggested_start_point;
    #define fail_process_find_global_prop() *h_closest = -1; *unknown=0; return false;
    if(get_type!=1 and get_type!=2 and get_type!=3) {cout<<"not known parameter value of get_type(in func find_closest_with_mmax_known): "<<get_type<<endl; fail_process_find_global_prop();}
    if(verbose){
        cout<<endl<<"find closest aim (known maximum mass): "<<known_aim<<", user specified start point: "<<h_i<<", suggested start point: "<<suggested_start_point;
        cout<<", minimum_allowed_h: "<<minimum_allowed_h<<", use user start point: "<<use_user_start_point<<", h start: "<<h_init<<", h end: "<<h_max<<endl;
    }
    h = h_init, h_save = h_init;
    while(h>minimum_allowed_h and h<h_max and jump_size>=1e-4){
        try{
            x_result = int_whole_star(h);
            if(get_type==1 or get_type==2) test_known = x_result[0];
            else if(get_type==3) test_known = cal_lambda(x_result[0]/x_result[1], x_result[2]);
            else {}
        }
        catch (exception &){
            cout<<"Invalid value encountered in find nearest mass: h="<<h<<endl;
            test_known=0;
        }
        if (verbose) cout<<"step:"<<iter<<"    ,h: "<<h<<"     ,M: "<<test_known<<endl;
        if(get_type==3) jump_back = (test_known==0 or std::isnan(test_known) or test_known<known_aim);
        else jump_back = (test_known==0 or std::isnan(test_known) or test_known>known_aim);
        if (h!=h_init) jump_back |= (h+jump_size>h_max);
        else { if (jump_back and (iter==0)) jump_size = h-minimum_allowed_h;}
        h_save = h;
        if (jump_back) {
            h -= jump_size;
            jump_size /= 5.;
            h += jump_size;
        }
        else h += jump_size;
        iter += 1;
    }
    if (verbose) cout<<"aim: "<<known_aim<<", h: "<<h_save<<", used steps: "<<iter<<", test_known:"<<test_known<<endl;
    if (std::isnan(test_known) or h_save<=minimum_allowed_h) {fail_process_find_global_prop();}
    else {
        *h_closest = h_save;
        if (x_result[0]>0. and x_result[1]>0.) {
            if(get_type==1) {*unknown = cal_lambda(x_result[0]/x_result[1], x_result[2]);}
            else if(get_type==2) {*unknown = x_result[1]*r_trans;}
            else if(get_type==3) {*unknown = x_result[0];}
            else {}
        }
        else *unknown = 0.;
        return true;
    }
}

/**
* @brief Find nearest known properties, return corresponding unknown ones.
* @param known1 Known global 1 to find with.
* @param known2 Known global 2 to find with.
* @param unknown1 Unknown global 1 to be found.
* @param unknown2 Unknown global 2 to be found.
* @param h_max Pointer to stare maximum enthalpy.
* @param get_type 1. mass to lambda; 2. mass to radius; 3. lambda to mass.
* @retval bool Whether finding process correctly worked.
*/
bool get_unknowns_from_knowns(double known1, double known2, double *unknown1, double *unknown2, double *h_max, int get_type){
    double h1=0., h2=0., M_max=0., L_max=0.;
    *unknown1 = 0.; *unknown2 = 0.; *h_max=0.;
    if(not check_mmax(h_max, &M_max)) return false;
    if(get_type==1 or get_type==2){
        if ((known1>M_max) or (known2>M_max)) return false;
    }
    else if (get_type==3){
        get_mrl(*h_max);
        L_max = mrl_result[2];
        if ((known1<L_max) or (known2<L_max)) return false;
    }
    else{
        cout<<"not known parameter value of get_type(in func get_unkowns_from_knowns): "<<get_type<<endl;
        return false;
    }
    if (known1<known2){
        if (find_closest_with_maxm_known(known1, suggested_start_point, *h_max, &h1, unknown1, get_type)) find_closest_with_maxm_known(known2, h1, *h_max, &h2, unknown2, get_type);
        else find_closest_with_maxm_known(known2, suggested_start_point, *h_max, &h2, unknown2, get_type);
    }
    else if(known1>known2){
        if (find_closest_with_maxm_known(known2, suggested_start_point, *h_max, &h2, unknown2, get_type)) find_closest_with_maxm_known(known1, h2, *h_max, &h1, unknown1, get_type);
        else find_closest_with_maxm_known(known1, suggested_start_point, *h_max, &h1, unknown1, get_type);
    }
    else{
        find_closest_with_maxm_known(known1, suggested_start_point, *h_max, &h1, unknown1, get_type);
        *unknown2 = *unknown1;
    }
    if (*unknown1>0. and *unknown2>0.) return true;
    else return false;
}


#endif
