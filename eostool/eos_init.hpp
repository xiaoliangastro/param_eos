#ifndef EOS_INIT_HPP
#define EOS_INIT_HPP

#include<string>
#include"global_variable_constants.hpp"
#include"quark_star_utils.hpp"


struct control_params{
    double min_tov_mass, max_tov_mass;
    bool check_causal; int param_method, verbose_level;
    double eos_table_max_h; char init_interpolation_function_type[3];
    double const_inter_step, dev_G;
    bool consid_const_inter_step, consid_dev_G;
};

extern "C"{
    void init_interp_coeff(const char *e_name);
    void init_interp_coeff_through_ep_array(double ee[], double pp[], const int len_eos);
    void init_control_params(control_params parameters);
    bool change_pars(double ppar[], int len_par, int len_border, double h_max=eos_table_max_h, const char init_function_type[]=init_interpolation_function_type);
}


void init_interp_coeff(const char *eosf_name){
//read table to get interpolation coefficients
    string ifname = string(eosf_name);
    if (verbose) cout<<"Load coeff file: "<<ifname<<endl;
    //string ifname = "/Users/jiangjinliang/work/try/Spectral_EoS/data/std_lowdense_coeff.txt";
    ifstream fp(ifname, ifstream::in);
    double val1, val2, val3, val4;
    state_type coeff(3);
    coeftb.clear();
    while(fp>>val1>>val2>>val3>>val4){
        coeff[0] = val1; coeff[1] = val2;
        coeff[2] = val3; coeftb[val4] = coeff;
    }
    fp.close();
    //map<double, state_type>::iterator id;
    //for (id = coeftb.begin(); id!= coeftb.end(); id++) cout<<"indexes:"<<id->second[0]<<"\t\t"<<id->second[1]<<"\t\t"<<id->second[2]<<"\t\t"<<id->first<<endl;
}

void init_interp_coeff_through_ep_array(double ee[], double pp[], const int len_eos){
//transform ep table to get interpolation coefficients
    state_type cc(len_eos), hh(len_eos);
    coeftb.clear();
    coeftb[0.0] = state_type{pp[0], ee[0], 5.0/3.0};
    hh[0] = (5.0/2.0)*log((ee[0]+pp[0])/ee[0]);
    for (int i=0; i<len_eos-1; i++){
        cc[i] = log(pp[i+1]/pp[i])/log(ee[i+1]/ee[i]);
        hh[i+1] = hh[i]+(cc[i]/(cc[i]-1))*log((ee[i]*(ee[i+1]+pp[i+1]))/(ee[i+1]*(ee[i]+pp[i])));
        coeftb[hh[i]] = state_type{pp[i], ee[i], cc[i]};
    }
    cc[len_eos-1] = cc[len_eos-2];
    coeftb[hh[len_eos-1]] = state_type{pp[len_eos-1], ee[len_eos-1], cc[len_eos-1]};
    //map<double, state_type>::iterator id;
    //for (id = coeftb.begin(); id!= coeftb.end(); id++) cout<<"indexes:"<<id->second[0]<<"\t\t"<<id->second[1]<<"\t\t"<<id->second[2]<<"\t\t"<<id->first<<endl;
}

void init_G_related_parameters(){
    deviation_G_plus1 = deviation_G+1.;
    G = GSL_CONST_CGS_GRAVITATIONAL_CONSTANT*deviation_G_plus1;
    gEicm3_to_dmls = pow(G, 3)*pow(Ms, 2)/pow(C, 6);
    ergEicm3_to_dmls = pow(G, 3)*pow(Ms, 2)/pow(C, 8);
    MeV4_to_dmls = MeV3_to_ifm3*MeVEifm3_to_ergEicm3*ergEicm3_to_dmls;
    length_trans = G*Ms/pow(C, 2);
    rho_trans = pow(C, 6)/(pow(G, 3)*pow(Ms, 2));
    p_trans = rho_trans*pow(C, 2);
    r_trans = length_trans/1.e5;
    e_0 = 9.075720079516068750e+13/rho_trans;
    p_0 = 2.979148040306152863e+32/p_trans;
    rho_0 = (e_0+p_0)/exp(h_0);
    rho_sat = 2.7e+14/rho_trans;
    max_possible_e = 1000*rho_sat;
    rho_borders[0] = 1.*rho_sat; rho_borders[1] = 1.85*rho_sat;
    rho_borders[2] = 3.7*rho_sat; rho_borders[3] = 7.4*rho_sat;
}

void init_control_params(control_params pars){
    int param_method = pars.param_method, verbose_level = pars.verbose_level;// parameterization method and verbose level
    minm_tov=pars.min_tov_mass, maxm_tov = pars.max_tov_mass, check_causal=pars.check_causal;// maximum mass constraint
    eos_table_max_h = pars.eos_table_max_h; for(int i=0; i<3; i++) init_interpolation_function_type[i] = pars.init_interpolation_function_type[i];// init cal_eos_table related parameters
    sg_const_step = pars.const_inter_step; consid_const_inter_step = pars.consid_const_inter_step;// constant integrate step
    //cout<<"TOV: "<<minm_tov<<"  "<<maxm_tov<<" "<<check_causal<<" "<<param_method<<" "<<verbose_level<<" "<<eos_table_max_h<<" "<<sg_const_step<<" "<<pars.dev_G<<" "<<consid_const_inter_step<<" "<<pars.consid_dev_G<<endl;
    if (pars.consid_dev_G) {deviation_G = pars.dev_G; init_G_related_parameters();}// deviation of gravitational constant
    //initiate
    interp_only = 0; piecewise = 0; piecewise_p = 0; spectral = 0; spectral_causal = 0;
    phase_trans = 0; quark_star = 0; adapt_piecewise = 0; cons_cs = 0; piece_spec_phtr_css = 0;
    //parameterization method
    if (param_method==1) {interp_only=1;}
    else if (param_method==2) {piecewise=1;}
    else if (param_method==3) {piecewise_p=1;}
    else if (param_method==4) {spectral=1;}
    else if (param_method==5) {spectral_causal=1;}
    else if (param_method==6) {phase_trans=1;}
    else if (param_method==7) {quark_star=1;}
    else if (param_method==8) {adapt_piecewise=1;}
    else if (param_method==9) {cons_cs=1;}
    else if (param_method==10) {piece_spec_phtr_css=1;}
    else {cout<<"Unknown parameterization method: "<<param_method<<endl; exit(1);}
    // verbose level
    if (verbose_level==0) {}
    else if (verbose_level==1) {verbose=1; vverbose=0; vvverbose=0;}
    else if (verbose_level==2) {verbose=1; vverbose=1; vvverbose=0;}
    else if (verbose_level==3){verbose=1; vverbose=1; vvverbose=1;}
    else {cout<<"Unknown verbose level: "<<verbose_level<<endl; exit(1);}
}


//----------------------------------------------------------------------------------------------
//parameter utilities
//----------------------------------------------------------------------------------------------



void cal_eborders(state_type params){
    e_borders.clear();
    state_type gamma, pb, rb;
    if (piecewise_p or adapt_piecewise or piece_spec_phtr_css){
        for(int i=0; i<params.size(); i++){
            if (i==0) gamma.push_back(log(params[0]/p_0)/log(rho_borders[0]/rho_0));
            else gamma.push_back(log(params[i]/params[i-1])/log(rho_borders[i]/rho_borders[i-1]));
        }
    }
    else {
        for(int i=0; i<params.size(); i++) gamma.push_back(params[i]);
    }
    int n_pieces = gamma.size();
    double ai_plus1;
    e_borders.push_back(e_0);
    rb.push_back(rho_0);
    pb.push_back(p_0);
    for (int i=0; i<n_pieces; i++){
        rb.push_back(rho_borders[i]);
        pb.push_back(pb[i]*pow(rb[i+1]/rb[i],gamma[i]));
        if (gamma[i]==1){
            ai_plus1 = e_borders[i]/rb[i]-pb[i]*log(rb[i])/rb[i];
            e_borders.push_back(rb[i+1]*ai_plus1+pb[i+1]*log(rb[i+1]));
        }
        else {
            ai_plus1 = e_borders[i]/rb[i]-pb[i]/(rb[i]*(gamma[i]-1));
            e_borders.push_back(rb[i+1]*ai_plus1+pb[i+1]/(gamma[i]-1));
        }
    }
    if (phase_trans) {
        e_borders.push_back(rho_borders[n_pieces]);
        e_borders.push_back(max_possible_e);
    }
    if(verbose){
        cout<<"gamma: ";
        for(int i=0; i<gamma.size(); i++) cout<<gamma[i]<<"  ";
        cout<<endl;
        cout<<"e_borders: ";
        for(int i=0; i<e_borders.size(); i++) cout<<e_borders[i]<<"  ";
        cout<<endl;
    }
}

void cal_eos_table(double h_max, const char init_function_type[]){
    state_type eos_tb{p_0, e_0, rho_0};
    vector<state_type> x_o;
    eos_table_h.clear(); eos_table_p.clear(); eos_table_e.clear(); eos_table_rho.clear();
    try {
        integrate_adaptive(controlled_stepper_cal_eos, cal_eos, eos_tb, h_0, h_max, sg_step, push_back_state_and_time(x_o, eos_table_h));
        for (int i=0; i<eos_table_h.size(); i++){ eos_table_p.push_back(x_o[i][0]), eos_table_e.push_back(x_o[i][1]), eos_table_rho.push_back(x_o[i][2]); }
        if (init_function_type[0]=='1'){
            eos_table_function_h_base.clear();
            state_type h_bk1(eos_table_h), h_bk2(eos_table_h), h_bk3(eos_table_h);
            state_type p_bk(eos_table_p), e_bk(eos_table_e), rho_bk(eos_table_rho);
            auto function_h_p = pchip(std::move(h_bk1), std::move(p_bk));
            auto function_h_e = pchip(std::move(h_bk2), std::move(e_bk));
            auto function_h_rho = pchip(std::move(h_bk3), std::move(rho_bk));
            eos_table_function_h_base.push_back(function_h_p);
            eos_table_function_h_base.push_back(function_h_e);
            eos_table_function_h_base.push_back(function_h_rho);
        }
        if (init_function_type[1]=='1'){
            eos_table_function_e_base.clear();
            state_type e_bk1(eos_table_e), e_bk2(eos_table_e), e_bk3(eos_table_e);
            state_type h_bk(eos_table_h), p_bk(eos_table_p), rho_bk(eos_table_rho);
            auto function_e_h = pchip(std::move(e_bk1), std::move(h_bk));
            auto function_e_p = pchip(std::move(e_bk2), std::move(p_bk));
            auto function_e_rho = pchip(std::move(e_bk3), std::move(rho_bk));
            eos_table_function_e_base.push_back(function_e_h);
            eos_table_function_e_base.push_back(function_e_p);
            eos_table_function_e_base.push_back(function_e_rho);
        }
        if (init_function_type[2]=='1'){
            eos_table_function_rho_base.clear();
            state_type rho_bk1(eos_table_rho), rho_bk2(eos_table_rho), rho_bk3(eos_table_rho);
            state_type h_bk(eos_table_h), p_bk(eos_table_p), e_bk(eos_table_e);
            auto function_rho_h = pchip(std::move(rho_bk1), std::move(h_bk));
            auto function_rho_p = pchip(std::move(rho_bk2), std::move(p_bk));
            auto function_rho_e = pchip(std::move(rho_bk3), std::move(e_bk));
            eos_table_function_rho_base.push_back(function_rho_h);
            eos_table_function_rho_base.push_back(function_rho_p);
            eos_table_function_rho_base.push_back(function_rho_e);
        }
        if (verbose) {cout<<"init_interpolation_function_type: "<<string(init_function_type)<<endl;}
    }
    catch (exception & except) {
        cout<<"error encountered in cal_eos_table: "<<except.what()<<endl;
        throw;
    }
}

inline bool check_p_at_rho185(){ return eos_table_function_rho_base[1](1.85*rho_sat)*p_trans>1.21e34;}

bool check_gamma_piece_spec_phtr_css(double start_h=h_0, double max_check_rho=5.){
    double h = start_h;
    double rho, p, e, gamma;
    double p_start = eos_table_function_rho_base[1](1.*rho_sat), e_start = eos_table_function_rho_base[2](1.*rho_sat); 
    double gamma_start = (e_start+p_start)/(p_start*(1.+exp(eos_params[4])));
    bool valid = ((gamma_start>1.4) and (gamma_start<10.));
    if (verbose and (not valid)) cout<<"error in change_pars: gamma is not allowed at the 1rho_sat point!"<<endl;
    while(valid){
        p = eos_table_function_h_base[0](h);
        e = eos_table_function_h_base[1](h);
        rho = eos_table_function_h_base[2](h);
        if (rho>=max_check_rho*rho_sat) break;
        else if (rho<=rho_sat) {}
        else {
            gamma = gamma_which(h, p, e, eos_params);
            valid &= ((gamma>1.4) and (gamma<10.));
        }
        h += 0.01;
    }
    return valid;
}

bool change_pars(double ppar[], int len_par, int extra_par, double h_max, const char init_function_type[]){
    eos_params.clear();
    if (piecewise or piecewise_p) {for (int i=0; i<len_par; i++) eos_params.push_back(ppar[i]); cal_eborders(eos_params);}
    else if (cons_cs) { eos_params.push_back(ppar[0]); }//transfered ppar: cs
    else if (adapt_piecewise) {
        //transfered ppar: rho_borders_ratios*N+p_ratios*N
        int len_border = extra_par;
        rho_borders.clear();
        for (int i=0; i<len_border; i++){
            if (i==0){
                rho_borders.push_back(rho_0*ppar[i]);
                eos_params.push_back(p_0*ppar[len_border+i]);
            }
            else{
                rho_borders.push_back(rho_borders[i-1]*ppar[i]);
                eos_params.push_back(eos_params[i-1]*ppar[len_border+i]);
            }
        }
        cal_eborders(eos_params);
    }
    else if (piece_spec_phtr_css){
        //transfered ppar: p1, g0, g1, g2, rho_tr, deltar_tr, gamma_pt, cs
        //expected eos_params: rhob1, gamma1, rhob2, h0, g0, g1, g2, rhob3, gamma_pt, cs
        rho_borders.clear();
        double rhob_1 = rho_sat, pb_1 = ppar[0];
        double eb_1, h0, gamma1;
        gamma1 = log(pb_1/p_0)/log(rhob_1/rho_0);
        eos_params.push_back(rhob_1); eos_params.push_back(gamma1); eos_params.push_back(ppar[4]);
        state_type p_borders{pb_1};
        rho_borders.push_back(rhob_1);
        cal_eborders(p_borders);
        eb_1 = e_borders[1];
        h0 = log((eb_1+pb_1)/rhob_1);
        eos_params.push_back(h0);
        eos_params.push_back(ppar[1]); eos_params.push_back(ppar[2]); eos_params.push_back(ppar[3]);
        eos_params.push_back(ppar[4]+ppar[5]); eos_params.push_back(ppar[6]); eos_params.push_back(ppar[7]);
        if (verbose){
            cout<<"eos param:"<<endl<<endl;
            for (int j=0; j<eos_params.size(); j++) cout<<eos_params[j]<<"  ";
            cout<<endl<<endl;
        }
    }
    else if (phase_trans) {
        int len_border = extra_par;
        rho_borders.clear();
        for (int i=0; i<len_border; i++) rho_borders.push_back(ppar[i]);
        for (int i=len_border; i<len_par; i++) eos_params.push_back(ppar[i]);
        cal_eborders(eos_params);
    }
    else {for (int i=0; i<len_par; i++) eos_params.push_back(ppar[i]);}
    try{
        if (quark_star) {
            eos_params.push_back(solve_surface_h(-1.0, 0.7, eos_params));
            minimum_allowed_h = eos_params[3];
            suggested_start_point = minimum_allowed_h+0.1;
            if (verbose) cout<<"surf_h: "<<eos_params[3]<<", suggested h: "<<suggested_start_point<<endl;
        }
        if (not (quark_star or interp_only)){
            eos_table_max_h = h_max;
            cal_eos_table(h_max, init_function_type);
        }
        if (piece_spec_phtr_css){
            if (extra_par==1){
                if (not check_p_at_rho185()){
                    if (verbose) cout<<"error in change_pars: p="<<eos_table_function_rho_base[1](1.85*rho_sat)*p_trans<<" at rho=1.85*rho_sat not allowed!"<<endl;
                    return false;
                }
            }
            if (extra_par==2){
                double start_h = eos_table_function_rho_base[0]((1.+1e-10)*rho_sat);
                if (not check_gamma_piece_spec_phtr_css(start_h, 5.)){
                    if (verbose) cout<<"error in change_pars: gamma is not allowed in the 1-5 rho_sat range!"<<endl;
                    return false;
                }
            }
        }
    }
    catch (exception &except){
        cout<<"failed to change parameters in change_pars, "<<except.what()<<endl;
        return false;
    }
    return true;
}


#endif