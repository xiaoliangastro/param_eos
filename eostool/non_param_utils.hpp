#ifndef NON_PARAM_UTILS_HPP
#define NON_PARAM_UTILS_HPP


#include"param_utils.hpp"
#include"integrate_core.hpp"
#include<boost/math/interpolators/makima.hpp>

using boost::math::interpolators::makima;

void transform_eostable(state_type e, state_type p){
    int len_eos = e.size();
    state_type cc, hh;
    double key = coeftb.end()->first;
    hh.push_back(key);
    for(int i=0; i<len_eos-1; i++){
        cc.push_back(log(p[i+1]/p[i])/log(e[i+1]/e[i]));
        hh.push_back(hh[i]+(cc[i]/(cc[i]-1))*log((e[i]*(e[i+1]+p[i+1]))/(e[i+1]*(e[i]+p[i]))));
    }
    cc.push_back(cc[len_eos-2]);
    for(int i=0; i<len_eos; i++){ 
         state_type add{p[i], e[i], cc[i]};
         coeftb[hh[i]] = add;
    }
    //for(auto it=coeftb.begin(); it!=coeftb.end(); it++){
    //    cout<<"coef table: "<<it->first<<":";
    //    for(int itt=0; itt<it->second.size(); itt++) cout<<"  "<<it->second[itt];
    //    cout<<endl;
    //}
}

double non_param_like(int step_i, double frac, bool mktb){
    auto check_nonparam_causal = [=] (double e1, double p1, double e2, double p2){ return (p2-p1)/(e2-e1)>0. and (p2-p1)/(e2-e1)<1; };
    state_type e_upd, p_upd, m_org, r_org, m_c, r_c;
    double likelihood=0, like_product=1, safe_hmax_error=0.05,  min_m=1.0, min_h=0.05;
    double m_max, hc_max, h_step, m_step, r_stand, r_comp, safe_hcm, safe_mm, val1, val2;
    int n_table=20, table_rows, len_eos, si;
    bool iscausal, dummy;
    //initiate tables and eos
    table_rows = sizeof(MR_TABLE)/sizeof(MR_TABLE[0][0])/2;
    for(int i=0; i<table_rows; i++) {m_org.push_back(MR_TABLE[i][0]); r_org.push_back(MR_TABLE[i][1]);}
    table_rows = sizeof(EP_TABLE)/sizeof(EP_TABLE[0][0])/2;
    for(int i=0; i<table_rows; i++) {e_upd.push_back(EP_TABLE[i][0]/rho_trans); p_upd.push_back(EP_TABLE[i][1]/p_trans);}
    table_rows = sizeof(COEFF_TABLE)/sizeof(COEFF_TABLE[0][0])/4;
    for(int i=0; i<table_rows; i++) {coeftb[COEFF_TABLE[i][3]] = state_type({COEFF_TABLE[i][0], COEFF_TABLE[i][1], COEFF_TABLE[i][2]});}
    len_eos = e_upd.size();
    if(mktb){
        ifstream fp_chooselist("../result/choose_list.txt", ifstream::in);
        while(fp_chooselist>>val1>>val2) {
            if(val2>0) p_upd[val2] *= 1.+frac;
            else if(val2<0) p_upd[abs(val2)] *= 1.-frac;
            else {}
            if(val1==step_i) break;
        }
        fp_chooselist.close();
    }
    else{
        if(step_i==0) {
            iscausal = true;
        }
        else{
            if(step_i<=len_eos) {si = step_i-1; p_upd[si] *= (1.+frac);}
            else {si = step_i-len_eos-1; p_upd[si] /= (1.-frac);}
            //judge causal
            if(si==0) iscausal = check_nonparam_causal(e_upd[si], p_upd[si], e_upd[si+1], p_upd[si+1]);
            else if(si==(len_eos-1)) iscausal = check_nonparam_causal(e_upd[si-1], p_upd[si-1], e_upd[si], p_upd[si]);
            else iscausal = check_nonparam_causal(e_upd[si], p_upd[si], e_upd[si+1], p_upd[si+1]) and \
                            check_nonparam_causal(e_upd[si-1], p_upd[si-1], e_upd[si], p_upd[si]);
        }
    }
    if(iscausal){
        transform_eostable(e_upd, p_upd);
        dummy = check_mmax(&hc_max, &m_max);
        if (dummy){
            safe_hcm = hc_max-safe_hmax_error; //safe_mm = get_mrl(safe_hcm)[0];
            h_step = (safe_hcm-min_h)/double(n_table);
            m_step = 0.03;//(std::min(safe_mm, *std::max_element(m_org.begin(), m_org.end()))-min_m)/double(n_table);
            for(int i=0; i<n_table; i++){
                get_mrl(min_h+h_step*i);
                m_c.push_back(mrl_result[0]);
                r_c.push_back(mrl_result[1]);
                //cout<<min_h+h_step*i<<"   "<<m_c[i]<<"  "<<r_c[i]<<endl;
            }
            double max_morg = *std::max_element(m_org.begin(), m_org.end());
            double max_mc = *std::max_element(m_c.begin(), m_c.end());
            double max_m = std::min(max_morg, max_mc);
            auto f_compare = makima(std::move(m_c), std::move(r_c));
            auto f_standard = makima(std::move(m_org), std::move(r_org));
            int step_m = 0;
            double mm = min_m;
            if(mktb) system(":> ../result/cpmr_table.txt");
            while(mm<max_m){
                r_stand = f_standard(mm);
                r_comp = f_compare(mm);
                likelihood += exp(-pow(r_stand-r_comp, 2)/r_stand);
                if(mktb) system(("echo \'"+std::to_string(mm)+"   "+std::to_string(r_comp)+"   "+std::to_string(r_stand)+"\' >> ../result/cpmr_table.txt").c_str()); 
                step_m += 1;
                mm = min_m+step_m*m_step;
            }
            //likelihood = like_product;
        }
        else{ if(verbose) cout<<"maximum mass not allowed in step_i:"<<step_i<<endl; }
    }
    else{ if(verbose) cout<<"non causal encountered in step_i:"<<step_i<<endl; }
    //cout<<"step: "<<step_i<<"  "<<likelihood<<endl;
    return likelihood;
}


#endif
