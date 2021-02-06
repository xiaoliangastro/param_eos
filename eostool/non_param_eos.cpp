#include<iostream>
#include<fstream>
#include<string>
#include<algorithm>
#include<cstdlib>
#include<pthread.h>
#include"param_utils.hpp"
#include"non_param_utils.hpp"

///opt/rh/devtoolset-8/root/usr/bin/g++ -I/Software/boost/include -lgsl -lpthread  -std=c++17 -o run non_param_run.cpp
using namespace std;


struct non_param_struct{
    int step_i;
    double frac;
};

inline bool my_compare(map<int, double>::value_type &it1, map<int, double>::value_type &it2){return it1.second<it2.second;}

void* cal_like(void *args){
    struct non_param_struct *data = static_cast<non_param_struct*>(args);
    const string exe_str = "./do nonp "+std::to_string(data->step_i)+" "+std::to_string(data->frac);
    system(exe_str.c_str());
    pthread_exit(NULL);
}

void create_data(const string f_data, const string f_mr, const string f_coeff, state_type e_org, state_type p_org){
    int table_rows;
    double val0, val1, val2, val3;
    state_type m_org, r_org, p_tb, e_tb, c_tb, h_tb;
    ofstream fp_data(f_data, ofstream::out);
    //set coeff table
    ifstream fp_coeff(f_coeff, ifstream::in);
    while(fp_coeff>>val0>>val1>>val2>>val3){
        p_tb.push_back(val0); e_tb.push_back(val1); 
        c_tb.push_back(val2); h_tb.push_back(val3); 
    }
    fp_coeff.close();
    table_rows = e_tb.size();
    fp_data<<"#ifndef DATA_HPP"<<endl<<"#define DATA_HPP"<<endl<<endl<<endl;
    fp_data<<"const double COEFF_TABLE["<<table_rows<<"][4] = {{"<<p_tb[0]<<", "<<e_tb[0]<<", "<<c_tb[0]<<", "<<h_tb[0]<<"}, \\"<<endl;
    for(int i=1; i<table_rows-1; i++) fp_data<<string(37, ' ')<<'{'<<p_tb[i]<<", "<<e_tb[i]<<", "<<c_tb[i]<<", "<<h_tb[i]<<"}, \\"<<endl;
    fp_data<<string(37, ' ')<<'{'<<p_tb[table_rows-1]<<", "<<e_tb[table_rows-1]<<", "<<c_tb[table_rows-1]<<", "<<h_tb[table_rows-1]<<"}};"<<endl<<endl;
    //set m-r table
    ifstream fp_mr(f_mr, ifstream::in);
    fp_mr.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while(fp_mr>>val0>>val1>>val2>>val3){ m_org.push_back(val1); r_org.push_back(val2); }
    fp_mr.close();
    table_rows = m_org.size();
    fp_data<<"const double MR_TABLE["<<table_rows<<"][2] = {{"<<m_org[0]<<", "<<r_org[0]<<"}, \\"<<endl;
    for(int i=1; i<table_rows-1; i++) fp_data<<string(33, ' ')<<'{'<<m_org[i]<<", "<<r_org[i]<<"}, \\"<<endl;
    fp_data<<string(33, ' ')<<'{'<<m_org[table_rows-1]<<", "<<r_org[table_rows-1]<<"}};"<<endl<<endl;
    //set eos table
    table_rows = e_org.size();
    fp_data<<"const double EP_TABLE["<<table_rows<<"][2] = {{"<<e_org[0]<<", "<<p_org[0]<<"}, \\"<<endl;
    for(int i=1; i<table_rows-1; i++) fp_data<<string(32, ' ')<<'{'<<e_org[i]<<", "<<p_org[i]<<"}, \\"<<endl;
    fp_data<<string(32, ' ')<<'{'<<e_org[table_rows-1]<<", "<<p_org[table_rows-1]<<"}};"<<endl<<endl;
    fp_data<<endl<<"#endif"<<endl;
    fp_data.close();
}


int main(int argc, char *argv[]){
    int steps = 2;
    double frac = 0.01, val0, val1, val2, val3, like_base;
    state_type e_org, p_org;
    map<int, double> likelihood_table;
    //files to be used
    const string f_base1 = "../";
    const string f_base2 = "../";
    const string f_mr = f_base1+"data/ap3_hmrl.txt";
    const string f_coeff = f_base1+"eos_table/std_ebase_lowdense_eos.txt";
    const string f_data = f_base2+"code/data.h";
    const string fi_eos = f_base2+"result/res_eos_bk.txt";
    const string fo_eos = f_base2+"result/res_eos.txt";
    const string f_like = f_base2+"result/likelihood_table.txt";
    //load eos table
    ifstream fpi_eos(fi_eos, ifstream::in);
    while(fpi_eos>>val1>>val2){
        e_org.push_back(val1);
        p_org.push_back(val2);
    }
    fpi_eos.close();
    int len_eos = e_org.size();
    state_type e_upd(e_org), p_upd(p_org);
    //create threads and initiate
    const int num_threads = 2*len_eos+1;
    pthread_t threads[num_threads];
    pthread_attr_t attr;
    void *status;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    non_param_struct args[num_threads];
    for(int i=0; i<num_threads; i++) {args[i].step_i = i; args[i].frac = frac;}
    //run
    for(int i=0; i<steps; i++){
        cout<<endl<<"============================"<<endl<<"=========step "<<i+1<<"========="<<endl<<"============================"<<endl;
        create_data(f_data, f_mr, f_coeff, e_upd, p_upd);
        likelihood_table.clear();
        system((":> "+f_like).c_str());
        system(("/opt/rh/devtoolset-8/root/usr/bin/g++ -I/Software/boost/include -lgsl -std=c++17 -o do "+f_base1+"code/rml_rkf78.cpp").c_str());
        for(int j=0; j<num_threads; j++){
            int rc = pthread_create(&threads[j], &attr, cal_like, &args[j]);
            if(rc){
                cout<<"Error: unable to create thread "<<j<<", "<<rc<<endl;
                exit(-1);
            }
        }
         // -shared -fPIC  -lgsl -std=c++17 -o rkf78_eos.so rml_rkf78.cpp 
        for(int j=0; j<num_threads; j++){
           int rc = pthread_join(threads[j], &status);
           if(rc){
              cout<<"Error: unable to join thread "<<j<<", "<<rc<<endl;
              exit(-1);
           }
        }
        ifstream fp_like(f_like, ifstream::in);
        while(fp_like>>val0>>val1) {likelihood_table[val0] = val1;}
        fp_like.close();
        //find max change and adjust
        like_base = likelihood_table[0];
        int idx = std::distance(likelihood_table.begin(), std::max_element(likelihood_table.begin(), likelihood_table.end(), my_compare));
        if(idx==0) {cout<<"No better choice, better remain unchanged at step "<<i<<endl; break;}
        else if(idx<=len_eos) p_upd[idx-1] *= (1.+frac);
        else p_upd[idx-len_eos-1] /= (1.-frac);
        cout<<"chosen idx: "<<pow(-1, idx/len_eos)*(idx%len_eos)<<endl;
    }
    //output eos
    ofstream fpo_eos(fo_eos, ofstream::out);
    for(int i=0; i<len_eos; i++){fpo_eos<<e_upd[i]<<"     "<<p_upd[i]<<endl;}
    fpo_eos.close();
}
