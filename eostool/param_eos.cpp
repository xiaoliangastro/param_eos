#include<iostream>
#include<fstream>
#include<cmath>
#include<exception>
#include<iomanip>
#include<string>
#include<ctime>
#include"optimize_posterior.hpp"
#include"integrate_core.hpp"
#include"param_utils.hpp"
#include"non_param_utils.hpp"
#include"eos_init.hpp"

using namespace std;

//g++ -lgsl -std=c++17 -shared -fPIC -o param_eos.so param_eos.cpp
///opt/rh/devtoolset-8/root/usr/bin/g++ -shared -fPIC -I/Software/boost/include -lgsl -std=c++17 -o param_eos.so param_eos.cpp 

//aim of the program
bool non_param = 0;//non-param method
bool non_param_mktb = 0;//non-param method make mr table
bool routine = 0;//routine(command: rt)
bool checkeos = 0;//check if eos is repeated correctly(command: ckeos)
bool mk_tool_table = 0;//make tool table(command: mktb)
bool opt_file = 0;//(command: optf)
bool opt_sg_point = 0;//optimize single point(command: opts)
bool opt_sg_pair = 0;//optimize single pair(command: optp)

string eos_fname = "/Users/jiangjinliang/work/try/Spectral_EoS/eos_table/std_ebase_lowdense_eos.txt";

//focus of the program
bool command_line_init(int argc, char *argv[]);


int main(int argc, char *argv[])
{

//basic initiate
    cout<<setprecision(12);
    double h_c, h_min, h_max, dh;
    int offset = 0;
    time_t t0 = time(0);
    command_line_init(argc, argv);

//initiate according to the aim of the program
    bool optimizeit = (opt_sg_point or opt_sg_pair) or opt_file;
    if (non_param){
        int step_i = stoi(argv[3]);
        double frac = stod(argv[4]);
        double like = non_param_like(step_i, frac, false);
        const string f_like = "../../result/likelihood_table.txt";
        system(("echo \'"+std::to_string(step_i)+" "+std::to_string(like)+"\' >> "+f_like).c_str());
    }
    else if(non_param_mktb){
        int step_i = stoi(argv[3]);
        double frac = stod(argv[4]);
        double like = non_param_like(step_i, frac, true);
    }
    else if (checkeos or mk_tool_table or routine or optimizeit){
        //initiate eos parameters
        if (piecewise or piecewise_p or spectral or spectral_causal) offset = 4;
        if (quark_star) offset = 3;
        if (interp_only) {offset=1; eos_fname = string(argv[3]);}
        if (piecewise or piecewise_p or spectral or spectral_causal or quark_star){
            double params[4];
            for (int i=0; i<offset; i++) params[i] = stod(argv[i+3]);
            change_pars(params, 4, 0);
        }
        if (piece_spec_phtr_css){
            offset = 8;
            double params[8];
            for (int i=0; i<offset; i++) params[i] = stod(argv[i+3]);
            change_pars(params, 8, 0);
        }
        if (piecewise or piecewise_p){
            cal_eborders(eos_params);
            if (verbose) { cout<<"e-borders: "; for (int i=0; i<e_borders.size(); i++) cout<<e_borders[i]<<"\t"; cout<<endl; }
        }
        //initiate interpolation coefficient or not
        if (not quark_star) init_interp_coeff(eos_fname.c_str());
    }
    else {}

//check if eos is right
    if (checkeos){
        h_min = stod(argv[3+offset]); h_max = stod(argv[4+offset]); dh = stod(argv[5+offset]);
        char fname[] = "../../result/repeat_eos.txt";
        make_eos_table(h_min, h_max, dh, 12, fname, "cgs");
    }

//make tool table for optimization
    if (mk_tool_table){
        h_min = stod(argv[3+offset]); h_max = stod(argv[4+offset]); dh = stod(argv[5+offset]);
        make_tool_table(h_min, h_max, dh);
    }

//routine
    if (routine){
        try {
            h_c = stod(argv[3+offset]);
            state_type x_result = int_whole_star(h_c);
    		double M = x_result[0];
    		double R = x_result[1];
            double C = M/R;
    		double LB = cal_lambda(C, x_result[2]);
    		cout<<"M:"<<M<<"M_sun     R:"<<R*r_trans<<"km"<<"     C:"<<C<<"     Lambda:"<<LB<<endl;        
            //cout<<h_c<<" "<<M<<" "<<R*r_trans<<endl;
        }
        catch (exception &){
            cout<<"Invalid value encountered in routine"<<endl;
            exit(0);
        }

    }

//optimization
    if (optimizeit){
        //read tool table to accelerate finding root
        double val1, val2, val3, val4;
        string tool_fname = "../../tool/eos_table/test_hmrl.txt";
        ifstream fp(tool_fname, ifstream::in);
        while(fp>>val1>>val2>>val3>>val4) toolm_h[val2] = val1;
        fp.close();
        //optimize
        if (opt_sg_point) {
            double mass = stod(argv[3+offset]);
            state_type mrlh = optimize_single_point(mass);
            cout<<"aim m:"<<mass<<"\t\t get m:"<<mrlh[0]<<"\t\t r:"<<mrlh[1]<<"\t\t l:"<<mrlh[2]<<"\t\t h:"<<mrlh[3]<<endl;
        }
        if (opt_sg_pair){
            double m1 = stod(argv[3+offset]), m2 = stod(argv[4+offset]);
            state_type l12t_m12 = optimize_single_pair(m1, m2);
            cout<<"m1:"<<l12t_m12[3]<<"\t\t m2:"<<l12t_m12[4]<<"\t\t l1:"<<l12t_m12[0]<<"\t\t l2:"<<l12t_m12[1]<<"\t\t lt:"<<l12t_m12[2]<<endl;
        }
        if (opt_file){
            //"../../data/ligo/low_ligo_ttest.txt", "../../result/test_ltm12l12.txt"
            string ifname = argv[3+offset];
            string ofname = argv[4+offset];
            //cout<<ifname<<" "<<ofname<<endl;
            optimize_file(ifname, ofname);
        }
    }

    if (verbose) cout<<"used time:"<<(time(0)-t0)/60.0<<"min(s)"<<endl;
}


//----------------------------------------------------------------------------------------------
//initiate functions
//----------------------------------------------------------------------------------------------


bool command_line_init(int argc, char *argv[]){
    string help = "usage: "+string(argv[0])+" aim param_method parameters [-v|-vv|-vvv]\n";
    help += "      (param_method 1:interp_only|2:piecewise|3:piecewise_p|4:spectral|5:spectral_causal|6:phase_trans|7:quark_star) \n";
    help += "      -h, --help   show this massage\n";
    help += "      rt           routine(pars: eos parameters + h_c)\n";
    help += "      ckeos        check if eos is ok(pars: eos parameters + h_min + h_max + dh)\n";
    help += "      mktb         make h-mr table(pars: eos parameters + h_min + h_max + dh)\n";
    help += "      opts         optimize a single mass(pars: eos parameters + m1)\n";
    help += "      optp         optimize a single pair(pars: eos parameters + m1 + m2)\n";
    help += "      optf         optimize a file with pairs of masses(pars: eos parameters + in_fname + out_fname)\n";
    help += "      nonp         non-param method(pars: none)\n";
    help += "      nonpmktb     non-param method make table(pars: none)";
    // aim of the executable
    string state = string(argv[1]);
    if (state=="--help") {cout<<help<<endl; return true;}
    else if (state=="-h") {cout<<help<<endl; return true;}
    else if (state=="rt") routine = true;
    else if (state=="ckeos") checkeos = true;
    else if (state=="mktb") mk_tool_table = true;
    else if (state=="optf") opt_file = true;
    else if (state=="opts") opt_sg_point = true;
    else if (state=="optp") opt_sg_pair = true;
    else if (state=="nonp") {non_param = true; interp_only=true;}
    else if (state=="nonpmktb") {non_param_mktb = true; interp_only=true;}
    else {cout<<"Unknown aim description, exit"<<endl;return false;}
    // parameterization method
    int param_method = stoi(argv[2]);
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
    else {cout<<"Unknown parameterization method: "<<param_method<<endl; return false;}
    for (int i=0; i<argc; i++){
        if (string(argv[i])=="-v") verbose = true;
        if (string(argv[i])=="-vv") {verbose = true; vverbose = true;}
        if (string(argv[i])=="-vvv") {vvverbose = true; vverbose = true; verbose = true;}
    }
    return true;
}

//debuger

//cout<<"let's go"<<endl;
//cout<<"I am fucked!"<<endl;
//cout<<"bye"<<endl;
// cout<<"I am pretty fucking far from ok!"<<endl;
// cout<<"I am far from ok!"<<endl;
// cout<<"I am pretty far from ok!"<<endl;