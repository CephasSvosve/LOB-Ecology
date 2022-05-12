


//
// Created by Cephas Svosve on 5/1/2022.
//

#ifndef UNTITLED32_WIENER_H
#define UNTITLED32_WIENER_H

#include "stochastic_math.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>

using namespace Eigen;
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enum processes{Earnings, Dividends, Free_Cash_Flows};

static tuple<VectorXd,VectorXd,VectorXd,MatrixXd, map<int,vector<double>>, double> params(processes process){
//trade period
    double trade_period = 10000;
    double int_rate = 0.02;
    double int_rate_12 = pow((1+int_rate),1./12) - 1;



//cross-sectional correlations
    MatrixXd div_cross_corr(3, 3);
    div_cross_corr <<   1, 0.129529841, 0.035547681,
                        0.129529841, 1, -0.079728063,
                        0.035547681,  -0.079728063,1;


    MatrixXd earnings_cross_corr(3, 3);
    earnings_cross_corr << 1, 0.251, 0.034,
            0.251, 1, 0.3,
            0.034,  0.3,1;

    MatrixXd fcf_cross_corr(3, 3);
    fcf_cross_corr << 1, 0.251, 0.034,
            0.251, 1, 0.3,
            0.034,  0.3,1;


//auto-correlations
    VectorXd div_aut = VectorXd(3);
    div_aut(0) =        0.819283923;
    div_aut(1) =        -0.406900449;
    div_aut(2) =        0.674315309;



    VectorXd earnings_aut = VectorXd(3);
    earnings_aut(0) =   0.02;
    earnings_aut(1) =   0.13;
    earnings_aut(2) =   0.36;

    VectorXd fcf_aut = VectorXd(3);
    fcf_aut(0) =        0.02;
    fcf_aut(1) =        0.13;
    fcf_aut(2) =        0.36;

//mu

//equivalent daily div growth

    double _g1 = pow((1+0.093472171),1./252)-1;
    double _g2 = pow((1+0.088976849),1./252)-1;
    double _g3 = pow((1+0.061685269),1./252)-1;


    VectorXd div_mu = VectorXd(3);
    div_mu(0) = _g1;//(1./252)* 0.021470588;
    div_mu(1) = _g2;//(1./252)* 0.011764706;
    div_mu(2) = _g3;//(1./252)* 0.002352941;




    VectorXd earnings_mu = VectorXd(3);
    earnings_mu(0) = (1./63) *0.031949098;
    earnings_mu(1) =  (1./63) *0.031869185;
    earnings_mu(2) =  (1./63) *0.028311525;

    VectorXd fcf_mu = VectorXd(3);
    fcf_mu(0) = 0.01;
    fcf_mu(1) = 0.02;
    fcf_mu(2) = 0.03;

//sigma
    VectorXd div_sig = VectorXd(3);
    div_sig(0) = sqrt(1./252)* 0.042191057;
    div_sig(1) = sqrt(1./252)* 0.175544468;
    div_sig(2) = sqrt(1./252)* 0.034144154;


    VectorXd earnings_sig = VectorXd(3);
    earnings_sig(0) = (1./63) * 0.121841232;
    earnings_sig(1) = (1./63) * 0.1921385;
    earnings_sig(2) = (1./63) * 0.307994611;

    VectorXd fcf_sig = VectorXd(3);
    fcf_sig(0) = 001.2;
    fcf_sig(1) = 002.3;
    fcf_sig(2) = 007.6;

    //Initial values

    double annuity_factor_1 = (1-pow((1+_g1)/(1+int_rate_12),252))/(int_rate_12-_g1);
    double annuity_factor_2 = (1-pow((1+_g2)/(1+int_rate_12),252))/(int_rate_12-_g2);
    double annuity_factor_3 = (1-pow((1+_g3)/(1+int_rate_12),252))/(int_rate_12-_g3);


    //equivalent daily payment annuity factor
    double _a1 = (pow((1+int_rate),-1.))/annuity_factor_1;
    double _a2 = (pow((1+int_rate),-1.))/annuity_factor_2;
    double _a3 = (pow((1+int_rate),-1.))/annuity_factor_3;


    vector<double> E_o = {10.23,9.16, 7.23}; //TODO add initial values for dividends and fcf
    vector<double> D_o = {0.07 * _a1,	0.0869 * _a2,	0.2775 * _a3};
    vector<double> D_o_1 = {0.1075,	0.32,	0.0833};
    vector<double> FCF_o = {45.39, 42.80, 24.72};

    //initial dividend growth rate
    vector<double> g_o = {0.23566879,	0,	0.183783784};

    map<int,vector<double>> Earnings_;
    map<int,vector<double>> Div;
    map<int,vector<double>> CashFlow;

    Earnings_.emplace(1,E_o);

    Div.emplace(0,g_o);
    Div.emplace(1,D_o);
    Div.emplace(2,D_o_1);

    CashFlow.emplace(1,FCF_o);



    tuple<VectorXd,VectorXd,VectorXd,MatrixXd, map<int,vector<double>>,double> parameters;

    if(process == Dividends){
        get<0>(parameters) =div_mu;
        get<1>(parameters) =div_sig;
        get<2>(parameters) =div_aut;
        get<3>(parameters) =div_cross_corr;
        get<4>(parameters) =Div;
        get<5>(parameters) = trade_period;

    }else
    if(process == Earnings){
        get<0>(parameters) =earnings_mu;
        get<1>(parameters) =earnings_sig;
        get<2>(parameters) =earnings_aut;
        get<3>(parameters) =earnings_cross_corr;
        get<4>(parameters) =Earnings_;
        get<5>(parameters) = trade_period;
    }else
    if(process == Free_Cash_Flows){
        get<0>(parameters) =fcf_mu;
        get<1>(parameters) =fcf_sig;
        get<2>(parameters) =fcf_aut;
        get<3>(parameters) =fcf_cross_corr;
        get<4>(parameters) =CashFlow;
        get<5>(parameters) = trade_period;

    }
return parameters;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static vector<double> load_process(map<int,vector<double>> P_o, int trade_period, int tickerID, VectorXd mu, VectorXd sig, double dt, MatrixXd dU){
    vector<double> process1;
    vector<double> div_process;
    vector<double> div_growth;
std::ofstream noise;
noise.open("noisefile.csv");
    noise<<"growth"<<","<<"D_t_1"<<","<<"D_t"<<","<< "dU(tickerIDt)"<<","<<"sig"<<"mu"<<","<<","<<"altered_sig"<<std::endl;




    int size = ceil(trade_period);

    for(int t = 0; t < size; t++) {
        if (process1.empty()) {
            double a = P_o.find(1)->second.at(tickerID);

            process1.push_back(a);



                div_process.push_back(a);


        }else {





            double y  = (process1[t-1]) + ((mu(tickerID))*process1[t-1]) + (sig(tickerID) * process1[t-1] * dU(tickerID, t));


            noise<<""<<","<<process1[t-1]<<","<<y<<","<< dU(tickerID, t)<<","<<sig(tickerID)<<","<<(mu(tickerID) - (pow(sig(tickerID),2.))/2.) * (1/252.)<<","<<(sig(tickerID) * sqrt(1/252.) * dU(tickerID, t))<<std::endl;


            process1.push_back(y);


            div_process.push_back(y);

                }
    }
    noise.close();
    return div_process;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static vector<vector<double>> generate(processes process){

    int trade_period = get<5>(params(process));
//we set cross-sectional correlation

    MatrixXd divCrossCorr(3, 3);
    divCrossCorr = get<3>(params(process));

//set the autocorrelations for the individual assets
    VectorXd Aut = VectorXd(3);
    Aut = get<2>(params(process));


//volatility
    VectorXd sig = VectorXd(3);
    sig = get<1>(params(process));



//mean
    VectorXd mu = VectorXd(3);
    mu = get<0>(params(process));

//initial value
    map<int,vector<double>> P_o = get<4>(params((process)));

    MatrixXd dU(3,trade_period);
    dU = stats::generateColoredNoise(trade_period, Aut, divCrossCorr);


    vector<double> x1 = load_process(P_o, trade_period,0,mu,sig,1,dU);
    vector<double> x2 = load_process(P_o,trade_period,1,mu,sig,1,dU);
    vector<double> x3 = load_process(P_o, trade_period,2,mu,sig,1,dU);

    vector<vector<double>> result;
    result.push_back(x1);
    result.push_back(x2);
    result.push_back(x3);

    return result;

}


#endif //UNTITLED32_WIENER_H
