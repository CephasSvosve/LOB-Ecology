//
// Created by Cephas Svosve on 5/1/2022.
//

#ifndef UNTITLED32_WIENER_H
#define UNTITLED32_WIENER_H

#include "stochastic_math.h"
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enum processes{Earnings, Dividends, Free_Cash_Flows};

static tuple<VectorXd,VectorXd,VectorXd,MatrixXd, double, double> params(processes process){
//trade period
    double trade_period = 1000;



//cross-sectional correlations
    MatrixXd div_cross_corr(3, 3);
    div_cross_corr << 1, 0.251, 0.034,
            0.034, 1, 0.3,
            0.3, 0.034,1;

    MatrixXd earnings_cross_corr(3, 3);
    earnings_cross_corr <<  1, 0.251, 0.034,
            0.034, 1, 0.3,
            0.3, 0.034,1;

    MatrixXd fcf_cross_corr(3, 3);
    fcf_cross_corr <<       1, 0.251, 0.034,
            0.034, 1, 0.3,
            0.3, 0.034,1;


//auto-correlations
    VectorXd div_aut = VectorXd(3);
    div_aut(0) =        0.02;
    div_aut(1) =        0.13;
    div_aut(2) =        0.36;

    VectorXd earnings_aut = VectorXd(3);
    earnings_aut(0) =   0.02;
    earnings_aut(1) =   0.13;
    earnings_aut(2) =   0.36;

    VectorXd fcf_aut = VectorXd(3);
    fcf_aut(0) =        0.02;
    fcf_aut(1) =        0.13;
    fcf_aut(2) =        0.36;

//mu
    VectorXd div_mu = VectorXd(3);
    div_mu(0) = 000.42;
    div_mu(1) = 000.53;
    div_mu(2) = 000.86;

    VectorXd earnings_mu = VectorXd(3);
    earnings_mu(0) = (1./63) * 3.23;
    earnings_mu(1) = (1./63) * 2.54;
    earnings_mu(2) = (1./63) * 1.86;

    VectorXd fcf_mu = VectorXd(3);
    fcf_mu(0) = 0.01;
    fcf_mu(1) = 0.02;
    fcf_mu(2) = 0.03;

//sigma
    VectorXd div_sig = VectorXd(3);
    div_sig(0) =  001.046319;
    div_sig(1) =  02.7746319;
    div_sig(2) =  03.7746319;



    VectorXd earnings_sig = VectorXd(3);
    earnings_sig(0) = (1./63) * 12.62;
    earnings_sig(1) = (1./63) * 9.23;
    earnings_sig(2) = (1./63) * 7.42;

    VectorXd fcf_sig = VectorXd(3);
    fcf_sig(0) = 001.2;
    fcf_sig(1) = 002.3;
    fcf_sig(2) = 007.6;

    //Initial values
    double E_o = 0.4048; //TODO add initial values for dividends and fcf
    double D_o = 0;//

    tuple<VectorXd,VectorXd,VectorXd,MatrixXd, double,double> parameters;

    if(process == Dividends){
        get<0>(parameters) =div_mu;
        get<1>(parameters) =div_sig;
        get<2>(parameters) =div_aut;
        get<3>(parameters) =div_cross_corr;
        get<4>(parameters) =D_o;
        get<5>(parameters) = trade_period;

    }else
    if(process == Earnings){
        get<0>(parameters) =earnings_mu;
        get<1>(parameters) =earnings_sig;
        get<2>(parameters) =earnings_aut;
        get<3>(parameters) =earnings_cross_corr;
        get<4>(parameters) =E_o;
        get<5>(parameters) = trade_period;
    }else
    if(process == Free_Cash_Flows){
        get<0>(parameters) =fcf_mu;
        get<1>(parameters) =fcf_sig;
        get<2>(parameters) =fcf_aut;
        get<3>(parameters) =fcf_cross_corr;
        get<4>(parameters) =E_o;
        get<5>(parameters) = trade_period;

    }
    return parameters;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static vector<double> load_process(double P_o, int trade_period, int tickerID, VectorXd mu, VectorXd sig, double dt, MatrixXd dU){
    vector<double> process1;
  // distribution in range [1, 6]

    for(int i = 0; i < trade_period; i++) {
        if (process1.empty()) {
            double a = P_o ;

            process1.push_back(a);
        }else {
            double y = (process1[i - 1] )+ mu(tickerID)  * dt + sig(tickerID) * dU(tickerID, i);
            process1.push_back(y);
        }
    }
    return process1;
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
    double P_o = get<4>(params((process)));

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
