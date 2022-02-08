//
// Created by Cephas Svosve on 21/12/2021.
//

#ifndef UNTITLED30_FUNDS_H
#define UNTITLED30_FUNDS_H
#include "value_investor.h"
#include "momentum_investor.h"
#include "noisy_investor.h"
#include <tuple>
#include <boost/random.hpp>
#include <map>
#include <Eigen/Dense>
#include "market_watch.h"
#include "stochastic_math.h"
#include <random>





using namespace std;

class funds
        :public fund
     {
private:
    int window_size = 200;

public:
funds();

enum trading_strategy{fundamental_value=0
//                        ,growth = 1
                            , noise=1
                                , momentum_investment=2

                                       };
trading_strategy fund_philosophy;

void trade_strategy(trading_strategy tradingStrategy);


map<int, order>
invest();


double
compute_earnings(company &stock, int &time, int &reb_period);


MatrixXd
generateWhiteNoise1(int seed, int rows, int columns);


double
lateralcorrcoef(VectorXd a);


std::map<int,order>
value_demand();


std::map<int,order>
growth_demand();

std::map<int,order>
momentum_demand();


std::map<int,order>
noise_demand();



};
#endif //UNTITLED30_FUNDS_H
