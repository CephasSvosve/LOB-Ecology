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
#include "order.h"





using namespace std;

class funds
        :public fund
     {
private:
    int window_size = 200;

public:
funds();

enum trading_strategy{aggressive_fundamental_value=0
                            , aggressive_noise=1
                                , aggressive_momentum_investment=2
                                    , aggressive_growth = 3
                                        , aggressive_index =4
                                            , fundamental_value =5
                                                , noise = 6
                                                    , momentum_investment = 7
                                                        , growth = 8
                                                            ,index = 9


                                       };




trading_strategy fund_philosophy;

void trade_strategy(trading_strategy tradingStrategy);


map<int, order>
invest();


static int concat(int a, int b);

static MatrixXd
generateWhiteNoise1(int seed, int rows, int columns);


static double
lateralcorrcoef(VectorXd a);


std::map<int,order>
ag_value_demand(),



ag_growth_demand(),


ag_momentum_demand(),



ag_noise_demand(),



ag_index_demand(),



value_demand(),



growth_demand(),


momentum_demand(),



noise_demand(),



index_demand();




};
#endif //UNTITLED30_FUNDS_H
