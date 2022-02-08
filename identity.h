//
// Created by Cephas Svosve on 20/12/2021.
//

#ifndef UNTITLED30_IDENTITY_H
#define UNTITLED30_IDENTITY_H
#include <map>

using namespace std;


enum identity_type{company_ = 0, asset= 1, trader= 2, fund_trader= 3, pricer= 4, strategy= 5};

static map<int,identity_type> companies;
static map<int,identity_type> assets;
static map<int,identity_type> traders;
static map<int,identity_type> fund_traders;
static map<int,identity_type> pricers;
static map<int,identity_type> trade_strategies;

int create_identity(identity_type a);


#endif //UNTITLED30_IDENTITY_H
