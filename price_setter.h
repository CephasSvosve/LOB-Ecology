//
// Created by Cephas Svosve on 5/1/2022.
//

#ifndef UNTITLED32_PRICE_SETTER_H
#define UNTITLED32_PRICE_SETTER_H
#include <map>
#include <math.h>
#include "funds.h"
#include "price.h"
#include "market_watch.h"
#include "order.h"
#include "company.h"
#include "wiener.h"
#include <algorithm>

using namespace std;


class price_setter
        :public funds {

private:
    map<int,company> assets;




    map<int,vector<order>> ask;
    map<int,vector<order>> market_orders;
    map<int,double> best_bid;
    map<int,double> best_ask;
    map<int,double> mid_price;


public:
    map<int, funds*> trading_institutions;
    map<int,company>
            &tradeable_assets();
    map<int,vector<order>> bid;
    void register_assets(int number_of_assets);
    void generate_fundamentals();
    void register_traders(vector<funds*> participants);
    void update_balances(int t);
    void set_initial_quotes();
    void receive_orders();
    void stack_market_orders(int id, order order_);
    void stack_bid(int id, order order_);
    void stack_ask(int id, order order_);
    void clear();
    void update_price(int k);


    bool not_similar(vector<order> &x);
    bool isValid(int asset_id);

};


#endif //UNTITLED32_PRICE_SETTER_H
