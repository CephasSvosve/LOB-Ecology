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
    map<int,vector<order>> bid;
    map<int,double> inventory_t_1;
    map<int,double> change_in_inventory;


public:
    map<int, funds*> trading_institutions;
    map<int,company>
    &tradeable_assets();

    double order_size = 100.;

    void register_assets(int number_of_assets);
    void generate_fundamentals();
    void register_traders(vector<funds*> participants);
    void update_balances(int t);
    void set_initial_quotes();
    void receive_orders();
    void stack_market_orders(int id, order order_);
    void stack_bid(int id, order order_);
    void stack_ask(int id, order order_);
    void clear(market_watch &watch);
    void update_price(int k, double reg);

    tuple<order,order>
    update_mm_orders(const map<double, vector<order>, greater<>> &bid_orders_
            ,const map<double, vector<order>> &ask_orders_
            ,double init_inventory
            ,double curr_inventory
            ,double curr_price
            ,double init_price
            ,int identifier
            ,int asset_id) const;

    bool not_similar(vector<order> &x);
    bool isValid(int asset_id);

    std::ofstream ordersfile;


};


#endif //UNTITLED32_PRICE_SETTER_H




////
//// Created by Cephas Svosve on 5/1/2022.
////
//
//#ifndef UNTITLED32_PRICE_SETTER_H
//#define UNTITLED32_PRICE_SETTER_H
//#include <map>
//#include <math.h>
//#include "funds.h"
//#include "price.h"
//#include "market_watch.h"
//#include "order.h"
//#include "company.h"
//#include "wiener.h"
//#include <algorithm>
//
//using namespace std;
//
//
//class price_setter
//        :public funds {
//
//private:
//    map<int,company> assets;
//
//
//
//
//    map<int,vector<order>> ask;
//    map<int,vector<order>> market_orders;
//    map<int,double> best_bid;
//    map<int,double> best_ask;
//    map<int,double> mid_price;
//    map<int,vector<order>> bid;
//
//
//public:
//    map<int, funds*> trading_institutions;
//    map<int,company>
//            &tradeable_assets();
//
//    void register_assets(int number_of_assets);
//    void generate_fundamentals(double interest_rate);
//    void register_traders(vector<funds*> participants);
//    void update_balances(int t);
//    void set_initial_quotes();
//    void receive_orders();
//    void stack_market_orders(int id, order order_);
//    void stack_bid(int id, order order_);
//    void stack_ask(int id, order order_);
//    void clear();
//    void update_price(int k, double reg);
//
//
//    bool not_similar(vector<order> &x);
//    bool isValid(int asset_id);
//
//};
//
//
//#endif //UNTITLED32_PRICE_SETTER_H
