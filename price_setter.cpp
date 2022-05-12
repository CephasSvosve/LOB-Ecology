




//
// Created by Cephas Svosve on 5/1/2022.
//

#include "price_setter.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <random>


template <typename T>
int sgn(T val) {
    return (0 < val) - (val < 0);
}


//method for accessing assets registered on the market
map<int,company>&
price_setter::tradeable_assets(){
    return this->assets;
}


//this lets us create stocks that will be traded on the market
void price_setter::register_assets(int number_of_assets){
    for(int i = 0; i < number_of_assets; i++){
        if(this->assets.size() < number_of_assets){
        company a;
        a.set_shares_outstanding(1.68e10*0.4/75);
        this->assets.emplace(a.get_identifier(),a);
        }
    }
}


//this generates fundamental processes and associates them with each asset
void price_setter::generate_fundamentals(){
    double day_cost_of_equity = pow((1.1),(1./252))-1;

//        auto earnings_processes  = generate(Earnings);
            auto dividends_processes = generate(Dividends);
//                auto free_cash_flow      = generate(Free_Cash_Flows);
                    int count = 0;

                //associate fundamental processes with assets

        for(auto &[i,j] : this->assets){
//            j.set_earnings_process(earnings_processes.at(count));
//                j.set_free_cash_flow_process(free_cash_flow.at(count));
                    j.set_dividends_process(dividends_processes.at(count));

                        //preset intrinsic value of the asset
                        //(1)....compute present value of each dividend
                        int size = j.get_dividend_process().size();


                        vector<double> iv_aggregate={};


                        iv_aggregate.reserve(size);

                        //growth rate g
                        double g =0;
                        for(int d=1; d < size; d++ ){

                         g += log(j.get_dividend_process().at(d))- log(j.get_dividend_process().at(d-1));

                        }


                         g = g/size;


                            for(int d = 0; d < size; d++){
                                double iv = 0;
                                iv = j.get_dividend(d) /(day_cost_of_equity-g);
                                iv_aggregate.push_back(iv);

                            }







                        j.set_intrinsic_value(iv_aggregate);


        count ++;
    }



}


//this notifies the price setter of the available traders
void price_setter::register_traders(vector<funds*> i) {
    for (auto &p : i) {
        this->trading_institutions.emplace(p->get_identifier(), p);
    }
}


//this updates traders' balances before sending new quotes
void price_setter::update_balances(int t)
{
    for(auto &[k,v]: this->trading_institutions){
        v->balance_cf(t);
    }
}


void price_setter::set_initial_quotes(){
    for(auto &[k,v] : this->assets){
        double init_price;
        std::cout<<" initial price for asset "<< k <<" is: USD";
        std::cin >> init_price;
        v.set_price(0,init_price,init_price);


//        this->bid.emplace(k,vector<order>());
//        this->ask.emplace(k,vector<order>());
//        this->market_orders.emplace(k,vector<order>());

        this->best_bid.emplace(k,init_price);
        this->best_ask.emplace(k,init_price);
        this->mid_price.emplace(k,init_price);
    }
    this->get_inform(this->assets);
}


void price_setter::receive_orders(){
    this->bid.clear();
        this->ask.clear();
            this->market_orders.clear();
                for(auto &[k,v] : this->trading_institutions){
                    if ((v->wealth) > 0) {
                        auto order_ = v->invest();
                            for (auto &[i, j] : order_) {
                                if (abs(j.get_order_size()) > 0) {
                                    if (j.get_order_type()) {//this checks whether the order is a limit(1) or market(0) order
                                        if (j.get_order_size() > 0) {//this checks whether order is a bid(+) or ask(-)
                                            stack_bid(j.get_ordered_asset(), j);
                                        } else {
                                            stack_ask(j.get_ordered_asset(), j);
                                        }
                                    } else {//if order is a market order we stack the market order vector
                                        stack_market_orders(j.get_ordered_asset(), j);
                                    }
                                }
                            }
                    }else{
            for(auto &[i,j]: v->stocks_at_hand){
                if(j<0) {
                    this->stocks_at_hand.find(i)->second += -j;
                    v->stocks_at_hand.find(i)->second = 0;
                }else{
                    this->stocks_at_hand.find(i)->second += j;
                    v->stocks_at_hand.find(i)->second = 0;
                }
            }
            v->bond_at_hand = 0;
            v->cash_at_hand = 0;
        }
    }


}


void price_setter::stack_market_orders(int id, order order_){
    auto k = id;
        auto i = this->market_orders.find(k);
            if(i != this->market_orders.end()){
                this->market_orders.find(k)->second.push_back(order_);
            }else{
            vector<order> a;
        a.push_back(order_);
    this->market_orders.emplace(k,a);
    }
}


void price_setter::stack_bid(int id, order order_){
    auto k = id;
        auto  u= this->bid.find(k);
            if(u != this->bid.end()){
                    this->bid.find(k)->second.push_back(order_);
                }else{
            vector<order> a;
         a.push_back(order_);
    this->bid.emplace(k,a);
    }
}


void price_setter::stack_ask(int id, order order_){
    auto k = id;
        auto c = this->ask.find(k);
            if(c != this->ask.end()){
                        this->ask.find(k)->second.push_back(order_);
                       }else{
            vector<order> a;
        a.push_back(order_);
    this->ask.emplace(k,a);
    }
}


//this method sorts compares orders by order prices in ascending order
bool ascending_order(order x, order y){
    return ((x.get_proposed_price()) < (y.get_proposed_price()));
}


bool descending_order(order x, order y){
    return ((x.get_proposed_price()) > (y.get_proposed_price()));
}


//this method sorts orders by order sizes in descending order
bool sort_by_size(order x, order y){
    return ((x.get_order_size()) > (y.get_order_size()));
}


//this method sorts orders by order id in decreasing order in descending order
bool sort_by_id(order x, order y){
    return ((x.get_id()) > (y.get_id()));
}


tuple<order,order>
match_market_orders(order &market_order, order &bid_order,order &ask_order){
tuple<order,order> result;


    if (market_order.get_order_size() < 0.){



                double bid_price = (floor(bid_order.get_proposed_price()*100))/100;


                    //do the standard market order to bid order
                    bool bid_is_active = bid_order.get_status();
                    bool mktask_is_active = market_order.get_status();

                    if (bid_is_active && mktask_is_active) {
                        //trade

                        double mkt_cash = market_order.get_order_size() * market_order.get_proposed_price();
                        double BID = bid_order.get_order_size();
                        double ASK = mkt_cash / bid_price;
                        std::cout<<market_order.get_id()<< " "<<bid_order.get_id()<<" price "<<bid_price <<" bid "<<BID<< " ask "<<ASK<<std::endl;
                        auto uncleared = BID + ASK;
                        std::cout<<"passed 1 "<<std::endl;
//solve for the 3 situations: uncleared <, > or = 0;
                        if (uncleared < 0.) {
                            //calculate executed bid
                            order exec_bid;
                            auto exec_bid_price = bid_price;
                            auto bid_trader_id = bid_order.get_id();

                            exec_bid.set_id(bid_trader_id);
                            exec_bid.set_status(order::active);
                            exec_bid.set_order_type(order::limit);
                            exec_bid.set_ordered_asset(bid_order.get_ordered_asset());
                            exec_bid.set_order_size(BID, exec_bid_price);

                            get<0>(result) = exec_bid;


                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = bid_price;
                            auto ask_trader_id = market_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(market_order.get_ordered_asset());
                            exec_ask.set_order_size(-BID, exec_ask_price);

                            get<1>(result) = exec_ask;



                            //update the ask orders remaining in the order vectors after execution of trade
                            market_order.set_order_size((uncleared * bid_price) / market_order.get_proposed_price(),
                                                        market_order.get_proposed_price());
                            //remove the fully filled bid order
                            bid_order.set_status(order::filled);

                        } else if (uncleared > 0.) {
                            //calculate executed bid
                            order exec_bid;
                            auto exec_bid_price = bid_price;
                            auto bid_trader_id = bid_order.get_id();

                            exec_bid.set_id(bid_trader_id);
                            exec_bid.set_status(order::active);
                            exec_bid.set_order_type(order::limit);
                            exec_bid.set_ordered_asset(bid_order.get_ordered_asset());
                            exec_bid.set_order_size((BID - uncleared), exec_bid_price);

                            get<0>(result) = exec_bid;



                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = bid_price;
                            auto ask_trader_id = market_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(market_order.get_ordered_asset());
                            exec_ask.set_order_size(-(BID - uncleared), exec_ask_price);

                            get<1>(result) = exec_ask;


                            bid_order.set_order_size(uncleared, bid_price);
                            market_order.set_status(order::filled);

                        } else {
                            //calculate executed bid
                            order exec_bid;
                            auto exec_bid_price = bid_price;
                            auto bid_trader_id = bid_order.get_id();

                            exec_bid.set_id(bid_trader_id);
                            exec_bid.set_status(order::active);
                            exec_bid.set_order_type(order::limit);
                            exec_bid.set_ordered_asset(bid_order.get_ordered_asset());
                            exec_bid.set_order_size(BID, exec_bid_price);

                            get<0>(result) = exec_bid;


                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = bid_price;
                            auto ask_trader_id = market_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(market_order.get_ordered_asset());
                            exec_ask.set_order_size(-BID, exec_ask_price);

                            get<1>(result) = exec_ask;


                            bid_order.set_status(order::filled);
                            market_order.set_status(order::filled);
                        }
                        bid_is_active = bid_order.get_status();
                        mktask_is_active = market_order.get_status();
                    }

        std::cout<<"passed 2 "<<std::endl;

            }else
                if (market_order.get_order_size() > 0.) {



                double ask_price = (floor(ask_order.get_proposed_price()*100))/100;


                    bool ask_is_active = ask_order.get_status();
                    bool mktbid_is_active = market_order.get_status();


                    if (ask_is_active && mktbid_is_active) {
                        //trade
                        double mkt_cash = market_order.get_order_size() * market_order.get_proposed_price();
                        double BID = mkt_cash / (ask_price);
                        double ASK = ask_order.get_order_size();
                        auto uncleared = BID + ASK;



//solve for the 3 situations: uncleared <, > or = 0;
                        if (uncleared < 0.) {
                            //calculate executed bid
                            order exec_bid;
                            auto exec_bid_price = ask_price;
                            auto bid_trader_id = market_order.get_id();

                            exec_bid.set_id(bid_trader_id);
                            exec_bid.set_status(order::active);
                            exec_bid.set_order_type(order::limit);
                            exec_bid.set_ordered_asset(market_order.get_ordered_asset());
                            exec_bid.set_order_size(BID, exec_bid_price);

                            get<0>(result) = exec_bid;

                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = ask_price;
                            auto ask_trader_id = ask_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(ask_order.get_ordered_asset());
                            exec_ask.set_order_size(-BID, exec_ask_price);

                            get<1>(result) = exec_ask;



                            //update the ask orders remaining in the order vectors after execution of trade
                            ask_order.set_order_size(uncleared, ask_price);
                            //remove the fully filled bid order
                            market_order.set_status(order::filled);

                        } else if (uncleared > 0.) {
                            //calculate executed bid
                            order exec_bid;
                            auto exec_bid_price = ask_price;
                            auto bid_trader_id = market_order.get_id();

                            exec_bid.set_id(bid_trader_id);
                            exec_bid.set_status(order::active);
                            exec_bid.set_order_type(order::limit);
                            exec_bid.set_ordered_asset(market_order.get_ordered_asset());
                            exec_bid.set_order_size((BID - uncleared), exec_bid_price);

                            get<0>(result) = exec_bid;



                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = ask_price;
                            auto ask_trader_id = ask_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(ask_order.get_ordered_asset());
                            exec_ask.set_order_size(-(BID - uncleared), exec_ask_price);

                            get<1>(result) = exec_ask;

                            market_order.set_order_size(
                                    (uncleared * ask_price) / market_order.get_proposed_price(),
                                    market_order.get_proposed_price());
                            ask_order.set_status(order::filled);


                        } else {
                            //calculate executed bid
                            order exec_bid;
                            auto exec_bid_price = ask_price;
                            auto bid_trader_id = market_order.get_id();

                            exec_bid.set_id(bid_trader_id);
                            exec_bid.set_status(order::active);
                            exec_bid.set_order_type(order::limit);
                            exec_bid.set_ordered_asset(market_order.get_ordered_asset());
                            exec_bid.set_order_size(BID, exec_bid_price);

                            get<0>(result) = exec_bid;




                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = ask_price;
                            auto ask_trader_id = ask_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(ask_order.get_ordered_asset());
                            exec_ask.set_order_size(-BID, exec_ask_price);

                            get<1>(result) = exec_ask;


                            market_order.set_status(order::filled);
                            ask_order.set_status(order::filled);
                        }
                    }
                    std::cout<<"passed 3 "<<std::endl;
        }

return result;
}


tuple<order,order>
match_limit_orders(order &bid_order,order &ask_order){
    tuple<order,order>  result;

    double bid_price = (floor(ask_order.get_proposed_price()*100))/100;//assuming the ask order arrived first
    double ask_price = (floor(ask_order.get_proposed_price()*100))/100;

    bool bid_is_active = bid_order.get_status();
    bool ask_is_active = ask_order.get_status();

    //non-stale orders
    if (bid_is_active && ask_is_active) {
        //trade
        std::cout<<"passed 4 "<<std::endl;
        double BID = bid_order.get_order_size();
        double ASK = ask_order.get_order_size();
        auto uncleared = BID + ASK;

        //solve for the 3 situations: uncleared <, > or = 0;
        if (uncleared < 0.) {
            //calculate executed bid
            order exec_bid;
            auto exec_bid_price = bid_price;
            auto bid_trader_id = bid_order.get_id();

            exec_bid.set_id(bid_trader_id);
            exec_bid.set_status(order::active);
            exec_bid.set_order_type(order::limit);
            exec_bid.set_ordered_asset(bid_order.get_ordered_asset());
            exec_bid.set_order_size(BID, exec_bid_price);

            get<0>(result) = exec_bid;



            //calculate executed ask
            order exec_ask;
            auto exec_ask_price = ask_price;
            auto ask_trader_id = ask_order.get_id();

            exec_ask.set_id(ask_trader_id);
            exec_ask.set_status(order::active);
            exec_ask.set_order_type(order::limit);
            exec_ask.set_ordered_asset(ask_order.get_ordered_asset());
            exec_ask.set_order_size(-BID, exec_ask_price);

            get<1>(result) = exec_ask;


            //update the ask orders remaining in the order vectors after execution of trade
            ask_order.set_order_size(uncleared, ask_price);
            //remove the fully filled bid order
            (&bid_order)->set_status(order::filled);
            std::cout<<"passed 5 "<<std::endl;
        } else if (uncleared > 0.) {
            //calculate executed bid
            order exec_bid;
            auto exec_bid_price = bid_price;
            auto bid_trader_id = bid_order.get_id();

            exec_bid.set_id(bid_trader_id);
            exec_bid.set_status(order::active);
            exec_bid.set_order_type(order::limit);
            exec_bid.set_ordered_asset(bid_order.get_ordered_asset());
            exec_bid.set_order_size((BID - uncleared), exec_bid_price);

            get<0>(result) = exec_bid;




            //calculate executed ask
            order exec_ask;
            auto exec_ask_price = ask_price;
            auto ask_trader_id = ask_order.get_id();

            exec_ask.set_id(ask_trader_id);
            exec_ask.set_status(order::active);
            exec_ask.set_order_type(order::limit);
            exec_ask.set_ordered_asset(ask_order.get_ordered_asset());
            exec_ask.set_order_size(-(BID - uncleared), exec_ask_price);

            get<1>(result) = exec_ask;


            bid_order.set_order_size(uncleared, bid_price);
            (&ask_order)->set_status(order::filled);
            std::cout<<"passed 6 "<<std::endl;
        } else {
            //calculate executed bid
            order exec_bid;
            auto exec_bid_price = bid_price;
            auto bid_trader_id = bid_order.get_id();

            exec_bid.set_id(bid_trader_id);
            exec_bid.set_status(order::active);
            exec_bid.set_order_type(order::limit);
            exec_bid.set_ordered_asset(bid_order.get_ordered_asset());
            exec_bid.set_order_size(BID, exec_bid_price);

            get<0>(result) = exec_bid;



            //calculate executed ask
            order exec_ask;
            auto exec_ask_price = ask_price;
            auto ask_trader_id = ask_order.get_id();

            exec_ask.set_id(ask_trader_id);
            exec_ask.set_status(order::active);
            exec_ask.set_order_type(order::limit);
            exec_ask.set_ordered_asset(ask_order.get_ordered_asset());
            exec_ask.set_order_size(-BID, exec_ask_price);

            get<1>(result) = exec_ask;

            //update cash earned from spread
            (&bid_order)->set_status(order::filled);
            (&ask_order)->set_status(order::filled);

        }
//        ask_is_active = ask_order.get_status();
//        bid_is_active =  bid_order.get_status();
        std::cout<<"passed 7 "<<std::endl;
    }

    return result;
}


map<int, vector<order>>
record_executed_orders(tuple<order,order> &settlements){

    map<int, vector<order>> executed_orders_;

    order order1 = get<0>(settlements);
    order order2 = get<1>(settlements);

    if (executed_orders_.find(order1.get_id()) != executed_orders_.end()){
        executed_orders_.find(order1.get_id())->second.push_back(order1);
    }else{
        executed_orders_.emplace(order1.get_id(), vector<order>{order1});
    }


    if (executed_orders_.find(order2.get_id()) != executed_orders_.end()){
        executed_orders_.find(order2.get_id())->second.push_back(order2);
    }else{
        executed_orders_.emplace(order2.get_id(), vector<order>{order2});
    }

    return executed_orders_;
}


order *get_best_bid(map<double, vector<order>, greater<>> &bids){
    order *result = nullptr;
    for(auto &[b,vec]: bids){
        for(auto &x:vec){
            if (x.get_status() == 1){
                result = &x;
                goto exit;
            }
        }
    }
    exit:;
    return result;
}


order *get_best_ask(map<double, vector<order>> &asks){
    order *result = nullptr;
    for(auto &[a,vec]: asks){
        for(auto &x:vec){
            if (x.get_status() == 1){
                result = &x;
                goto exit;
            }
        }
    }
    exit:;
    return result;
}


void cancel_mm_orders(int order_id
                    , map<double, vector<order>, greater<>> &bids
                                ,map<double, vector<order>> &asks){


for(auto &[price, vec]: bids){
    for(auto &x : vec){
        if(x.get_id() == order_id){
            x.set_status(order::filled);
        }
    }
}

for(auto &[price, vec]: asks){
    for(auto &x : vec){
        if(x.get_id() == order_id){
            x.set_status(order::filled);
        }
    }
}

}


void add_mm_orders(tuple<order,order> new_orders
                  ,map<double, vector<order>, greater<>> &bids
                     ,map<double, vector<order>> &asks){

    double bid_price = get<0>(new_orders).get_proposed_price();
    bid_price = floor(bid_price*100)/100;

    double ask_price = get<1>(new_orders).get_proposed_price();
    ask_price = floor(ask_price*100)/100;


    //add bid
    if (bids.find(bid_price) != bids.end()){
        bids.find(bid_price)->second.push_back(get<0>(new_orders));
        }
    else{
            bids.emplace(bid_price, vector{get<0>(new_orders)});
            }

    //add ask
    if (asks.find(ask_price) != asks.end()){
        asks.find(ask_price)->second.push_back(get<1>(new_orders));
    }else{
        asks.emplace(ask_price, vector{get<1>(new_orders)});
    }


}


tuple<order,order>
price_setter::update_mm_orders(const map<double, vector<order>, greater<>> &bid_orders_
                    ,const map<double, vector<order>> &ask_orders_
                        ,double init_inventory
                            ,double curr_inventory
                                ,double curr_price
                                    ,double init_price
                                        ,int identifier
                                            ,int asset_id) const{

    tuple<order,order>  new_order;

    double best_ask_;
    double best_bid_;

    //search for best bid
    for(auto &[price, vec]: bid_orders_){
        int x = 0;
        for(auto &o:vec){
            if (!o.filled){
                x++;
            }
        }
        if(x>0){
            best_bid_ = price;
            break;
        }
    }

    //search for best ask
    for(auto &[price, vec]: ask_orders_){
        int y = 0;
        for(auto &o:vec){
            if (!o.filled){
                y++;
            }
        }
        if(y>0){
            best_ask_ = price;
            break;
        }
    }



    double target_bid;
    double target_ask;

    double new_bid;
    double new_ask;

    order sell_limit;
    order buy_limit;

//compute target bid and target ask
    double d_shares   =  curr_inventory - init_inventory;

    target_ask =  init_price-((d_shares/order_size)*0.01);


    target_bid =  target_ask-0.05 ;//market maker's bid price is always atleast 5 basis points below the ask price


//setting a new_ask price give a target ask


        new_ask = target_ask;



//setting a new bid_price given a target bid

        new_bid = target_bid;

//round order prices  downwards to two dp
        new_ask = (floor(new_ask*100))/100;
        new_bid = (floor(new_bid*100))/100;

//create order
    //double cash_ = cash_at_hand/3;

    //buy order
    buy_limit.set_id(identifier);
    buy_limit.set_order_type(order::limit);
    buy_limit.set_ordered_asset(asset_id);
    buy_limit.set_order_size(order_size, new_bid);
    buy_limit.set_status(order::active);
    get<0>(new_order) = buy_limit;

    //sell order
    sell_limit.set_id(identifier);
    sell_limit.set_order_type(order::limit);
    sell_limit.set_ordered_asset(asset_id);
    sell_limit.set_order_size(-order_size, new_ask);
    sell_limit.set_status(order::active);
    get<1>(new_order) = sell_limit;


    return new_order;
}


void
price_setter::clear(market_watch &market_clock){


//Params definition
    map<int, order> ask_result_ = {};
    map<int, order> bid_result_ = {};


//map<asset, map<price,order>> market_maker_orders
    map<int, map<int, order>> mm_orders = {};


    auto quotes = this->assets;


    auto t = clock->current_time();
    std::cout << std::endl;
    std::cout << "time " << t << std::endl;
    if (int(t) == t) {
        this->balance_cf(t);
    }

    double short_div = 0.;
    for (auto &[j, a]:this->trading_institutions) {
        for (auto &[i, b]: stocks_on_market) {
            if ((a->stocks_at_hand.find(i) != a->stocks_at_hand.end())
                && a->stocks_at_hand.find(i)->second < 0.) {
                short_div += a->stocks_at_hand.find(i)->second * -b.get_dividend(t);
            }
        }
    }

    this->cash_at_hand += short_div;


    {
//compute the ratio of each stock's signal relative to sum described in (1)
//    for(auto &[s, v] : quotes){
//
//
//
//        order buy_limit;
//        order sell_limit;
//        double non_overlapping_ask =0;
//        double non_overlapping_bid =0;
//        double bid_ =0;
//        double ask_ =0;
//
////get the highest bid price,however, if there are no new better bid quotes we maintain bid at the recent last bid
//        if(this->bid.find(s)!= this->bid.end()){
//
//            sort(this->bid.find(s)->second.begin(),this->bid.find(s)->second.end(),descending_order);
//            bid_ = this->bid.find(s)->second.begin()->get_proposed_price();
//            non_overlapping_bid = bid_;
//
//        }else{
//
//            if (this->ask.find(s) != this->ask.end()) {
//                sort(this->ask.find(s)->second.begin(), this->ask.find(s)->second.end(), ascending_order);
//                bid_ = this->ask.find(s)->second.begin()->get_proposed_price()-0.01;
//                non_overlapping_bid = bid_;
//
//            } else{
//                bid_= get<1>(v.get_price());
//                non_overlapping_bid =bid_;
//            }
//
//        }
//
////get the lowest ask price, if there are no new asks, maintain ask at the recent last ask
//        if (this->ask.find(s) != this->ask.end()) {
//            sort(this->ask.find(s)->second.begin(), this->ask.find(s)->second.end(), ascending_order);
//            ask_ = this->ask.find(s)->second.begin()->get_proposed_price();
//            non_overlapping_ask= ask_;
//
//        }else {
//
//            if(this->bid.find(s)!= this->bid.end()){
//
//                sort(this->bid.find(s)->second.begin(),this->bid.find(s)->second.end(),descending_order);
//                ask_ = this->bid.find(s)->second.begin()->get_proposed_price()+0.01;
//                non_overlapping_ask = ask_;
//
//            } else{
//                ask_= get<2>(v.get_price());
//                non_overlapping_ask = ask_;
//            }
//        }
//
//
////iterate through asks to find the first ask that is greater than the highest bid price
//        if (this->ask.find(s) != this->ask.end()) {
//            for (auto &i:this->ask.find(s)->second) {
//                if (i.get_proposed_price() >= bid_) {
//                    non_overlapping_ask = i.get_proposed_price()+0.01;
//                    break;
//                }
//                //if all asks are greater than highest bid, (i.e. the 'for-loop' does not break)
//                //then the ask is set to one basis point above the highest bid
//                if(this->bid.find(s)!= this->bid.end()) {
//                    non_overlapping_ask = this->bid.find(s)->second.begin()->get_proposed_price() + 0.01;
//                }else {
//                    non_overlapping_ask = bid_+0.01;
//                }
//            }
//        }
//
////iterate through bids to find the first bid that is smaller than the lowest ask price
//        if (this->bid.find(s) != this->bid.end()) {
//            for (auto &i:this->bid.find(s)->second) {
//                if (i.get_proposed_price() < ask_) {
//                    non_overlapping_bid = i.get_proposed_price()-0.01;
//                    break;
//                }
//                //if all bids are smaller than lowest ask, (i.e. the 'for-loop' does not break)
//                //then the bid is set to one basis point below the lowest ask
//                if(this->ask.find(s)!= this->ask.end()) {
//                    non_overlapping_bid = this->ask.find(s)->second.begin()->get_proposed_price()-0.01;
//                }else {
//                    non_overlapping_bid = ask_-0.01;
//                }
//            }
//        }
//
//
//
//        int short_stocks = 0;
//            for(auto &[j,k]:this->trading_institutions){
//                if ((k->stocks_at_hand.find(s) != k->stocks_at_hand.end())
//                        && k->stocks_at_hand.find(s)->second<0){
//
//                    short_stocks += k->stocks_at_hand.find(s)->second;
//                }
//            }
//
//
//        double init_price = 5;
//        const double q    = this->initial_invetory.find(s)->second;
//        double inventory_ = (this->stocks_at_hand.find(s)->second + short_stocks);
//        double d_shares   =  inventory_ - q;
//        double bid_price_ =  non_overlapping_bid;
//        double ask_price_ = max(non_overlapping_ask, 3 - d_shares*10e-09);//((pow(q,2))/(pow(q,2)-pow(d_shares,2)));
//
//        sell_limit.set_id(this->get_identifier());
//        sell_limit.set_order_type(order::limit);
//        sell_limit.set_ordered_asset(s);
//        sell_limit.set_order_size(-inventory_, ask_price_);
//        sell_limit.set_status(order::active);
//        //stack_ask(sell_limit.get_ordered_asset(), sell_limit);
//
//
//        double tot_price = 0;
//            for(auto &[x,v] : this->assets){
//                tot_price += v.get_midprice();
//            }
//        double this_price = v.get_midprice();
//        double k_proportion = this_price/tot_price;
//        double cash_ = this->cash_at_hand/3;
//
//
//        buy_limit.set_id(this->get_identifier());
//        buy_limit.set_order_type(order::limit);
//        buy_limit.set_ordered_asset(s);
//        buy_limit.set_order_size(cash_/bid_price_, bid_price_);
//        buy_limit.set_status(order::active);
//        //stack_bid(buy_limit.get_ordered_asset(), buy_limit);
//
//
    }


    for (auto &[i, j] : this->stocks_at_hand) {

        int short_stocks = 0;
        for (auto &[j, k]:this->trading_institutions) {
            if ((k->stocks_at_hand.find(i) != k->stocks_at_hand.end())
                && k->stocks_at_hand.find(i)->second < 0) {

                short_stocks += k->stocks_at_hand.find(i)->second;
            }
        }

        if (inventory_t_1.find(i) != inventory_t_1.end()) {
            inventory_t_1.find(i)->second = j + short_stocks;
        } else {
            inventory_t_1.emplace(i, j);
        }
    }



//receive orders

//furthest range a market maker is allowed to charge away from the NBBO price

    map<int, vector<order>> executed_orders_;
    for (auto &[k, stock]:assets) {
        double prevailing_bid = get<1>(stock.get_price());//*(1-regulation_range_);
        double prevailing_ask = get<2>(stock.get_price());//*(1+regulation_range_);



        int short_stocks = 0;
        for (auto &[j, x]:this->trading_institutions) {
            if ((x->stocks_at_hand.find(k) != x->stocks_at_hand.end())
                && x->stocks_at_hand.find(k)->second < 0) {

                short_stocks += x->stocks_at_hand.find(k)->second;
            }
        }



//create price maps
//buy-limits
        map<double, vector<order>, greater<>> bid_orders = {};
        if (!this->bid.empty()) {
            if (this->bid.find(k) != this->bid.end()) {
                for (auto &x : this->bid.find(k)->second) {
                    double _bid_price = x.get_proposed_price();
                    _bid_price = (floor(_bid_price * 100)) / 100;

                    if (bid_orders.find(_bid_price) != bid_orders.end()) {
                        bid_orders.find(_bid_price)->second.push_back(x);
                    } else {
                        bid_orders.emplace(_bid_price, vector({x}));
                    }


                }
            }
        }



//sell-limits
        map<double, vector<order>> ask_orders = {};
        if (!this->ask.empty()) {
            if (this->ask.find(k) != this->ask.end()) {
                for (auto &x : this->ask.find(k)->second) {
                    double _ask_price = x.get_proposed_price();
                    _ask_price = (floor(_ask_price * 100)) / 100;

                    if (ask_orders.find(_ask_price) != ask_orders.end()) {
                        ask_orders.find(_ask_price)->second.push_back(x);
                    } else {
                        ask_orders.emplace(_ask_price, vector({x}));
                    }
                }
            }
        }



//market orders

//sum of orders
        double sum_ = 0;
        if (!this->market_orders.empty()) {
            if (this->market_orders.find(k) != this->market_orders.end()) {
                for (auto &x : this->market_orders.find(k)->second) {
                   if (x.get_status() ==1){
                       sum_ += abs(x.get_order_size());
                   }
                }
            }
        }

//probability of first execution
    map<double, int> probability_ = {};
        if (!this->market_orders.empty()) {
            if (this->market_orders.find(k) != this->market_orders.end()) {
                for (auto &x : this->market_orders.find(k)->second){
                    if(x.get_status() ==1){
                        probability_.emplace(abs(x.get_order_size())/sum_, x.get_id());
                }
            }
        }
    }

//cumulative probabilities
    map<int, double> cum_probability_ = {};
        double cum_sum = 0;
            for (auto &[p, id] : probability_){
                    cum_sum += p;
                        cum_probability_.emplace(id, cum_sum);
            }

//TODO pick executed transaction by cum_probability
        std::random_device rd;
        std::mt19937 e2(rd());
        std::uniform_real_distribution<> dist(0, 1);
        double random_number = dist(e2);

        int first_comer = 0;
        for (auto &[id, p] : cum_probability_){
            if (p>= random_number){
                first_comer = id;
                break;
                }
        }


        map<double, vector<order>> market_bids_asks = {};
        if (!this->market_orders.empty()) {
            if (this->market_orders.find(k) != this->market_orders.end()) {
                for (auto &x : this->market_orders.find(k)->second) {
                    double market_price = x.get_proposed_price();
                    market_price = (floor(market_price * 100)) / 100;
                    if (market_bids_asks.find(market_price) != market_bids_asks.end()) {
                        market_bids_asks.find(market_price)->second.push_back(x);
                    } else {
                        market_bids_asks.emplace(market_price, vector({x}));
                    }

                }
            }
        }



//match market orders
        for (auto &[price, market_orders]:market_bids_asks) {
            for (auto &x:market_orders) {
                if (x.get_id() == first_comer && x.get_status()==1) {
                    tuple<order, order> new_orders1 =
                            update_mm_orders(bid_orders, ask_orders, this->initial_invetory.find(k)->second,
                                             (this->stocks_at_hand.find(k)->second), stock.get_midprice(), 3,
                                             this->get_identifier(), k);
                    double kh = x.get_order_size();

                    add_mm_orders(new_orders1, bid_orders, ask_orders);
//                    std::cout << "mmsk0 " << get<1>(new_orders1).get_proposed_price() << std::endl;
                    //match market orders with limit orders
                    tuple<order, order>
                            matched_orders =
                            match_market_orders(x, *get_best_bid(bid_orders), *get_best_ask(ask_orders));


                    executed_orders_ = {};
                    executed_orders_ = record_executed_orders(matched_orders);


                    //update balances after each match
                    if (executed_orders_.find(this->get_identifier()) != executed_orders_.end()) {
                        balance_bd(executed_orders_.find(
                                this->get_identifier()
                        )->second);
                    }

                    for (auto &[i, v]:this->trading_institutions) {
                        bool trader_participated = executed_orders_.find(i) != executed_orders_.end();
                        if (trader_participated) {
                            auto exec_order = executed_orders_.find(i)->second;
                            v->balance_bd(exec_order);
                        }
                    }

                }
           }
    }

int active_mmbid_orders;
int active_mmask_orders;
//update market maker's orders
int active_market_orders =0;
        for (auto &[price, orders]:market_bids_asks) {
            for (auto &x:orders) {
                if(x.get_status() == 1){
                    active_market_orders++;
                }
            }
        }



do {
//    std::cout<<"passed by 0"<<std::endl;
    active_mmbid_orders =0;
    active_mmask_orders =0;

//update market maker's orders
    cancel_mm_orders(
            this->get_identifier(), bid_orders, ask_orders);

    tuple<order, order> new_orders =
            update_mm_orders(bid_orders, ask_orders, this->initial_invetory.find(k)->second,
                             (this->stocks_at_hand.find(k)->second), stock.get_midprice(), 3,
                             this->get_identifier(), k);

            add_mm_orders(new_orders, bid_orders, ask_orders);



//match limit orders after all market orders have been executed
    order *best_bid_order = get_best_bid(bid_orders);
    order *best_ask_order = get_best_ask(ask_orders);

//    std::cout<<" best_bid_price "<<best_bid_order->get_proposed_price()<<" size " << best_bid_order->get_order_size()<<std::endl;
//    std::cout<<" best_ask_price "<<best_ask_order->get_proposed_price()<<" size " << best_ask_order->get_order_size()<<std::endl;
//    std::cout<<" mmask "<<get<1>(new_orders).get_proposed_price()<<std::endl;


    if (best_bid_order->get_proposed_price() >= best_ask_order->get_proposed_price()) {
//        std::cout<<"passed by 2"<<std::endl;
        tuple<order, order>
                matched_limits2 =
                match_limit_orders(*get_best_bid(bid_orders), *get_best_ask(ask_orders));


        executed_orders_ = {};
        executed_orders_ = record_executed_orders(matched_limits2);

        //update balances after each match
        if (executed_orders_.find(this->get_identifier()) != executed_orders_.end()) {
            balance_bd(executed_orders_.find(
                    this->get_identifier()
            )->second);
        }

        for (auto &[i, v]:this->trading_institutions) {
            bool trader_participated = executed_orders_.find(i) != executed_orders_.end();
            if (trader_participated) {
                auto exec_order = executed_orders_.find(i)->second;
                v->balance_bd(exec_order);
            }
        }
    }

    for(auto &[price, vec]:bid_orders){
        for(auto &ord:vec){
            if (ord.get_id() == this->get_identifier()
                    && ord.get_status() ==1 ){
                active_mmbid_orders ++;
//                std::cout<<"passed by 3"<<std::endl;
                goto exit;
            }
        }
    }exit:;

    for(auto &[price, vec]:ask_orders){
        for(auto &ord:vec){
            if (ord.get_id() == this->get_identifier()
                && ord.get_status() ==1 ){
                active_mmask_orders ++;
//                std::cout<<"passed by 4"<<std::endl;
                goto exit_;
            }
        }
    }exit_:;
//    prevailing_ask = floor(get_best_ask(ask_orders)->get_proposed_price()*100)/100;
//    prevailing_bid = floor(get_best_bid(bid_orders)->get_proposed_price()*100)/100;
//std::cout<<" passed 8 "<<std::endl;
//    std::cout<<" mmbid "<<(active_mmbid_orders==0)<<std::endl;
//    std::cout<<" mmask "<<(active_mmask_orders==0)<<std::endl;
}while ((active_mmbid_orders==0) || (active_mmask_orders==0));
//
//else{
//    do {
//    active_mmbid_orders =0;
//    active_mmask_orders =0;
////update market maker's orders
//    cancel_mm_orders(
//            this->get_identifier(), bid_orders, ask_orders);
//
//    tuple<order, order> new_orders =
//            update_mm_orders(bid_orders, ask_orders, this->initial_invetory.find(k)->second,
//                             (this->stocks_at_hand.find(k)->second+short_stocks), stock.get_midprice(), 3,
//                             this->get_identifier(), k);
//
//    add_mm_orders(new_orders, bid_orders, ask_orders);
//
////    std::cout<<"mm_ask "<<get<1>(new_orders).get_proposed_price()<<" mm_bid "<< get<0>(new_orders).get_proposed_price()<<std::endl;
////match limit orders after all market orders have been executed
//    double best_bid_price = floor(get_best_bid(bid_orders)->get_proposed_price() * 100) / 100;
//    double best_ask_price = floor(get_best_ask(ask_orders)->get_proposed_price() * 100) / 100;
//
//    if (best_bid_price == best_ask_price) {
//        tuple<order, order>
//                matched_limits2 =
//                match_limit_orders(*get_best_bid(bid_orders), *get_best_ask(ask_orders));
//
//
//        executed_orders_ = {};
//        executed_orders_ = record_executed_orders(matched_limits2);
//
//        //update balances after each match
//        if (executed_orders_.find(this->get_identifier()) != executed_orders_.end()) {
//            balance_bd(executed_orders_.find(
//                    this->get_identifier()
//            )->second);
//        }
//
//        for (auto &[i, v]:this->trading_institutions) {
//            bool trader_participated = executed_orders_.find(i) != executed_orders_.end();
//            if (trader_participated) {
//                auto exec_order = executed_orders_.find(i)->second;
//                v->balance_bd(exec_order);
//            }
//        }
//    }
//
//    for(auto &[price, vec]:bid_orders){
//        for(auto &ord:vec){
//            if (ord.get_id() == this->get_identifier()
//                && ord.get_status() ==1 ){
//                active_mmbid_orders ++;
//                goto exitpoint;
//            }
//        }
//    }exitpoint:;
//
//    for(auto &[price, vec]:ask_orders){
//        for(auto &ord:vec){
//            if (ord.get_id() == this->get_identifier()
//                && ord.get_status() ==1 ){
//                active_mmask_orders ++;
//                goto exitpoint_;
//            }
//        }
//    }exitpoint_:;
////    prevailing_ask = floor(get_best_ask(ask_orders)->get_proposed_price()*100)/100;
////    prevailing_bid = floor(get_best_bid(bid_orders)->get_proposed_price()*100)/100;
//
//}while ((active_mmbid_orders==0) || (active_mmask_orders==0));
//}

        prevailing_ask = floor(get_best_ask(ask_orders)->get_proposed_price()*100)/100;
        prevailing_bid = floor(get_best_bid(bid_orders)->get_proposed_price()*100)/100;

//        std::cout<<"after mmUpdate: ask "<<prevailing_ask<<" "<<"bid "<<prevailing_bid<<"midprice "<<(prevailing_ask+prevailing_bid)/2 <<std::endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//inform traders about new prices
    for (auto &[i, v]: this->trading_institutions) {
        v->stocks_on_market.find(k)->second.set_price(
                int(this->clock->current_time()), prevailing_bid, prevailing_ask);

    }

        this->assets.find(k)->second.set_price(
                int(this->clock->current_time())
                    , prevailing_bid
                    , prevailing_ask);



    this->best_bid.find(k)->second = prevailing_bid;
    this->best_ask.find(k)->second = prevailing_ask;

    std::cout<<"ask "<<prevailing_ask<<" "<<"bid "<<prevailing_bid<<"midprice "<<(prevailing_ask+prevailing_bid)/2 <<std::endl;
}





//this->clock->tick();
std::cout<<"inner time "<< this->clock->current_time()<<std::endl;
std::cout<<"outer time "<< market_clock.current_time()<<std::endl;
}


