//
// Created by Cephas Svosve on 5/1/2022.
//

#include "price_setter.h"


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
        a.set_shares_outstanding(100000);
        this->assets.emplace(a.get_identifier(),a);
        }
    }
}

//this generates fundamental processes and associates them with each asset
void price_setter::generate_fundamentals(){
    auto earnings_processes  = generate(Earnings);
        auto dividends_processes = generate(Dividends);
            auto free_cash_flow      = generate(Free_Cash_Flows);
                int count = 0;

        for(auto &[i,j] : this->assets){
            j.set_earnings_process(earnings_processes.at(count));
                j.set_dividends_process(dividends_processes.at(count));
                    j.set_free_cash_flow_process(free_cash_flow.at(count));

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
        v.set_price(watch.current_time(),init_price,init_price);


        this->bid.emplace(k,vector<order>());
        this->ask.emplace(k,vector<order>());
        this->market_orders.emplace(k,vector<order>());

        this->best_bid.emplace(k,init_price);
        this->best_ask.emplace(k,init_price);
        this->mid_price.emplace(k,init_price);
    }
    this->get_inform(this->assets);
}

void price_setter::receive_orders(){
    for(auto &[k,v] : this->trading_institutions){
        auto order_ = v->invest();
        for(auto &[i,j] : order_){


           // std::cout<<"ordre "<<i<<" "<<x<<std::endl;

            if(j.get_order_type() != 0){//this checks whether the order is a limit(1) or market(0) order
                if(j.get_order_size() > 0){//this checks whether order is a bid(+) or ask(-)
                    stack_bid(j.get_ordered_asset(), j);
            } else{
                stack_ask(j.get_ordered_asset(),j);
                }
            }else{//if order is a market order we stack the market order vector
               stack_market_orders(j.get_ordered_asset(),j);

            }
        }
    }

}


void price_setter::stack_market_orders(int id, order order_){
    auto k = id;
        auto i = this->market_orders.find(k);
            if(i != this->market_orders.end()){
                //identify orders replaced and and cancel them
                for(auto &u :this->market_orders.find(k)->second){
                    if(order_.get_id() == u.get_id()
                       && (sgn(order_.get_order_size()) == sgn(u.get_order_size())))
                    {
                        u.set_order_size(order_.get_order_size(),order_.get_proposed_price());
                        u.set_id(order_.get_id());
                        goto exit_point;
                    }
                }
                this->market_orders.find(k)->second.push_back(order_);
            }else{
            vector<order> a;
        a.push_back(order_);
    this->market_orders.emplace(k,a);
    }
             exit_point:;
}




void price_setter::stack_bid(int id, order order_){
    auto k = id;
        auto  u= this->bid.find(k);
            if(u != this->bid.end()){
                //identify orders replaced and and cancel them
                for(auto &i :this->bid.find(k)->second){
                    if(order_.get_id() == i.get_id()
                       && sgn(order_.get_order_size()) == sgn(i.get_order_size()))
                    {
                        i.set_order_size(order_.get_order_size(),order_.get_proposed_price());
                        i.set_id(order_.get_id());
                        goto exit_point;
                    }
                }
                    this->bid.find(k)->second.push_back(order_);
                }else{
            vector<order> a;
         a.push_back(order_);
    this->bid.emplace(k,a);
    }
            exit_point:;
}




void price_setter::stack_ask(int id, order order_){
    auto k = id;
        auto c = this->ask.find(k);
            if(c != this->ask.end()){
                //identify orders replaced and and cancel them from the queue
                        for(auto &i :this->ask.find(k)->second){
                            if(sgn(order_.get_order_size()) == sgn(i.get_order_size())
                                && order_.get_id() == i.get_id())
                            {
                               i.set_order_size(order_.get_order_size(),order_.get_proposed_price());
                               i.set_id(order_.get_id());
                                goto exit_point;
                            }
                        }
                        this->ask.find(k)->second.push_back(order_);
                       }else{
            vector<order> a;
        a.push_back(order_);
    this->ask.emplace(k,a);
    }
            exit_point:;
}



//this method sorts compares orders by order prices
bool order_sort(order x, order y){
    return ((x.get_proposed_price()) < (y.get_proposed_price()));
}

//this method sorts compares orders by order sizes
bool sort_by_size(order x, order y){
    return ((x.get_order_size()) > (y.get_order_size()));
}


void price_setter::update_price(int k) {

    double bid_price;
    double ask_price;

    auto i = this->bid.find(k)->second.empty();
    auto j = this->ask.find(k)->second.empty();

    //if the bids are not empty, set the bid_price to the best bid proposed price
    //else set bid_price to the last traded price
    if (i != 1 ){
        bid_price = this->bid.find(k)->second.rbegin()->get_proposed_price();
    }else{
        bid_price = this->best_bid.find(k)->second;
    }

    //we do the same for ask_price
    if (j != 1 ){
        ask_price = this->ask.find(k)->second.begin()->get_proposed_price();
    }else{
        ask_price = this->best_ask.find(k)->second;
    }


    this->best_bid.find(k)->second = bid_price;
    this->best_ask.find(k)->second = ask_price;
    this->mid_price.find(k)->second = (bid_price + ask_price)/2.;

}


//this method assesses whether a vector still has both buy and sell orders
bool
price_setter::not_similar(vector<order> &x){
    bool status;
    int buy_count = 0;
    int sell_count = 0;
    for (auto &i : x) {
        if (i.get_order_size() > 0) {
            buy_count++;
        } else if (i.get_order_size() < 0){
            sell_count++; }
    }
    if(buy_count >0 && sell_count>0){
        status = true;
    } else{status=false;}

    return status;
}

bool
price_setter::isValid(int asset_id){
    int result=0;
    for(auto &[i,v]:this->assets){
        if(v.get_identifier() == asset_id){
            result ++;
        }
    }

    return (result > 0);
}

void
remove_empty_orders(int k
                        ,map<int,vector<order>> &bid_map
                            , map<int,vector<order>> &ask_map
                                , map<int,vector<order>> &market_order_map){

    if(!bid_map.find(k)->second.empty()) {
        std::sort(bid_map.find(k)->second.begin(), bid_map.find(k)->second.end(), sort_by_size);
        int size = bid_map.find(k)->second.size();
        for(int i = size; i>0; i--){
            if (abs(bid_map.find(k)->second.rbegin()->get_order_size()) == 0){
                bid_map.find(k)->second.pop_back();
            }else{
                break;
            }
        }
    }

    if(!ask_map.find(k)->second.empty()) {
        std::sort(ask_map.find(k)->second.begin(), ask_map.find(k)->second.end(), sort_by_size);
        int size = ask_map.find(k)->second.size();
        for(int i = size; i>0; i--){
            if (abs(ask_map.find(k)->second.rbegin()->get_order_size()) == 0){
                ask_map.find(k)->second.pop_back();
            }else{
                break;
            }
        }
    }

    if(!market_order_map.find(k)->second.empty()) {
        std::sort(market_order_map.find(k)->second.begin(), market_order_map.find(k)->second.end(), sort_by_size);
        int size = market_order_map.find(k)->second.size();
        for(int i = size; i>0; i--){
            if (abs(market_order_map.find(k)->second.rbegin()->get_order_size()) == 0){
                market_order_map.find(k)->second.pop_back();
            }else{
                break;
            }
        }
    }

}

void
price_setter::clear(){

//furthest range a market maker is allowed to charge away from the NBBO price
double regulation_range = 0.08;


map<int,vector<order>> executed_orders_;
//create a template for return result
    for(auto &[k,v] : this->trading_institutions) {
            executed_orders_.emplace(k, vector<order>());
    }


for (auto &[k, v]:assets) {

    //check if balance sheet is sufficient
//    if(this->wealth > 0){}
//    {
//        if (this->cash_at_hand < 0) {}
//        if (this->stocks_at_hand.find(k)->second < 0) {}
//    }else{
//
//    }


    // remove empty orders
    if(!this->bid.find(k)->second.empty()) {
        std::sort(bid.find(k)->second.begin(), bid.find(k)->second.end(), sort_by_size);
        int size = this->bid.find(k)->second.size();
        for(int i = size; i>0; i--){
            if (abs(this->bid.find(k)->second.rbegin()->get_order_size()) == 0){
                this->bid.find(k)->second.pop_back();
            }else{
                break;
            }
        }
    }

    if(!this->ask.find(k)->second.empty()) {
        std::sort(ask.find(k)->second.begin(), ask.find(k)->second.end(), sort_by_size);
        int size = this->ask.find(k)->second.size();
        for(int i = size; i>0; i--){
            if (abs(this->ask.find(k)->second.rbegin()->get_order_size()) == 0){
                this->ask.find(k)->second.pop_back();
            }else{
                break;
            }
        }
    }

    if(!this->market_orders.find(k)->second.empty()) {
        std::sort(market_orders.find(k)->second.begin(), market_orders.find(k)->second.end(), sort_by_size);
        int size = this->market_orders.find(k)->second.size();
        for(int i = size; i>0; i--){
            if (abs(this->market_orders.find(k)->second.rbegin()->get_order_size()) == 0){
                this->market_orders.find(k)->second.pop_back();
            }else{
                break;
            }
        }
    }



    // sort the orders in ascending order by order price;
    if(!this->bid.find(k)->second.empty()) {
        std::sort(bid.find(k)->second.begin(), bid.find(k)->second.end(), order_sort);
    }

    if(!this->ask.find(k)->second.empty()) {
        std::sort(ask.find(k)->second.begin(), ask.find(k)->second.end(), order_sort);
    }
    if(!this->market_orders.find(k)->second.empty()) {
        std::sort(market_orders.find(k)->second.begin(), market_orders.find(k)->second.end(),
                  order_sort);
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{////calculate weighted average price of limit orders
//
//
//    //std::cout<<"*******stocks at hand"<<k<<" is " << this->stocks_at_hand.find(k)->second <<std::endl;
//    //std::cout<<"*******cash at hand"<<" is " << this->cash_at_hand <<std::endl;
//    double demanded_stocks = 0.;
//    double count = 0.;
//    for (auto &i : this->bid.find(k)->second) {
//        demanded_stocks += i.get_proposed_price() * i.get_order_size();
//        count += i.get_order_size();
//    }
////market maker's ask/offer price (price at which market maker is willing to sell stocks) is determined as a weighted average of
//// the traders' demands (stocks that traders are willing to buy)
//    double MM_ASK = demanded_stocks / count;
//
//
//    double supplied_stocks = 0.;
//    double total = 0.;
//    for (auto &i : this->ask.find(k)->second) {
//        supplied_stocks += i.get_proposed_price() * i.get_order_size();
//        total += i.get_order_size();
//    }
//
////market maker's bid price (price at which market maker is willing to buy stocks) is determined as a weighted average of
//// the traders' supply (stocks that traders are willing to sell)
//    double MM_BID = supplied_stocks / total;
//
//
////given all orders fall on one side, we assume a bid-ask spread of 0
//    if (isnan(MM_BID) && isnan(MM_ASK)){
//        MM_BID = this->last_bid.find(k)->second;
//        MM_ASK = this->last_ask.find(k)->second;
//    }else if ( isnan(MM_ASK)){
//        MM_ASK = MM_BID;
//    } else if (isnan(MM_BID)){
//        MM_BID = MM_ASK;
//    }
//
//    //std::cout<<"bidiing" << MM_BID<<std::endl;
//    this->last_bid.find(k)->second = MM_BID;
//    this->last_ask.find(k)->second = MM_ASK;
//    this->last_price.find(k)->second = (MM_BID + MM_ASK)/2.;
//
//
//    //clear bids
//    while (!this->bid.find(k)->second.empty() &&
//            this->bid.find(k)->second.rbegin()->get_proposed_price() >= MM_ASK) {
//
//
//        auto trader_demand = this->bid.find(k)->second.rbegin();
//        auto uncleared = this->stocks_at_hand.find(k)->second - (trader_demand->
//                get_order_size());
//
//        //check whether market maker is able to fill the order
//        if (uncleared > 0.) {
//
//            //calculate executed bid
//            auto size = this->bid.find(k)->second.size() - 1;
//            auto exec_bid = this->bid.find(k)->second[size];
//            auto exec_bid_price = this->bid.find(k)->second.rbegin()->get_proposed_price();
//            auto bid_trader_id = exec_bid.get_id();
//            exec_bid.set_order_size(trader_demand->get_order_size(), exec_bid_price);
//            executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
//
//
//
//
//            //update stocks remaining
//            this->stocks_at_hand.find(k)->second -= (*trader_demand).get_order_size();
//            this->cash_at_hand += trader_demand->get_order_size() * trader_demand->get_proposed_price();
//
//
//
//            //remove filled order
//            this->bid.find(k)->second.pop_back();
//
//
//        } else if (uncleared < 0.) {
//
//            //calculate executed bid
//            auto size = this->bid.find(k)->second.size() - 1;
//            auto exec_bid = this->bid.find(k)->second[size];
//            auto exec_bid_price = this->bid.find(k)->second.rbegin()->get_proposed_price();
//            auto bid_trader_id = exec_bid.get_id();
//            exec_bid.set_order_size(this->stocks_at_hand.find(k)->second, exec_bid_price);
//
//            executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
//
//
//
//
//
//            //calculate remaining stocks held
//            this->cash_at_hand += this->stocks_at_hand.find(k)->second * trader_demand->get_proposed_price();
//            this->stocks_at_hand.find(k)->second = 0;
//
//
//            //update the bid order remaining in the bid_order vector after execution of trade
//            this->bid.find(k)->
//                    second.rbegin()->
//                    set_order_size(-uncleared, exec_bid_price);
//
//            //std::cout<<"passed 1"<<this->stocks_at_hand.find(k)->second<<std::endl;
//        } else {
//
//            //calculate executed bid
//            auto size = this->bid.find(k)->second.size() - 1;
//            auto exec_bid = this->bid.find(k)->second[size];
//            auto exec_bid_price = this->bid.find(k)->second.rbegin()->get_proposed_price();
//            auto bid_trader_id = exec_bid.get_id();
//            exec_bid.set_order_size(this->stocks_at_hand.find(k)->second, exec_bid_price);
//
//            executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
//
//
//
//
//
//            //calculate executed ask
//            this->cash_at_hand += this->stocks_at_hand.find(k)->second* trader_demand->get_proposed_price();
//            this->stocks_at_hand.find(k)->second = 0;
//
//
//
//            //remove filled order
//            this->bid.find(k)->second.pop_back();
//            //std::cout<<"passed 2"<<this->stocks_at_hand.find(k)->second<<std::endl;
//
//        }
//    }
//
//
//
//
////clear asks
//    while (!this->ask.find(k)->second.empty() &&
//           this->ask.find(k)->second.begin()->get_proposed_price() <= MM_BID) {
//
//        auto trader_supply = this->ask.find(k)->second[0];
//
//
//        auto y = trader_supply.
//                get_order_size();
//
//        auto z = trader_supply.get_proposed_price();
//
//        auto uncleared = this->cash_at_hand + (trader_supply.
//                get_order_size() * trader_supply.get_proposed_price());
//
//        //check whether market maker is able to fill the order
//        if (uncleared > 0.) {
//
//            //calculate executed ask
//            auto exec_ask = this->ask.find(k)->second[0];
//            auto exec_ask_price = trader_supply.get_proposed_price();
//            auto ask_trader_id = exec_ask.get_id();
//            exec_ask.set_order_size(trader_supply.get_order_size(), exec_ask_price);
//            executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
//
//
//
//            //update cash at hand
//            this->cash_at_hand += trader_supply.get_order_size() * trader_supply.get_proposed_price();
//            this->stocks_at_hand.find(k)->second -= trader_supply.get_order_size();
//
//
//
//            //remove the filled order
//            this->ask.find(k)->second.erase(this->ask.find(k)->second.begin());
//            //std::cout<<"passed 3"<<this->stocks_at_hand.find(k)->second<<std::endl;
//
//        } else if (uncleared < 0.) {
//
//            //calculate executed ask
//            auto exec_ask = this->ask.find(k)->second[0];
//            auto exec_ask_price = trader_supply.get_proposed_price();
//            auto ask_trader_id = exec_ask.get_id();
//            exec_ask.set_order_size(this->cash_at_hand / exec_ask_price, exec_ask_price);
//            executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
//
//
//
//
//            //calculate executed ask
//            this->stocks_at_hand.find(k)->second += this->cash_at_hand / exec_ask_price;
//            this->cash_at_hand = 0;
//
//
//
//            //update the ask order remaining in the bid_order vector after execution of trade
//            this->ask.find(k)->
//                    second.begin()->
//                    set_order_size(uncleared / exec_ask_price, exec_ask_price);
//            //std::cout<<"passed 4"<<this->stocks_at_hand.find(k)->second<<std::endl;
//
//        } else {
//
//            //calculate executed ask
//            auto exec_ask = this->ask.find(k)->second[0];
//            auto exec_ask_price = trader_supply.get_proposed_price();
//            auto ask_trader_id = exec_ask.get_id();
//            exec_ask.set_order_size(trader_supply.get_order_size(), exec_ask_price);
//            executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
//
//
//
//
//
//            //update cash
//            this->cash_at_hand = 0;
//            this->stocks_at_hand.find(k)->second -= trader_supply.get_order_size();
//
//
//
//
//            //remove filled order
//            this->ask.find(k)->second.erase(this->ask.find(k)->second.begin());
//            //std::cout<<"passed 5"<<this->stocks_at_hand.find(k)->second<<std::endl;
//        }
//    }
//
//
//    //clear market orders
//    while (!this->market_orders.find(k)->second.empty() &&
//                isValid(this->market_orders.find(k)->second.rbegin()->get_ordered_asset())) {
//        double a = this->stocks_at_hand.find(k)->second;
//        double time = watch.current_time();
//      //std::cout<<"time "<<watch.current_time()<<"stocks "<<a<<std::endl;
//        auto market_order = this->market_orders.find(k)->second.rbegin();
//        int market_order2 = market_order->get_ordered_asset();
//        double market_order3 = this->market_orders.find(k)->second.size();
//        if(abs(market_order->get_order_size()) == 0){
//            this->market_orders.find(k)->second.pop_back();
//        }
//
//        //sell-side market order
//        if (market_order->get_order_size() < 0.){
//
//
//            auto trader_supply = market_order;
//            auto uncleared = this->cash_at_hand + (trader_supply->
//                    get_order_size() * MM_BID);
//            //std::cout<<"cash "<<this->cash_at_hand <<std::endl;
//            //std::cout<<"order size "<<trader_supply->
//                    get_order_size()<<std::endl;
//            //std::cout<<"MMBID "<<MM_BID<<std::endl;
//            //check whether market maker is able to fill the order
//            if (uncleared > 0.) {
//
//                //calculate executed ask
//                auto size = this->market_orders.find(k)->second.size() - 1;
//                auto exec_ask = this->market_orders.find(k)->second[size];
//                auto exec_ask_price = MM_BID;
//                auto ask_trader_id = exec_ask.get_id();
//                exec_ask.set_order_size(trader_supply->get_order_size(), exec_ask_price);
//                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
//
//
//                //update cash at hand
//                this->cash_at_hand += trader_supply->get_order_size() * MM_BID;
//                this->stocks_at_hand.find(k)->second -= trader_supply->get_order_size();
//
//
//                //remove the filled order
//                this->market_orders.find(k)->second.pop_back();
//                //std::cout<<"passed 6"<<this->stocks_at_hand.find(k)->second<<std::endl;
//
//            } else if (uncleared < 0.) {
//
//                //calculate executed ask
//                auto size = this->market_orders.find(k)->second.size() - 1;
//                auto exec_ask = this->market_orders.find(k)->second[size];
//                auto exec_ask_price = MM_BID;
//                auto ask_trader_id = exec_ask.get_id();
//                exec_ask.set_order_size(this->cash_at_hand/ exec_ask_price, exec_ask_price);
//                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
//
//
//
//
//                //calculate executed ask
//                this->stocks_at_hand.find(k)->second += this->cash_at_hand / exec_ask_price;
//                this->cash_at_hand = 0;
//
//
//                //update the ask order remaining in the bid_order vector after execution of trade
//                this->market_orders.find(k)->
//                        second.rbegin()->
//                        set_order_size(uncleared / exec_ask_price, exec_ask_price);
//                //std::cout<<"time "<<watch.current_time()<<std::endl;
//                //std::cout<<"passed 7"<<this->stocks_at_hand.find(k)->second<<std::endl;
//
//            } else {
//
//                //calculate executed ask
//                auto size = this->market_orders.find(k)->second.size() - 1;
//                auto exec_ask = this->market_orders.find(k)->second[size];
//                auto exec_ask_price = MM_BID;
//                auto ask_trader_id = exec_ask.get_id();
//                exec_ask.set_order_size(this->cash_at_hand / exec_ask_price, exec_ask_price);
//                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
//
//
//
//
//
//                //update cash
//                this->stocks_at_hand.find(k)->second += this->cash_at_hand / exec_ask_price;
//                this->cash_at_hand = 0;
//
//
//
//                //remove filled order
//                this->market_orders.find(k)->second.pop_back();
//                //std::cout<<"passed 8"<<this->stocks_at_hand.find(k)->second<<std::endl;
//            }
//        } else
//            //buy-side market order
//        if (market_order->get_order_size() > 0.) {
//
//
//            auto trader_demand = this->market_orders.find(k)->second.rbegin();
//            auto uncleared = this->stocks_at_hand.find(k)->second - (trader_demand->
//                    get_order_size());
//            //std::cout<<"stocks at hand" << this->stocks_at_hand.find(k)->second <<std::endl;
//            //std::cout<<"cash at hand" << this->cash_at_hand <<std::endl;
//            //std::cout<<"uncleared" << uncleared <<std::endl;
//            //check whether market maker is able to fill the order
//            if (uncleared > 0.) {
//
//                //calculate executed bid
//                auto size = this->market_orders.find(k)->second.size() - 1;
//                auto exec_bid = this->market_orders.find(k)->second[size];
//                auto exec_bid_price = MM_ASK;
//                auto bid_trader_id = exec_bid.get_id();
//                exec_bid.set_order_size(trader_demand->get_order_size(), exec_bid_price);
//                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
//
//
//
//
//
//                //update stocks remaining
//                this->stocks_at_hand.find(k)->second -= trader_demand->get_order_size();
//                this->cash_at_hand += trader_demand->get_order_size()*trader_demand->get_proposed_price();
//
//
//
//                //remove filled order
//                this->market_orders.find(k)->second.pop_back();
//
//                //std::cout<<"passed 9 *"<<this->stocks_at_hand.find(k)->second<<std::endl;
//
//            } else if (uncleared < 0.) {
//
//                //std::cout<<"time "<<watch.current_time()<<std::endl;
//                auto od = this->market_orders.find(k)->second.rbegin();
//                //std::cout<<"id  *"<<this->market_orders.find(k)->second.rbegin()->get_id()<<std::endl;
//
//
//                //calculate executed bid
//                auto size = this->market_orders.find(k)->second.size() - 1;
//                auto exec_bid = this->market_orders.find(k)->second[size];
//                auto exec_bid_price = MM_ASK;
//                auto bid_trader_id = exec_bid.get_id();
//                exec_bid.set_order_size(this->stocks_at_hand.find(k)->second, exec_bid_price);
//
//                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
//
//
//                //calculate remaining stocks held
//                this->cash_at_hand += this->stocks_at_hand.find(k)->second * trader_demand->get_proposed_price();
//                this->stocks_at_hand.find(k)->second = 0;
//
//
//                //update the bid order remaining in the bid_order vector after execution of trade
//                this->market_orders.find(k)->
//                        second.rbegin()->
//                        set_order_size(-uncleared, exec_bid_price);
//
//                //std::cout<<"passed 10"<<this->stocks_at_hand.find(k)->second<<std::endl;
//            } else {
//
//                //calculate executed bid
//                auto size = this->market_orders.find(k)->second.size() - 1;
//                auto exec_bid = this->market_orders.find(k)->second[size];
//                auto exec_bid_price = MM_ASK;
//                auto bid_trader_id = exec_bid.get_id();
//                exec_bid.set_order_size(this->stocks_at_hand.find(k)->second, exec_bid_price);
//
//                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
//
//
//
//
//
//                //calculate executed ask
//                this->cash_at_hand += this->stocks_at_hand.find(k)->second * trader_demand->get_proposed_price();
//                this->stocks_at_hand.find(k)->second = 0;
//
//
//
//                //remove filled order
//                this->market_orders.find(k)->second.pop_back();
////std::cout<<"passed 11"<<this->stocks_at_hand.find(k)->second<<std::endl;
//            }
//        }
//
//    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//clear limit orders when there is an overlap of orders
    //(1) check whether best bid > best ask and clear otherwise proceed to (2)
    while (!this->bid.find(k)->second.empty() && !this->ask.find(k)->second.empty() &&
           this->bid.find(k)->second.back().get_proposed_price() >=
           this->ask.find(k)->second.front().get_proposed_price()) {


        //..........(1.2)
        auto BID = this->bid.find(k)->second.back().get_order_size();
        auto ASK = this->ask.find(k)->second.front().get_order_size();
        auto uncleared = BID + ASK;




//situation (a).... bid order is partially filled where as ask order if fully filled;
// we adjust the bid order sizes and delete the ask order
        if (uncleared > 0.) {
            //calculate executed bid;
            auto exec_bid = this->bid.find(k)->second.back();
            auto exec_bid_price = this->bid.find(k)->second.back().get_proposed_price();
            auto bid_trader_id = exec_bid.get_id();
            exec_bid.set_order_size(BID - uncleared, exec_bid_price);
            executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);

//calculate executed ask
            auto exec_ask = this->ask.find(k)->second.front();
            auto exec_ask_price = this->ask.find(k)->second.front().get_proposed_price();
            auto ask_trader_id = exec_ask.get_id();
            exec_ask.set_order_size(-(BID - uncleared), exec_ask_price);
            executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);


//adjust marketmaker's inventory

            auto spread = exec_bid_price - exec_ask_price;
            this->cash_at_hand += (BID - uncleared) * spread;


//update the bid order remaining in the bid_order vector after execution of trade
            this->bid.find(k)->
                    second.back().
                    set_order_size(uncleared, this->bid.find(k)->
                    second.back().
                    get_proposed_price());


//remove the fully filled ask order
            this->ask.find(k)->
                    second.erase(this->ask.find(k)->
                    second.begin());
        }



//situation (b).... bid order is partially filled where as ask order if fully filled;
// we adjust the bid order sizes and delete the ask order

        else if (uncleared < 0.) {
            //calculate executed bid
            auto exec_bid = this->bid.find(k)->second.back();
            auto exec_bid_price = this->bid.find(k)->second.back().get_proposed_price();
            auto bid_trader_id = exec_bid.get_id();
            exec_bid.set_order_size(BID, exec_bid_price);
            executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);


            //calculate executed ask
            auto exec_ask = this->ask.find(k)->second.front();
            auto exec_ask_price = this->ask.find(k)->second.front().get_proposed_price();
            auto ask_trader_id = exec_ask.get_id();
            exec_ask.set_order_size(-(BID), exec_ask_price);
            executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);



//adjust marketmaker's inventory

            auto spread = exec_bid_price - exec_ask_price;
            this->cash_at_hand += (BID - uncleared) * spread;



//update the ask orders remaining in the order vectors after execution of trade
            this->ask.find(k)->
                    second.front().
                    set_order_size(uncleared, this->ask.find(k)->
                    second.front().
                    get_proposed_price());


//remove the fully filled bid order
            this->bid.find(k)->second.pop_back();
        }



//situation (c).... both bid and ask orders are fully filled;
// we remove both orders from the order vectors
        else {

            //calculate executed bid
            auto exec_bid = this->bid.find(k)->second.back();
            auto exec_bid_price = this->bid.find(k)->second.back().get_proposed_price();
            auto bid_trader_id = exec_bid.get_id();
            exec_bid.set_order_size(BID, exec_bid_price);
            executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);


            //calculate executed ask
            auto exec_ask = this->ask.find(k)->second.front();
            auto exec_ask_price = this->ask.find(k)->second.front().get_proposed_price();
            auto ask_trader_id = exec_ask.get_id();
            exec_ask.set_order_size(-(BID), exec_ask_price);
            executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);


//adjust marketmaker's inventory
            auto spread = exec_bid_price - exec_ask_price;
            this->cash_at_hand += (BID - uncleared) * spread;


            //remove the fully filled orders
            this->ask.find(k)->second.erase(this->ask.find(k)->second.begin());
            this->bid.find(k)->second.pop_back();
        }
        this->update_price(k);
        //std::cout<<"process 1"<< std::endl;
    }



//situation (d.1.0).... no more matching limit orders but market orders need to be cleared
//(2) if best bid < best ask move on to clear market orders
//in the case of a bid market order that executes against an ask limit order


    while (!this->market_orders.find(k)->second.empty()) {

        //check whether this is a bid market order
        if (this->market_orders.find(k)->second.back().get_order_size() > 0.
            && !this->ask.find(k)->second.empty()) {


            //follow steps in .....(1.2)
            auto MKT_BID = this->market_orders.find(k)->second.back().get_order_size();
            auto ASK = this->ask.find(k)->second.front().get_order_size();
            auto uncleared = MKT_BID + ASK;



//situation (d.1.1).... market order partially filled, ask limit fully filled
//if order is partially filled we adjust the order sizes
//proceed as in the previous case in 1.2
            if (uncleared > 0.) {

                //calculate executed mkt_bid
                auto exec_mktbid = this->market_orders.find(k)->second.back();
                auto exec_mktbid_price = this->market_orders.find(k)->second.back().get_proposed_price();
                auto mkt_bid_trader_id = exec_mktbid.get_id();
                exec_mktbid.set_order_size(MKT_BID - uncleared, exec_mktbid_price);
                executed_orders_.find(mkt_bid_trader_id)->second.push_back(exec_mktbid);


                //calculate executed ask
                auto exec_ask = this->ask.find(k)->second.front();
                auto exec_ask_price = this->ask.find(k)->second.front().get_proposed_price();
                auto ask_trader_id = exec_ask.get_id();
                exec_ask.set_order_size(-(MKT_BID - uncleared), exec_ask_price);
                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);



//update the mkt orders remaining in the order vectors after execution of trade
                this->market_orders.find(k)->
                        second.back().
                        set_order_size(uncleared, this->market_orders.find(k)->
                        second.back().
                        get_proposed_price());


//remove the fully filled ask orders
                this->ask.find(k)->second.erase(this->ask.find(k)->second.begin());
            }





//situation (e).... market order fully filled, ask limit partially filled
            else if (uncleared < 0.) {

                //calculate executed mkt_bid
                auto exec_mktbid = this->market_orders.find(k)->second.back();
                auto exec_mktbid_price = this->market_orders.find(k)->second.back().get_proposed_price();
                auto mkt_bid_trader_id = exec_mktbid.get_id();
                exec_mktbid.set_order_size(MKT_BID, exec_mktbid_price);
                executed_orders_.find(mkt_bid_trader_id)->second.push_back(exec_mktbid);



                //calculate executed ask
                auto exec_ask = this->ask.find(k)->second.front();
                auto exec_ask_price = this->ask.find(k)->second.front().get_proposed_price();
                auto ask_trader_id = exec_ask.get_id();
                exec_ask.set_order_size(-(MKT_BID), exec_ask_price);
                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);




//update the ask orders remaining in the order vectors after execution of trade
                this->ask.find(k)->
                        second.front().
                        set_order_size(uncleared, this->ask.find(k)->
                        second.front().
                        get_proposed_price());


//remove the fully filled market orders
                this->market_orders.find(k)->second.pop_back();
            }


//situation (e).... both market and ask order are fully filled
            else {


                //calculate executed bid
                auto exec_mktbid = this->market_orders.find(k)->second.back();
                auto exec_mktbid_price = this->market_orders.find(k)->second.back().get_proposed_price();
                auto mktbid_trader_id = exec_mktbid.get_id();
                exec_mktbid.set_order_size(MKT_BID, exec_mktbid_price);
                executed_orders_.find(mktbid_trader_id)->second.push_back(exec_mktbid);




                //calculate executed ask
                auto exec_ask = this->ask.find(k)->second.front();
                auto exec_ask_price = this->ask.find(k)->second.front().get_proposed_price();
                auto ask_trader_id = exec_ask.get_id();
                exec_ask.set_order_size(-(MKT_BID), exec_ask_price);
                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);




//remove the fully filled ask orders
                this->ask.find(k)->
                        second.erase(this->ask.find(k)->
                        second.begin());

//remove the fully filled market orders
                this->market_orders.find(k)->second.pop_back();
            }
        }


//situation (f.1.0)....we have an ask market order in place of a bid market order
//in the case of an ask market order
        else if (this->market_orders.find(k)->second.back().get_order_size() < 0. &&
                 !this->bid.find(k)->second.empty()) {


            auto MKT_ASK = this->market_orders.find(k)->second.back().get_order_size();
            auto BID = this->bid.find(k)->second.back().get_order_size();
            auto uncleared = MKT_ASK + BID;


//if order is partially filled we adjust the order sizes
            if (uncleared < 0.) {

                //calculate executed mkt_ask
                auto exec_mktask = this->market_orders.find(k)->second.back();
                auto exec_mktask_price = this->market_orders.find(k)->second.back().get_proposed_price();
                auto mkt_ask_trader_id = exec_mktask.get_id();
                exec_mktask.set_order_size(MKT_ASK - uncleared, exec_mktask_price);
                executed_orders_.find(mkt_ask_trader_id)->second.push_back(exec_mktask);


                //calculate executed bid
                auto exec_bid = this->bid.find(k)->second.back();
                auto exec_bid_price = this->bid.find(k)->second.back().get_proposed_price();
                auto bid_trader_id = exec_bid.get_id();
                exec_bid.set_order_size(-(MKT_ASK - uncleared), exec_bid_price);
                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);


//update the partly filled market ask order
                this->market_orders.find(k)->
                        second.back().
                        set_order_size(uncleared, this->market_orders.find(k)->
                        second.back().
                        get_proposed_price());


//remove the fully filled bid order
                this->bid.find(k)->second.pop_back();
            }

//if order is fully filled we adjust the order sizes
            else if (uncleared > 0.) {

                //calculate executed mkt_ask
                auto exec_mktask = this->market_orders.find(k)->second.back();
                auto exec_mktask_price = this->market_orders.find(k)->second.back().get_proposed_price();
                auto mkt_ask_trader_id = exec_mktask.get_id();
                exec_mktask.set_order_size(MKT_ASK, exec_mktask_price);
                executed_orders_.find(mkt_ask_trader_id)->second.push_back(exec_mktask);


                //calculate executed bid
                auto exec_bid = this->bid.find(k)->second.back();
                auto exec_bid_price = this->bid.find(k)->second.back().get_proposed_price();
                auto bid_trader_id = exec_bid.get_id();
                exec_bid.set_order_size(-(MKT_ASK), exec_bid_price);
                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);



//update the partly filled bid order
                this->bid.find(k)->
                        second.back().
                        set_order_size(uncleared, this->bid.find(k)->
                        second.back().get_proposed_price());

//remove the fully filled market_ask order
                this->market_orders.find(k)->second.pop_back();
            }



//both market and bid order are fully filled
            else {

                //calculate executed mkt_ask
                auto exec_mktask = this->market_orders.find(k)->second.back();
                auto exec_mktask_price = this->market_orders.find(k)->second.back().get_proposed_price();
                auto mkt_ask_trader_id = exec_mktask.get_id();
                exec_mktask.set_order_size(MKT_ASK, exec_mktask_price);
                executed_orders_.find(mkt_ask_trader_id)->second.push_back(exec_mktask);




                //calculate executed bid
                auto exec_bid = this->bid.find(k)->second.back();
                auto exec_bid_price = this->bid.find(k)->second.back().get_proposed_price();
                auto bid_trader_id = exec_bid.get_id();
                exec_bid.set_order_size(-(MKT_ASK), exec_bid_price);
                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);


//remove the fully filled orders
                this->bid.find(k)->second.pop_back();
                this->market_orders.find(k)->second.pop_back();
            }
        }

//when only market orders are available with no limit orders to fill them,
//they are matched to corresponding market orders

        else if (this->market_orders.find(k)->second.back().get_order_size() < 0. &&
                 this->bid.find(k)->second.empty()) {



                double current_order = this->market_orders.find(k)->second.back().get_order_size();


                //calculate executed mkt_ask
                auto exec_mktask = this->market_orders.find(k)->second.back();
                auto exec_mktask_price = this->best_bid.find(k)->second*(1-regulation_range);
                auto mkt_ask_trader_id = exec_mktask.get_id();
                exec_mktask.set_order_size(current_order, exec_mktask_price);
                executed_orders_.find(mkt_ask_trader_id)->second.push_back(exec_mktask);


                //update market maker's inventory
                this->stocks_at_hand.find(k)->second += current_order;
                this->cash_at_hand -= current_order * exec_mktask_price;

                //remove the filled order from the list
                this->market_orders.find(k)->second.pop_back();

            }


        else if (this->market_orders.find(k)->second.back().get_order_size() > 0. &&
                 this->ask.find(k)->second.empty()){


                double current_cash = this->market_orders.find(k)->second.back().get_order_size() *
                                        this->market_orders.find(k)->second.back().get_proposed_price();


                //calculate executed mkt_ask
                auto exec_mktbid = this->market_orders.find(k)->second.back();
                auto exec_mktbid_price = this->best_ask.find(k)->second*(1+regulation_range);
                auto mkt_bid_trader_id = exec_mktbid.get_id();
                exec_mktbid.set_order_size(current_cash/exec_mktbid_price, exec_mktbid_price);
                executed_orders_.find(mkt_bid_trader_id)->second.push_back(exec_mktbid);


                //update market maker's inventory
                this->cash_at_hand += current_cash;
                this->stocks_at_hand.find(k)->second -= current_cash/exec_mktbid_price;

                //remove the filled order from the list
                this->market_orders.find(k)->second.pop_back();



        } else {
            //remove market orders with order size = 0
            std::cout<<"Time: "<< this->clock.current_time()<<" "<<
            "order id: "<<this->market_orders.find(k)->second.rbegin()->get_ordered_asset()<<" "<<
            "size: "<<this->market_orders.find(k)->second.rbegin()->get_order_size()<<" "<<
            "trader: "<<this->market_orders.find(k)->second.rbegin()->get_id()<<" "<<
            "price: "<<this->market_orders.find(k)->second.rbegin()->get_proposed_price()<<"**********"<<std::endl;
            this->market_orders.find(k)->second.pop_back();
        }
this->update_price(k);

        //std::cout<<"process 2"<< std::endl;
        //std::cout<<"size"<< this->market_orders.find(k)->second.begin()->get_order_size()<<std::endl;

    }





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//inform traders about new prices
    for (auto &[i, v]: this->trading_institutions) {
        v->stocks_on_market.find(k)->second.set_price(
                watch.current_time()
                    , this->best_bid.find(k)->second
                        , this->best_ask.find(k)->second);



        this->assets.find(k)->second.set_price(
                watch.current_time()
                    , this->best_bid.find(k)->second
                        , this->best_ask.find(k)->second);

    }
    this->balance();
}


this->clock.tick();
this->balance_cf(this->clock.current_time());
//std::cout<<"time "<<this->clock.current_time()<<std::endl;

    for(auto &[i,v]:this->trading_institutions){
//        //std::cout<<"div trader "<<v->stocks_on_market.find(1)->second.get_dividend(0)<<std::endl;
    auto exec_order = executed_orders_.find(i)->second;
        v->balance_bd(this->clock.current_time(),exec_order);
        v->clock.tick();
    }



}
