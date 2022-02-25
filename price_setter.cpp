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
    for(auto &[k,v] : this->trading_institutions) {
        if (v->wealth > 0) {
            auto order_ = v->invest();
            for (auto &[i, j] : order_) {
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
bool order_sort(order x, order y){
    return ((x.get_proposed_price()) < (y.get_proposed_price()));
}

//this method sorts orders by order sizes in descending order
bool sort_by_size(order x, order y){
    return ((x.get_order_size()) > (y.get_order_size()));
}


//this method sorts orders by order id in decreasing order in descending order
bool sort_by_id(order x, order y){
    return ((x.get_id()) > (y.get_id()));
}


void price_setter::update_price(int k, double reg) {

    double bid_price;
    double ask_price;

    bool bid_empty = this->bid.find(k)->second.empty();
    bool ask_empty = this->ask.find(k)->second.empty();

    //if the bids are not empty, set the bid_price to the best bid proposed price
    //else set bid_price to the last traded price
    if (!bid_empty ){
        bid_price = this->bid.at(k).rbegin()->get_proposed_price();
    }else{
        bid_price = this->best_bid.at(k)*(1+reg);
    }

    //we do the same for ask_price
    if (!ask_empty){
        ask_price = this->ask.at(k).begin()->get_proposed_price();
    }else{
        ask_price = this->best_ask.at(k)*(1-reg);
    }


    this->best_bid.at(k) = bid_price;
    this->best_ask.at(k) = ask_price;
    this->mid_price.at(k)= (bid_price + ask_price)/2.;

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
price_setter::clear(){


//Params definition
map<int,order> result_;
result_.clear();
auto quotes = this->stocks_on_market;
double regulation_range_ = 8;
auto t = clock.current_time();

//compute the ratio of each stock's signal relative to sum described in (1)
        for(auto &[s, v] : quotes){

            auto bid_ =this->best_bid.find(s)->second;
            auto ask_ =this->best_ask.find(s)->second;
if(!this->bid.empty()){
    if(this->bid.find(s)!= this->bid.end()){
    sort(this->bid.find(s)->second.begin(),this->bid.find(s)->second.end(),order_sort);
    bid_ = this->bid.find(s)->second.rbegin()->get_proposed_price();
    }
}

if(!this->ask.empty()) {
    if (this->ask.find(s) != this->ask.end()) {
        sort(this->ask.find(s)->second.begin(), this->ask.find(s)->second.end(), order_sort);
        ask_ = this->ask.find(s)->second.begin()->get_proposed_price();
    }
}









//investment in the subject stock

            order buy_limit;
            order sell_limit;
            auto order_num = concat(t, s);



//place sell-limits
            for(int i=0; i<floor(regulation_range_*ask_); i++) {


                order_num = concat(s, t);
                order_num = concat(-order_num, 100+i);
                double N = this->stocks_at_hand.find(s)->second;
                double slots = regulation_range_ * ask_;
                double quoted_price_ = ask_ + (i * 0.01);//move price by one ticker
                auto demand = -N / (slots);
                sell_limit.set_id(this->get_identifier());
                sell_limit.set_order_type(order::limit);
                sell_limit.set_ordered_asset(s);
                sell_limit.set_order_size(demand, quoted_price_);
                sell_limit.set_status(order::active);
                result_.emplace(order_num, sell_limit);

            }


//place buy-limits
//first calculate how much money is to be allocated t make a market on each stock
            double tot_price = 0;
            for(auto &[x,v] : this->stocks_at_hand){
                tot_price += (get<1>(this->stocks_on_market.find(x)->second.get_price())
                              + get<2>(this->stocks_on_market.find(x)->second.get_price()))/2.;
            }


            for(int ki=0; ki<floor(regulation_range_*bid_); ki++) {

//compute cash proportion for each particular stock
                order_num = concat(s, t);
                order_num = concat(order_num, 100+ki);
                double this_price = (get<1>(v.get_price())+get<2>(v.get_price()))/2.;
                double k_proportion = this_price/tot_price;
                double cash = k_proportion * this->cash_at_hand;
                double price = bid_ - (ki*0.01);
                double slots_ = regulation_range_ * bid_;
                auto demand = cash / (slots_*price);


//move price by one ticker
                buy_limit.set_id(this->get_identifier());
                buy_limit.set_order_type(order::limit);
                buy_limit.set_ordered_asset(s);
                buy_limit.set_order_size(demand, price);
                buy_limit.set_status(order::active);
                result_.emplace(order_num, buy_limit);

            }


        }


    for (auto &[i, j] : result_) {
        if (j.get_order_type()) {//this checks whether the order is a limit(1) or market(0) order
            if (j.get_order_size() > 0) {//this checks whether order is a bid(+) or ask(-)
                stack_bid(j.get_ordered_asset(), j);
            } else {
                stack_ask(j.get_ordered_asset(), j);
            }
        }

    }



//furthest range a market maker is allowed to charge away from the NBBO price

map<int,vector<order>> executed_orders_;
for (auto &[k, stock]:assets) {

    if(!this->bid.empty()) {
        if (this->bid.find(k) != this->bid.end()) {
            sort(this->bid.find(k)->second.begin(), this->bid.find(k)->second.end(), order_sort);

        }
    }

    if(!this->ask.empty()) {
        if (this->ask.find(k) != this->ask.end()) {
            sort(this->ask.find(k)->second.begin(), this->ask.find(k)->second.end(), order_sort);

        }
    }

    if(!this->market_orders.empty()){
        if(this->bid.find(k)!= this->bid.end()) {
            sort(this->market_orders.find(k)->second.begin(), this->market_orders.find(k)->second.end(), order_sort);
        }
    }


    {
        if(!this->market_orders.empty()) {
            std::cout << t << std::endl;
            std::cout << k << ".market " << "|";
            int size = this->market_orders.find(k)->second.size();
            for (int i = size - 1; i >= 0; i--) {
                std::cout << this->market_orders.find(k)->second.at(i).get_order_size() << "(" <<
                          this->market_orders.find(k)->second.at(i).get_id() << ")" << "(" <<
                          this->market_orders.find(k)->second.at(i).get_proposed_price() << ")" << "|";
            }
            std::cout << std::endl;
        }
        if(!this->bid.empty()) {
            std::cout << k << ".bid " << "|";
            int size = this->bid.find(k)->second.size();
            for (int i = size - 1; i >= 0; i--) {
                std::cout << this->bid.find(k)->second.at(i).get_order_size() << "(" <<
                          this->bid.find(k)->second.at(i).get_id() << ")" << "(" <<
                          this->bid.find(k)->second.at(i).get_proposed_price() << ")" << "|";
            }
            std::cout << std::endl;
        }

        if(!this->ask.empty()) {
            std::cout << k << ".ask " << "|";
            int size = this->ask.find(k)->second.size();
            for (int i = 0; i < size; i++) {
                std::cout << this->ask.find(k)->second.at(i).get_order_size() << "(" <<
                          this->ask.find(k)->second.at(i).get_id() << ")" << "(" <<
                          this->ask.find(k)->second.at(i).get_proposed_price() << ")" << "|";
            }
            std::cout << std::endl;
        }

}
//
//for(auto &[i,v]:this->trading_institutions){
//    std::cout << i<<" wealth "<<v->wealth << std::endl;
//}


//create price maps
//buy-limits
    map<double, vector<order>, greater<>> bid_orders = {};
    if(!this->bid.empty()) {
        if (this->bid.find(k) != this->bid.end()) {
            for (auto &x : this->bid.find(k)->second) {
                double best_bid_price = x.get_proposed_price();
                best_bid_price = (floor(best_bid_price * 100)) / 100;

                if (bid_orders.find(best_bid_price) != bid_orders.end()) {
                    bid_orders.find(best_bid_price)->second.push_back(x);
                } else {
                    bid_orders.emplace(best_bid_price, vector({x}));
                }


            }
        }
    }



    //limit sells
    map<double, vector<order>> ask_orders = {};
    if(!this->ask.empty()) {
        if (this->ask.find(k) != this->ask.end()) {
            for (auto &x : this->ask.find(k)->second) {
//        std::cout<<"price"<< <<std::endl;
                double best_ask_price = x.get_proposed_price();
                best_ask_price = (ceil(best_ask_price * 100)) / 100;

                if (ask_orders.find(best_ask_price) != ask_orders.end()) {
                    ask_orders.find(best_ask_price)->second.push_back(x);
                } else {
                    ask_orders.emplace(best_ask_price, vector({x}));
                }

            }
        }
    }



//match overlapping limit orders

// Price           bid | | ask
//  100          ......| |
//  101          ......| |
//  102            ****| |**
//  103             ***| |****
//  104                | |.....
//  105                | |......
//  106                | |.........

//bids are arranged in descending order while asks are in ascending order by price
    for(auto &[bid_price, bid_vec] : bid_orders) {
        for (auto &bid_order: bid_vec) {
            for (auto &[ask_price, ask_vec]: ask_orders) {
                for (auto &ask_order : ask_vec) {
                    if (bid_price >= ask_price) {

                        bool bid_is_unfilled = bid_order.get_status();
                        bool ask_is_unfilled = ask_order.get_status();

                        //non-stale orders
                        if (bid_is_unfilled && ask_is_unfilled) {
                            //trade

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
                                exec_bid.set_ordered_asset(stock.get_identifier());
                                exec_bid.set_order_size(BID, exec_bid_price);

                                bool buyer_exists = executed_orders_.find(bid_trader_id) != executed_orders_.end();
                                if(buyer_exists){
                                    executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
                                }else{
                                    executed_orders_.emplace(bid_trader_id, vector({exec_bid}));
                                }



                                //calculate executed ask
                                order exec_ask;
                                auto exec_ask_price = ask_price;
                                auto ask_trader_id = ask_order.get_id();

                                exec_ask.set_id(ask_trader_id);
                                exec_ask.set_status(order::active);
                                exec_ask.set_order_type(order::limit);
                                exec_ask.set_ordered_asset(stock.get_identifier());
                                exec_ask.set_order_size(-BID, exec_ask_price);

                                bool seller_exists = executed_orders_.find(ask_trader_id) != executed_orders_.end();
                                if(seller_exists){
                                    executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
                                }else{
                                    executed_orders_.emplace(ask_trader_id, vector({exec_ask}));
                                }

                                double bid_ask_spread = bid_price - ask_price;
                                this->cash_at_hand += bid_ask_spread * BID;

                                //update the ask orders remaining in the order vectors after execution of trade
                                ask_order.set_order_size(uncleared, ask_price);
                                //remove the fully filled bid order
                                bid_order.set_status(order::filled);

                            }
                            else
                            if (uncleared > 0.) {
                                //calculate executed bid
                                order exec_bid;
                                auto exec_bid_price = bid_price;
                                auto bid_trader_id = bid_order.get_id();

                                exec_bid.set_id(bid_trader_id);
                                exec_bid.set_status(order::active);
                                exec_bid.set_order_type(order::limit);
                                exec_bid.set_ordered_asset(stock.get_identifier());
                                exec_bid.set_order_size((BID - uncleared), exec_bid_price);

                                bool buyer_exists = executed_orders_.find(bid_trader_id) != executed_orders_.end();
                                if(buyer_exists){
                                    executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
                                }else{
                                    executed_orders_.emplace(bid_trader_id, vector({exec_bid}));
                                }




                                //calculate executed ask
                                order exec_ask;
                                auto exec_ask_price = ask_price;
                                auto ask_trader_id = ask_order.get_id();

                                exec_ask.set_id(ask_trader_id);
                                exec_ask.set_status(order::active);
                                exec_ask.set_order_type(order::limit);
                                exec_ask.set_ordered_asset(stock.get_identifier());
                                exec_ask.set_order_size(-(BID - uncleared), exec_ask_price);

                                bool seller_exists = executed_orders_.find(ask_trader_id) != executed_orders_.end();
                                if(seller_exists){
                                    executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
                                }else{
                                    executed_orders_.emplace(ask_trader_id, vector({exec_ask}));
                                }

                                //update cash earned from spread
                                double bid_ask_spread = bid_price - ask_price;
                                this->cash_at_hand += bid_ask_spread * BID;

                                bid_order.set_order_size(uncleared, bid_price);
                                ask_order.set_status(order::filled);
                            } else {
                                //calculate executed bid
                                order exec_bid;
                                auto exec_bid_price = bid_price;
                                auto bid_trader_id = bid_order.get_id();

                                exec_bid.set_id(bid_trader_id);
                                exec_bid.set_status(order::active);
                                exec_bid.set_order_type(order::limit);
                                exec_bid.set_ordered_asset(stock.get_identifier());
                                exec_bid.set_order_size(BID, exec_bid_price);

                                bool buyer_exists = executed_orders_.find(bid_trader_id) != executed_orders_.end();
                                if(buyer_exists){
                                    executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
                                }else{
                                    executed_orders_.emplace(bid_trader_id, vector({exec_bid}));
                                }



                                //calculate executed ask
                                order exec_ask;
                                auto exec_ask_price = ask_price;
                                auto ask_trader_id = ask_order.get_id();

                                exec_ask.set_id(ask_trader_id);
                                exec_ask.set_status(order::active);
                                exec_ask.set_order_type(order::limit);
                                exec_ask.set_ordered_asset(stock.get_identifier());
                                exec_ask.set_order_size(-BID, exec_ask_price);

                                bool seller_exists = executed_orders_.find(ask_trader_id) != executed_orders_.end();
                                if(seller_exists){
                                    executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
                                }else{
                                    executed_orders_.emplace(ask_trader_id, vector({exec_ask}));
                                }

                                //update cash earned from spread
                                double bid_ask_spread = bid_price - ask_price;
                                this->cash_at_hand += bid_ask_spread * BID;

                            }
                        ask_is_unfilled = bid_order.get_status();
                        bid_is_unfilled = ask_order.get_status();

                        }

                        if(ask_is_unfilled) {
                            //roll to next bid row but keep count of the ask we were now at to reduce computation burden
                            goto next_bid;
                        }
                    }else {
                        goto exit_loop;
                    }

                }
            }
            next_bid:;
        }
    }
    exit_loop:;


//match market orders
if(!this->market_orders.empty()){
    for (auto &market_order : this->market_orders.find(k)->second) {
        if (market_order.get_order_size() < 0.) {
            for (auto &[bid_price, bid_vec] : bid_orders) {
                for (auto &bid_order: bid_vec) {

                    bool bid_is_unfilled = bid_order.get_status();
                    bool mktask_is_unfilled = market_order.get_status();

                    if (bid_is_unfilled && mktask_is_unfilled) {
                        //trade

                        double mkt_cash = market_order.get_order_size() * market_order.get_proposed_price();
                        double BID = bid_order.get_order_size();
                        double ASK = mkt_cash / bid_price;
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
                            exec_bid.set_ordered_asset(stock.get_identifier());
                            exec_bid.set_order_size(BID, exec_bid_price);

                            bool buyer_exists = executed_orders_.find(bid_trader_id) != executed_orders_.end();
                            if (buyer_exists) {
                                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
                            } else {
                                executed_orders_.emplace(bid_trader_id, vector({exec_bid}));
                            }


                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = bid_price;
                            auto ask_trader_id = market_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(stock.get_identifier());
                            exec_ask.set_order_size(-BID, exec_ask_price);

                            bool seller_exists = executed_orders_.find(ask_trader_id) != executed_orders_.end();
                            if (seller_exists) {
                                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
                            } else {
                                executed_orders_.emplace(ask_trader_id, vector({exec_ask}));
                            }



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
                            exec_bid.set_ordered_asset(stock.get_identifier());
                            exec_bid.set_order_size((BID - uncleared), exec_bid_price);

                            bool buyer_exists = executed_orders_.find(bid_trader_id) != executed_orders_.end();
                            if (buyer_exists) {
                                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
                            } else {
                                executed_orders_.emplace(bid_trader_id, vector({exec_bid}));
                            }



                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = bid_price;
                            auto ask_trader_id = market_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(stock.get_identifier());
                            exec_ask.set_order_size(-(BID - uncleared), exec_ask_price);

                            bool seller_exists = executed_orders_.find(ask_trader_id) != executed_orders_.end();
                            if (seller_exists) {
                                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
                            } else {
                                executed_orders_.emplace(ask_trader_id, vector({exec_ask}));
                            }


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
                            exec_bid.set_ordered_asset(stock.get_identifier());
                            exec_bid.set_order_size(BID, exec_bid_price);

                            bool buyer_exists = executed_orders_.find(bid_trader_id) != executed_orders_.end();
                            if (buyer_exists) {
                                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
                            } else {
                                executed_orders_.emplace(bid_trader_id, vector({exec_bid}));
                            }



                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = bid_price;
                            auto ask_trader_id = market_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(stock.get_identifier());
                            exec_ask.set_order_size(-BID, exec_ask_price);

                            bool seller_exists = executed_orders_.find(ask_trader_id) != executed_orders_.end();
                            if (seller_exists) {
                                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
                            } else {
                                executed_orders_.emplace(ask_trader_id, vector({exec_ask}));
                            }


                            bid_order.set_status(order::filled);
                            market_order.set_status(order::filled);
                        }
                        bid_is_unfilled = bid_order.get_status();
                        mktask_is_unfilled = market_order.get_status();
                    }

                    //assess which of the orders is filled
                    if (bid_is_unfilled) {
                        //roll to next ask same row
                        goto next_mkt_order;
                    }
                }
            }
        } else if (market_order.get_order_size() > 0.) {
            for (auto &[ask_price, ask_vec] : ask_orders) {
                for (auto &ask_order: ask_vec) {

                    bool ask_is_unfilled = ask_order.get_status();
                    bool mktbid_is_unfilled = market_order.get_status();


                    if (ask_is_unfilled && mktbid_is_unfilled) {
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
                            exec_bid.set_ordered_asset(stock.get_identifier());
                            exec_bid.set_order_size(BID, exec_bid_price);

                            bool buyer_exists = executed_orders_.find(bid_trader_id) != executed_orders_.end();
                            if (buyer_exists) {
                                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
                            } else {
                                executed_orders_.emplace(bid_trader_id, vector({exec_bid}));
                            }


                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = ask_price;
                            auto ask_trader_id = ask_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(stock.get_identifier());
                            exec_ask.set_order_size(-BID, exec_ask_price);

                            bool seller_exists = executed_orders_.find(ask_trader_id) != executed_orders_.end();
                            if (seller_exists) {
                                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
                            } else {
                                executed_orders_.emplace(ask_trader_id, vector({exec_ask}));
                            }



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
                            exec_bid.set_ordered_asset(stock.get_identifier());
                            exec_bid.set_order_size((BID - uncleared), exec_bid_price);

                            bool buyer_exists = executed_orders_.find(bid_trader_id) != executed_orders_.end();
                            if (buyer_exists) {
                                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
                            } else {
                                executed_orders_.emplace(bid_trader_id, vector({exec_bid}));
                            }



                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = ask_price;
                            auto ask_trader_id = ask_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(stock.get_identifier());
                            exec_ask.set_order_size(-(BID - uncleared), exec_ask_price);

                            bool seller_exists = executed_orders_.find(ask_trader_id) != executed_orders_.end();
                            if (seller_exists) {
                                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
                            } else {
                                executed_orders_.emplace(ask_trader_id, vector({exec_ask}));
                            }

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
                            exec_bid.set_ordered_asset(stock.get_identifier());
                            exec_bid.set_order_size(BID, exec_bid_price);

                            bool buyer_exists = executed_orders_.find(bid_trader_id) != executed_orders_.end();
                            if (buyer_exists) {
                                executed_orders_.find(bid_trader_id)->second.push_back(exec_bid);
                            } else {
                                executed_orders_.emplace(bid_trader_id, vector({exec_bid}));
                            }




                            //calculate executed ask
                            order exec_ask;
                            auto exec_ask_price = ask_price;
                            auto ask_trader_id = ask_order.get_id();

                            exec_ask.set_id(ask_trader_id);
                            exec_ask.set_status(order::active);
                            exec_ask.set_order_type(order::market);
                            exec_ask.set_ordered_asset(stock.get_identifier());
                            exec_ask.set_order_size(-BID, exec_ask_price);

                            bool seller_exists = executed_orders_.find(ask_trader_id) != executed_orders_.end();
                            if (seller_exists) {
                                executed_orders_.find(ask_trader_id)->second.push_back(exec_ask);
                            } else {
                                executed_orders_.emplace(ask_trader_id, vector({exec_ask}));
                            }


                            market_order.set_status(order::filled);
                            ask_order.set_status(order::filled);
                        }
                        ask_is_unfilled = ask_order.get_status();
                        mktbid_is_unfilled = market_order.get_status();
                    }
                    //assess which of the orders is filled
                    if (ask_is_unfilled) {
                        //roll to next ask same row
                        goto next_mkt_order;
                    }
                }
            }
        }
        next_mkt_order:;
    }
}


    double prevailing_bid = 0;
    double prevailing_ask = 0;

    int x=0;
    for(auto &[price,vec]:bid_orders){
        for(auto &od : vec){
            bool bid_available = od.get_status();
            if(bid_available){
                prevailing_bid = price;
                x++;
                goto exit_point;
            }
        }
        exit_point:;
        if(x>0){
            break;
        }else{
            prevailing_bid= this->best_bid.find(k)->second;
        }
    }


    int y = 0;
    for(auto &[price,vec]:ask_orders){

        for(auto &ax_ : vec){
            bool ask_available = ax_.get_status();
            if(ask_available){
                prevailing_ask = price;
                y++;
                goto exit_point_;
            }
        }
        exit_point_:;
        if(y>0){
            break;
        }else{
                prevailing_ask = this->best_ask.find(k)->second;

        }
    }



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//inform traders about new prices

    for (auto &[i, v]: this->trading_institutions) {
        v->stocks_on_market.find(k)->second.set_price(
                this->clock.current_time()
                    , prevailing_bid
                        , prevailing_ask);



        this->assets.find(k)->second.set_price(
                this->clock.current_time()
                    , prevailing_bid
                    , prevailing_ask);

    }

    this->best_bid.find(k)->second = prevailing_bid;
    this->best_ask.find(k)->second = prevailing_ask;



}



    this->clock.tick();

if(executed_orders_.find(this->get_identifier()) != executed_orders_.end()){
    this->balance_bd(this->clock.current_time(),
                     executed_orders_.find(
                             this->get_identifier()
                     )->second);
}

    for (auto &[i, v]:this->trading_institutions) {
        bool trader_participated = executed_orders_.find(i) != executed_orders_.end();
        if (trader_participated) {
            auto exec_order = executed_orders_.find(i)->second;
            v->balance_bd(this->clock.current_time(), exec_order);
        }
        v->clock.tick();
    }



}


