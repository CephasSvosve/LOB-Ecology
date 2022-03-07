//
// Created by Cephas Svosve on 22/12/2021.
//

#include "funds.h"
#include <string>

funds::funds() {

}




map<int, order>
        funds::invest(){
    int a = this->clock.current_time();
    map<int, order> result;
    switch (this->fund_philosophy) {

        case 0:
            result= this->ag_value_demand();
            break;

        case 1:
            result= this->ag_noise_demand();
            break;
        case 2:
            if(a>window_size){
            result= this->ag_momentum_demand();}
            break;
        case 3:
            if (a>252){
            result= this->ag_growth_demand();}
            break;
        case 4:
            result= this->ag_index_demand();
            break;
        case 5:
            result= this->value_demand();
            break;
        case 6:
            result= this->noise_demand();
            break;
        case 7:
            if(a>window_size){
                result= this->momentum_demand();}
            break;
        case 8:
            if (a>252){
            result= this->growth_demand();}
            break;
        case 9:
            result= this->index_demand();
            break;
        default:
            break;

    }
    return result;

}


void funds::trade_strategy(trading_strategy tradingStrategy){
    this->fund_philosophy = tradingStrategy;
}


int funds::concat(int a, int b)
{

    // Convert both the integers to string
    string s1 = to_string(a);
    string s2 = to_string(b);

    // Concatenate both strings
    string s = s1 + s2;

    // Convert the concatenated string
    // to integer
    int c = stoi(s);

    // return the formed integer
    return c;
}

double min(double a, double b){

    double result = 0;
    bool a_is_minimum = a<b;
    if(a_is_minimum){
        result = a;
    }else{
        result = b;
    }

    return result;

}

double max(double a, double b){

    double result = 0;
    bool a_is_maximum = a>b;
    if(a_is_maximum){
        result = a;
    }else{
        result = b;
    }

    return result;

}
std::map<int,order>
funds::ag_value_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    double Beta = 9;//dummy threshold TODO: link this variable with actual earnings threshold
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();






//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid,ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto earnings = v.get_earnings(t);
        phi = -(log((quoted_price_)/earnings)-log(Beta));//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes){
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != i){current_position = i->second;}
        const auto &[time, bid,ask] = v.get_price();
        auto earnings = v.get_earnings(t);
        auto quoted_price_ = (bid+ask)/2.;


        phi = -(log((quoted_price_)/earnings)-log(Beta));
        auto allocation = phi/sum_of_signals; //stock wealth allocation




//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        int order_num;
        order_num = concat(t,k);

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;

//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;
            bool liquidity_required;

            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);



                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }
                //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_cash_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_cash_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings;this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings;this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
                    }

        //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else
        {
            auto disposable_cash = free_cash - ((1-portfolio_alloc)*investable_wealth);
            if((disposable_cash > 0) && (free_cash > 0)) {
                excess_demand_market = (disposable_cash * allocation * portfolio_alloc / (quoted_price_)) - current_position;

                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }

            //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_cash_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_cash_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings;this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings;this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }

        }
    }
    return result_;
}


std::map<int,order>
funds::ag_growth_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    //double Beta = 14;//dummy threshold TODO: link this variable with actual earnings threshold
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();




//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid,ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto E = v.get_earnings(t);
        auto E_o = v.get_earnings(t-252);
         auto E_growth = max(0,100*(E-E_o)/E_o);
        auto P_E      = quoted_price_/E;
        (E_growth) > 0 ?: E_growth = P_E;
        phi = -(log((P_E)/E_growth));//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes) {
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != i) { current_position = i->second; }
        const auto &[time, bid, ask] = v.get_price();
        auto quoted_price_ = (bid + ask) / 2.;
        auto E = v.get_earnings(t);
        auto E_o = v.get_earnings(t-252);
         auto E_growth = max(0,100*(E-E_o)/E_o);
        auto P_E      = quoted_price_/E;
        (E_growth) > 0 ?: E_growth = P_E;
        phi = -(log((P_E)/E_growth));//stock signal
        auto allocation = phi / sum_of_signals; //stock wealth allocation




//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        int order_num;
        order_num = concat(t,k);

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;

//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;
            bool liquidity_required;

            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);


                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }
            //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_income_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = E* E_growth;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = E* E_growth;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else
        {
            auto disposable_income = free_cash - ((1-portfolio_alloc)*investable_wealth);
            if((disposable_income > 0) && (free_cash > 0)) {
                excess_demand_market = (disposable_income * allocation * portfolio_alloc / (quoted_price_)) - current_position;

                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }

//            //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_income_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = E* E_growth;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = E* E_growth;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }

        }
    }
    return result_;
}


std::map<int,order>
funds::ag_momentum_demand() {

//Params definition
    map<int, order> result_;
    auto quotes = this->stocks_on_market;
    //double Beta = 9;//dummy threshold TODO: link this variable with actual earnings
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();




//compute a sum of signals, see McFadden choice function...........................(1)
    for (auto &[k, v] : quotes) {
        if (t > window_size) {
            auto quoted_price_ = (get<1>(v.get_price()) + get<2>(v.get_price()))/2.;
            int window_size_MA1 = 50;
            int window_size_MA2 = 200;

            vector<double> hist_prices_MA1 = v.get_price_range(window_size_MA1 - 1);
            vector<double> hist_prices_MA2 = v.get_price_range(window_size_MA2 - 1);
            double MA1 = std::accumulate(hist_prices_MA1.begin(), hist_prices_MA1.end(), 0.);
            double MA2 = std::accumulate(hist_prices_MA2.begin(), hist_prices_MA2.end(), 0.);

            double trend1 = (MA1 + (quoted_price_)) / window_size_MA1;
            double trend2 = (MA2 + (quoted_price_)) / window_size_MA2;
            if (MA2> 0.) { phi = ((trend1 / trend2) - 1); }
            else { phi = 0.; }//stock signal
            sum_of_signals = sum_of_signals + abs(phi);
        }
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for (auto &[k, v] : quotes) {
        double current_position = 0;
        auto i = stocks_at_hand.find(k);
        if (stocks_at_hand.end() != i) { current_position = i->second; }

        auto quoted_price_ = (get<1>(v.get_price()) + get<2>(v.get_price())) / 2.;
        int window_size_MA1 = 50;
        int window_size_MA2 = 200;

        vector<double> hist_prices_MA1 = v.get_price_range(window_size_MA1 - 1);
        vector<double> hist_prices_MA2 = v.get_price_range(window_size_MA2 - 1);
        double MA1 = std::accumulate(hist_prices_MA1.begin(), hist_prices_MA1.end(), 0.);
        double MA2 = std::accumulate(hist_prices_MA2.begin(), hist_prices_MA2.end(), 0.);

        double trend1 = (MA1 + (quoted_price_)) / window_size_MA1;
        double trend2 = (MA2 + (quoted_price_)) / window_size_MA2;
        if (MA2 > 0.) { phi = ((trend1 / trend2) - 1); }
        else { phi = 0.; }//stock signal


//stock wealth allocation
        auto allocation = phi / sum_of_signals;



//investment in the subject stock
        double excess_demand_market;
        order market_order;
        order limit_order;
        auto order_num = concat(t,k);

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;

//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;
            bool liquidity_required;

            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);


                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }
        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else
        {
            auto disposable_income = free_cash - ((1-portfolio_alloc)*investable_wealth);
            if((disposable_income > 0) && (free_cash > 0)) {
                excess_demand_market = (disposable_income * allocation * portfolio_alloc / (quoted_price_)) - current_position;

                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }

        }
    }
    return result_;
}


MatrixXd funds::generateWhiteNoise1(int seed,int rows, int columns){
    MatrixXd randoms(rows, columns);
    VectorXd a;


    //algorithm for generating random numbers that are seeded on changing time
    time_t now = seed;
    boost::random::mt19937 gen{static_cast<uint32_t>(now)};
    boost::normal_distribution<> nd(0.0, 1.0);
    boost::variate_generator<boost::mt19937 &,
    boost::normal_distribution<> > var_nor(gen, nd);



    //we generate the matrix of random numbers
    for (int rw = 0; rw < rows; rw++) {
        //we make sure the naturally occurring auto-correlation is sufficiently small by using a do-while loop
        //here we load each row with appropriate random numbers
        a = VectorXd(columns);
        for (int i = 0; i < columns; ++i) {
            a(i) = var_nor();}

        randoms.row(rw) = a;
    }
    return randoms;
}


double funds::lateralcorrcoef(VectorXd a){

    int n = a.size()-1;
    double x[n],y[n],xx[n],xy[n],yy[n];
    double sumx,sumy,sumxx,sumxy,sumyy;


    for (int i=0; i < n; i++)
    {
        x[i] = a[i+1];
        y[i] = a[i];
        xx[i] = pow(x[i],2);
        xy[i] = x[i] * y[i];
        yy[i] = pow(y[i],2);

        sumx = sumx + x[i];
        sumy = sumy + y[i];

        sumxx = sumxx + xx[i];
        sumxy = sumxy + xy[i];
        sumyy = sumyy + yy[i];
    }

    double r = (n*sumxy - sumx*sumy)/sqrt((n*sumxx-pow(sumx,2))*(n*sumyy-pow(sumy,2)));

    return r;
}


std::map<int,order>
funds::ag_noise_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    //dummy threshold TODO: link this variable with actual earnings
    double portfolio_alloc = 1;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();



//Ornstein Uhlenbeck params
    MatrixXd dX = generateWhiteNoise1(6,3,10000);

    //noise with a half life of 6 years, to match empirical evedence
    double mean_reversion_rate = 1 - pow(0.5, 1/(6*252.0));
  double sigma = 0.12;

    //Ornstein Uhlenbeck noise generation
    double  X_t ;
    int i =0;
    int c =0;
    map<int, vector<double>> noise;

    double Beta = 9;

//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid, ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto earnings = v.get_earnings(t);
        double var;
        double mean;
        var = (0.5*pow(sigma,2)) * (1/mean_reversion_rate)*
           (1-exp(-2*mean_reversion_rate*int(t)));
        mean =  0*(1-exp(-mean_reversion_rate * int(t)));

        X_t = mean + sqrt(var) * dX(i%3,t);
        phi = -(log((quoted_price_)/earnings)-log(Beta)) + X_t;//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
        i++;

    }


c=0;

//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes){
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto it = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != it){current_position = it->second;}
        const auto &[time, bid, ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto earnings = v.get_earnings(t);

        double var;
        double mean;
        var = (0.5*pow(sigma,2)) * (1/mean_reversion_rate)*
           (1-exp(-2*mean_reversion_rate*int(t)));
        mean =  0*(1-exp(-mean_reversion_rate * int(t)));
        X_t = mean + sqrt(var) * dX(c%3,t);
        c++;
        phi = -(log((quoted_price_)/earnings)-log(Beta)) + X_t;
        auto allocation = phi/sum_of_signals;//stock wealth allocation



//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;
        auto order_num = concat(t,k);

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;


//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;
            bool liquidity_required;

            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);


                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }
            //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_income_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings * exp(-X_t);this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings* exp(-X_t);this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else
        {
            auto disposable_income = free_cash - ((1-portfolio_alloc)*investable_wealth);
            if((disposable_income > 0) && (free_cash > 0)) {
                excess_demand_market = (disposable_income * allocation * portfolio_alloc / (quoted_price_)) - current_position;

                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }

            //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_income_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings* exp(-X_t);this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings* exp(-X_t);this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }

        }
    }
    return result_;
}


std::map<int, order> funds::ag_index_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    double portfolio_alloc = 1;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();






//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid,ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto shares_outstanding = v.get_shares_outstanding();
        phi = shares_outstanding * quoted_price_;//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes) {
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != i) { current_position = i->second; }
        const auto &[time, bid, ask] = v.get_price();
        auto shares_outstanding = v.get_shares_outstanding();
        auto quoted_price_ = (bid + ask) / 2.;


        phi = shares_outstanding * quoted_price_;
        auto allocation = phi / sum_of_signals; //stock wealth allocation




//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        // std::cout<<"earnings "<<earnings<<std::endl;
        int order_num;
        order_num = concat(t,k);
        //this->reb_period = 1;

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;

//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;
            bool liquidity_required;

            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,0);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);


                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }

        }
    }
    return result_;
}


std::map<int, order>
funds::value_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    double Beta = 9;//dummy threshold TODO: link this variable with actual earnings threshold
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();






//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid,ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto earnings = v.get_earnings(t);
        phi = -(log((quoted_price_)/earnings)-log(Beta));//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes){
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != i){current_position = i->second;}
        const auto &[time, bid,ask] = v.get_price();
        auto earnings = v.get_earnings(t);
        auto quoted_price_ = (bid+ask)/2.;


        phi = -(log((quoted_price_)/earnings)-log(Beta));
        auto allocation = phi/sum_of_signals; //stock wealth allocation




//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        int order_num;
        order_num = concat(t,k);

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;


//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;
            bool liquidity_required;

            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);


                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }
           // put limit orders to close positions when target is reached
            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {

                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
                order_num = concat(order_num, 1);
                excess_demand_limit = disposable_income_/quoted_price_;
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                auto target_price = Beta * earnings;this->target_price_ = target_price;
                limit_order.set_order_size(excess_demand_limit, target_price);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

            }else if(this->stocks_at_hand.find(k)->second>0){
                order_num = concat(order_num, 1);
                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                auto target_price = Beta * earnings;this->target_price_ = target_price;
                limit_order.set_order_size(excess_demand_limit, target_price);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);
            }
        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else
        {
            auto disposable_income = free_cash - ((1-portfolio_alloc)*investable_wealth);
            if((disposable_income > 0) && (free_cash > 0)) {
                excess_demand_market = (disposable_income * allocation * portfolio_alloc / (quoted_price_)) - current_position;

                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }

           // put limit orders to close positions when target is reached
            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {

                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
                order_num = concat(order_num, 1);
                excess_demand_limit = disposable_income_/quoted_price_;
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                auto target_price = Beta * earnings;this->target_price_ = target_price;
                limit_order.set_order_size(excess_demand_limit, target_price);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

            }else if(this->stocks_at_hand.find(k)->second>0){
                order_num = concat(order_num, 1);
                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                auto target_price = Beta * earnings;this->target_price_ = target_price;
                limit_order.set_order_size(excess_demand_limit, target_price);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);
            }

        }
    }
    return result_;
}


std::map<int,order>
funds::growth_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    //double Beta = 14;//dummy threshold TODO: link this variable with actual earnings threshold
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();




//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid,ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto E = v.get_earnings(t);
        auto E_o = v.get_earnings(t-252);
        auto E_growth = 100*(E-E_o)/E_o;
        auto P_E      = quoted_price_/E;
        (E_growth) > 0 ?: E_growth = P_E;
        phi = -(log((P_E)/E_growth));//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes) {
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != i) { current_position = i->second; }
        const auto &[time, bid, ask] = v.get_price();
        auto quoted_price_ = (bid + ask) / 2.;
        auto E = v.get_earnings(t);
        auto E_o = v.get_earnings(t-252);
         auto E_growth = max(0,100*(E-E_o)/E_o);
        auto P_E      = quoted_price_/E;
        (E_growth) > 0 ?: E_growth = P_E;
        phi = -(log(P_E/E_growth));//stock signal
        auto allocation = phi / sum_of_signals; //stock wealth allocation




//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        int order_num;
        order_num = concat(t,k);
        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;



//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;
            bool liquidity_required;

            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);


                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }
            //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_income_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = E* E_growth;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = E* E_growth;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else
        {
            auto disposable_income = free_cash - ((1-portfolio_alloc)*investable_wealth);
            if((disposable_income > 0) && (free_cash > 0)) {
                excess_demand_market = (disposable_income * allocation * portfolio_alloc / (quoted_price_)) - current_position;

                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }

            //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_income_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = E* E_growth;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = E* E_growth;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }

        }
    }
    return result_;
}


std::map<int,order>
funds::momentum_demand() {

//Params definition
    map<int, order> result_;
    auto quotes = this->stocks_on_market;
    //double Beta = 9;//dummy threshold TODO: link this variable with actual earnings
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();




//compute a sum of signals, see McFadden choice function...........................(1)
    for (auto &[k, v] : quotes) {
        if (t > window_size) {
            auto quoted_price_ = (get<1>(v.get_price()) + get<2>(v.get_price()))/2.;
            int window_size_MA1 = 50;
            int window_size_MA2 = 200;

            vector<double> hist_prices_MA1 = v.get_price_range(window_size_MA1 - 1);
            vector<double> hist_prices_MA2 = v.get_price_range(window_size_MA2 - 1);
            double MA1 = std::accumulate(hist_prices_MA1.begin(), hist_prices_MA1.end(), 0.);
            double MA2 = std::accumulate(hist_prices_MA2.begin(), hist_prices_MA2.end(), 0.);

            double trend1 = (MA1 + (quoted_price_)) / window_size_MA1;
            double trend2 = (MA2 + (quoted_price_)) / window_size_MA2;
            if (MA2> 0.) { phi = ((trend1 / trend2) - 1); }
            else { phi = 0.; }//stock signal
            sum_of_signals = sum_of_signals + abs(phi);
        }
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for (auto &[k, v] : quotes) {
        double current_position = 0;
        auto i = stocks_at_hand.find(k);
        if (stocks_at_hand.end() != i) { current_position = i->second; }

        auto quoted_price_ = (get<1>(v.get_price()) + get<2>(v.get_price())) / 2.;
        int window_size_MA1 = 50;
        int window_size_MA2 = 200;

        vector<double> hist_prices_MA1 = v.get_price_range(window_size_MA1 - 1);
        vector<double> hist_prices_MA2 = v.get_price_range(window_size_MA2 - 1);
        double MA1 = std::accumulate(hist_prices_MA1.begin(), hist_prices_MA1.end(), 0.);
        double MA2 = std::accumulate(hist_prices_MA2.begin(), hist_prices_MA2.end(), 0.);

        double trend1 = (MA1 + (quoted_price_)) / window_size_MA1;
        double trend2 = (MA2 + (quoted_price_)) / window_size_MA2;
        if (MA2 > 0.) { phi = ((trend1 / trend2) - 1); }
        else { phi = 0.; }//stock signal


//stock wealth allocation
        auto allocation = phi / sum_of_signals;



//investment in the subject stock
        double excess_demand_market;
        order market_order;
        order limit_order;
        auto order_num = concat(t,k);

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;

//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;
            bool liquidity_required;

            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);


                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }
        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else
        {
            auto disposable_income = free_cash - ((1-portfolio_alloc)*investable_wealth);
            if((disposable_income > 0) && (free_cash > 0)) {
                excess_demand_market = (disposable_income * allocation * portfolio_alloc / (quoted_price_)) - current_position;

                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }

        }
    }
    return result_;
}


std::map<int,order>
funds::noise_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    //dummy threshold TODO: link this variable with actual earnings
    double portfolio_alloc = 1;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();



//Ornstein Uhlenbeck params
    MatrixXd dX = generateWhiteNoise1(6,3,10000);

    //noise with a half life of 6 years, to match empirical evedence
    double mean_reversion_rate = 1 - pow(0.5, 1/(6*252.0));
   double sigma = 0.12;

    //Ornstein Uhlenbeck noise generation
    double  X_t ;
    int i =0;
    int c =0;
    map<int, vector<double>> noise;

    double Beta = 9;

//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid, ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto earnings = v.get_earnings(t);
        double var;
        double mean;
        var = (0.5*pow(sigma,2)) * (1/mean_reversion_rate)*
           (1-exp(-2*mean_reversion_rate*int(t)));
        mean =  0*(1-exp(-mean_reversion_rate * int(t)));

        X_t = mean + sqrt(var) * dX(i%3,t);


        phi = -(log((quoted_price_)/earnings)-log(Beta)) + X_t;//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
        i++;

    }


    c=0;
//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes){
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto it = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != it){current_position = it->second;}
        const auto &[time, bid, ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto earnings = v.get_earnings(t);

        double var;
        double mean;
        var = (0.5*pow(sigma,2)) * (1/mean_reversion_rate)*
           (1-exp(-2*mean_reversion_rate*int(t)));
        mean =  0*(1-exp(-mean_reversion_rate * int(t)));
        X_t = mean + sqrt(var) * dX(c%3,t);
        c++;
        phi = -(log((quoted_price_)/earnings)-log(Beta)) + X_t;
        auto allocation = phi/sum_of_signals;//stock wealth allocation



//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;
        auto order_num = concat(t,k);
        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;



//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;
            bool liquidity_required;

            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);


                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }
            //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_income_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings* exp(-X_t);this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings* exp(-X_t);this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else
        {
            auto disposable_income = free_cash - ((1-portfolio_alloc)*investable_wealth);
            if((disposable_income > 0) && (free_cash > 0)) {
                excess_demand_market = (disposable_income * allocation * portfolio_alloc / (quoted_price_)) - current_position;

                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }

            //put limit orders to close positions when target is reached
//            if((this->stocks_at_hand.find(k)->second<0) && (free_cash>0)) {
//
//                double disposable_income_ = min(abs(this->stocks_at_hand.find(k)->second*quoted_price_), free_cash);
//                order_num = concat(order_num, 1);
//                excess_demand_limit = disposable_income_/quoted_price_;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings* exp(-X_t);this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//            }else if(this->stocks_at_hand.find(k)->second>0){
//                order_num = concat(order_num, 1);
//                excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                auto target_price = Beta * earnings* exp(-X_t);this->target_price_ = target_price;
//                limit_order.set_order_size(excess_demand_limit, target_price);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }

        }
    }
    return result_;
}

std::map<int, order> funds::index_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    double portfolio_alloc = 1;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();






//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid,ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto shares_outstanding = v.get_shares_outstanding();
        phi = shares_outstanding * quoted_price_;//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes) {
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != i) { current_position = i->second; }
        const auto &[time, bid, ask] = v.get_price();
        auto shares_outstanding = v.get_shares_outstanding();
        auto quoted_price_ = (bid + ask) / 2.;


        phi = shares_outstanding * quoted_price_;
        auto allocation = phi / sum_of_signals; //stock wealth allocation




//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        int order_num;
        order_num = concat(t,k);
        //this->reb_period = 1;

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for(auto &[i,x]:this->stocks_at_hand){
            if(x<0){
                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else{
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = this->cash_at_hand - short_positions;
        investable_wealth = free_cash + long_positions;

//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if((this->reb_period*floor(t/this->reb_period) <= t) && (t < 5+this->reb_period*floor(t/this->reb_period) )){

            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
            excess_demand_market = target_position - current_position;


            //check if there are resources to make the desired trade
            if(0 <= target_position &&  target_position <= current_position){
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if((target_position< 0) && (current_position > 0)){
                //if target is a negative position but we currently have a positive number of stocks,
                // we first sell the stocks we have then short sell the remainder
                order_num = concat(order_num,1);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                market_order.set_ordered_asset(k);
                market_order.set_order_size(-current_position, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);


                double remaining_demand = min(abs(target_position*quoted_price_),free_cash);

                if(remaining_demand>0){
                    auto demand = (-remaining_demand/(quoted_price_));
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }else {
                double remaining_demand = min(abs(excess_demand_market*quoted_price_),free_cash);
                if(remaining_demand>0){
                    auto sign = excess_demand_market/abs(excess_demand_market);
                    auto demand = sign * remaining_demand/quoted_price_;
                    order_num = concat(order_num,1);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::limit);
                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }
            }

        }
    }
    return result_;
}
