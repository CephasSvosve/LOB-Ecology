//
// Created by Cephas Svosve on 22/12/2021.
//

#include "funds.h"
#include <string>
#include <fstream>

funds::funds() {

}




map<int, order>
funds::invest(){
    int a = this->clock->current_time();
    map<int, order> result;
    switch (this->fund_philosophy) {

        case 4:
            result= this->ag_value_demand();
            break;
        case 5:
            result= this->ag_noise_demand();
            break;
        case 6:
            if(a>window_size){
                result= this->ag_momentum_demand();}
            break;
        case 7:
            if (a>value_window){
                result= this->ag_growth_demand();}
            break;
        case 8:
            result= this->ag_index_demand();
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
funds::ag_value_demand() {

    //Params definition
    map<int, order> result_;
    auto quotes = this->stocks_on_market;
    // wealth allocated to stocks
    double portfolio_alloc = 0.6;
    //trading signal for the ith stock
    double phi;
    //sum of the trading signals
    double sum_of_signals;
    //time variable t
    int t = int(clock->current_time());
    //proportion of orders executed by limit entries
    double limit_order_proportion =0.5;


//computation of the sum of signals...........(1)
    for (auto &[k, v] : quotes) {
        const auto &[time, bid, ask] = v.get_price();
        //value parameter represents the intrinsic value of a stock k,
        double value = v.get_value(t);
        //a value investor will not purchase or sell stocks when price = intrinsic value
        double quoted_price_ = value;

                    //the decision to buy a stock is inferred from the ask price
                    //while the decision to sell a stock is inferred from the ask price
                    if(ask<=value){
                        quoted_price_ = ask;
                    }else
                        if(bid>=value){
                                quoted_price_ = bid;
                        }

        //phi = -(x + log(x) - 1)
        phi = -((quoted_price_/value) + log((quoted_price_)/value) -1);//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
    }



//computation of the ratio of each stock's signal relative to sum described in (1)
    for (auto &[k, v] : quotes) {
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != i) { current_position = i->second; }
        const auto &[time, bid, ask] = v.get_price();
        auto value = v.get_value(t);
        double quoted_price_ = value;
        if(ask<=value){
                quoted_price_ = ask;
            }else
            if(bid>=value){
                quoted_price_ = bid;
            }



        phi = -((quoted_price_/value) + log((quoted_price_)/value) -1);
        double  allocation = phi/sum_of_signals;




//investment in the subject stock
       //        if(allocation > 0){
//            quoted_price_ = ask;
//        }else{
//            quoted_price_ = bid;
//      }

        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        int order_num;
        order_num = concat(t, k);

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for (auto &[i, x]:this->stocks_at_hand) {
            if (x < 0) {
                short_positions += -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else {
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = max(0, (this->cash_at_hand - short_positions));



//rebalancing date
        if (rebalancing_count.find(k) == rebalancing_count.end()){
            rebalancing_count.emplace(k,0);
        }
        
        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {

            rebalancing_count.find(k)->second++;
            if(bond_at_hand < ((1-portfolio_alloc)*wealth)){
                double bond_purchase = min(((1-portfolio_alloc)*wealth)-bond_at_hand, free_cash);
                bond_at_hand += bond_purchase;
                this->cash_at_hand -= bond_purchase;
                free_cash -= bond_purchase;
            }else {
                double excess_bond = bond_at_hand - ((1-portfolio_alloc)*wealth);
                bond_at_hand -= excess_bond;
                this->cash_at_hand += excess_bond;
                free_cash = max(0, (this->cash_at_hand - short_positions));
            }


            auto T = free_cash * allocation*portfolio_alloc  / (quoted_price_);
            auto C = current_position;

            if(T!=0) {
                order_num = concat(1, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

                order_num = concat(2, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }
            if(abs(C+allocation)<abs(C)) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::market);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((1-limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

                order_num = concat(4, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);


            }else

            if(abs(C)>0) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size(-C, value);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);


            }

        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else {

            auto T = free_cash * allocation*portfolio_alloc   / (quoted_price_);
            auto C = current_position;

            if(T!=0) {
                order_num = concat(1, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

                order_num = concat(2, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }
            if(abs(C+allocation)<abs(C)) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::market);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((1-limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

                order_num = concat(4, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

            }else

            if(abs(C)>0) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size(-C, value);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

//                std::cout<<"normal_offload " <<"order "<<-C<<" price "<<value<<std::endl;
            }
        }
    }
    return result_;
}

MatrixXd funds::generateWhiteNoise1(int seed,int rows, int columns){
    MatrixXd randoms(rows, columns);
    VectorXd a;


    //algorithm for generating random numbers that are seeded on changing time
    time_t now = 6;
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
            a(i) = var_nor();
        }

        randoms.row(rw) = a;
    }
    return randoms;
}


std::map<int,order>
funds::ag_noise_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    //dummy threshold TODO: link this variable with actual value
    double portfolio_alloc = 1;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i

    int t = int(clock->current_time());
    double limit_order_proportion =0.5;



//Ornstein Uhlenbeck params


    //noise with a half life of 6 years, to match empirical evidence
    double mean_reversion_rate = 1 - pow(0.5, 1./(6*252.0));


    //Ornstein Uhlenbeck noise generation
    double  X_t ;

    int c =0;
    //map<int, vector<double>> noise;

    for(auto &[k, v] : quotes) {
                auto value = v.get_value(t);

                double var;
                double mean;
                double sigma;

                if (noise1.find(k) == noise1.end()) {
                    noise1.emplace(k, value / 2);
                }


                sigma = 0.12;
                var = (0.5 * pow(sigma, 2)) * (1 / mean_reversion_rate) * (1 - exp(-2 * mean_reversion_rate * 1));
                mean = noise1.find(k)->second * (exp(-mean_reversion_rate * 1)) + 0 * (1 - exp(-mean_reversion_rate * 1));
                X_t = mean + sqrt(var) * dX((k - 1), t);

std::cout<<"timeinnoise"<<t<<std::endl;
                if (day_count.find(k) == day_count.end()) {
                    day_count.emplace(k, 0);
                } else if (day_count.find(k)->second < t) {
                    noise1.find(k)->second = X_t;
                    day_count.find(k)->second++;
                }
    }
//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid, ask] = v.get_price();
        auto value = v.get_value(t);
        X_t =  noise1.find(k)->second;


        double quoted_price_ = value + X_t;

                //the decision to buy a stock is inferred from the ask price
                //while the decision to sell a stock is inferred from the ask price
                if(ask <= (value + X_t)){
                    quoted_price_ = ask;
                }else
                if(bid >= (value + X_t)){
                    quoted_price_ = bid;
                }






        phi = -(((quoted_price_)/(value + X_t))+log((quoted_price_)/(value + X_t))-1);//stock signal
        sum_of_signals = sum_of_signals + abs(phi);






    }




//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes){
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto it = stocks_at_hand.find(k);
        double current_position = 0;
        if (stocks_at_hand.end() != it){current_position = it->second;}
        const auto &[time, bid, ask] = v.get_price();
        auto value = v.get_value(t);
        X_t =  noise1.find(k)->second;

        double quoted_price_ = value + X_t;

            //the decision to buy a stock is inferred from the ask price
            //while the decision to sell a stock is inferred from the ask price
            if(ask<=(value + X_t)){
                quoted_price_ = ask;
            }else
            if(bid>=(value + X_t)){
                quoted_price_ = bid;
            }




        phi = -(((quoted_price_)/(value + X_t))+log((quoted_price_)/(value + X_t))-1);
        auto allocation = (phi/sum_of_signals);// = phi;//stock wealth allocation




        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        int order_num;
        order_num = concat(t, k);

        double free_cash = 0;
        double short_positions = 0;
        double long_positions = 0;
        double investable_wealth = 0;

        for (auto &[i, x]:this->stocks_at_hand) {
            if (x < 0) {
                short_positions += -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
            } else {
                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
            }
        }

        free_cash = max(0, (this->cash_at_hand - short_positions));


//rebalancing date
        if (rebalancing_count.find(k) == rebalancing_count.end()){
            rebalancing_count.emplace(k,0);
        }

        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {

            rebalancing_count.find(k)->second++;
            if(bond_at_hand < ((1-portfolio_alloc)*wealth)){
                double bond_purchase = min(((1-portfolio_alloc)*wealth)-bond_at_hand, free_cash);
                bond_at_hand += bond_purchase;
                this->cash_at_hand -= bond_purchase;
                free_cash -= bond_purchase;
            }else {
                double excess_bond = bond_at_hand - ((1-portfolio_alloc)*wealth);
                bond_at_hand -= excess_bond;
                this->cash_at_hand += excess_bond;
                free_cash = max(0, (this->cash_at_hand - short_positions));
            }


            auto T = free_cash * allocation*portfolio_alloc  / (quoted_price_);
            auto C = current_position;

            if(T!=0) {
                order_num = concat(1, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

                order_num = concat(2, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }
            if(abs(C+allocation)<abs(C)) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::market);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((1-limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

                order_num = concat(4, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);


            }else

            if(abs(C)>0) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size(-C, value + X_t);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);


            }

        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else {

            auto T = free_cash * allocation*portfolio_alloc   / (quoted_price_);
            auto C = current_position;

            if(T!=0) {
                order_num = concat(1, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

                order_num = concat(2, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }
            if(abs(C+allocation)<abs(C)) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::market);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((1-limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

                order_num = concat(4, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

            }else

            if(abs(C)>0) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size(-C, value + X_t);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);
            }
        }
    }
    return result_;
}

std::map<int,order>
funds::ag_growth_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    //double Beta = 14;//dummy threshold TODO: link this variable with actual value threshold
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    double limit_order_proportion =0.0;
    auto t = clock->current_time();





//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid,ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;


        vector<double> value = v.get_intrinsic_value();

        double MA1 = std::accumulate( value.begin() + (t-(value_window-252)),  value.begin() + (t), 0.)/(value_window-252);
        double MA2 = std::accumulate(  value.begin() + (t-(value_window)),  value.begin() + (t), 0.)
                /(value_window);


        double V = MA1;
        double V_o = MA2;
        phi = (V/V_o)+log(V/V_o) - 1;//stock signal
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


        vector<double> value = v.get_intrinsic_value();

        double MA1 = std::accumulate( value.begin() + (t-value_window),  value.begin() + (t), 0.)/value_window;
        double MA2 = std::accumulate(  value.begin() + (t-(value_window+252)),  value.begin() + (t), 0.)
                     /(value_window+252);


        double V = MA1;
        double V_o = MA2;
        phi = (V/V_o)+log(V/V_o) - 1; //stock signal
        auto allocation = (phi/sum_of_signals);// = phi / sum_of_signals; //stock wealth allocation

       //        if(allocation > 0){
//            quoted_price_ = ask;
//        }else{
//            quoted_price_ = bid;
//      }

        if (rebalancing_count.find(k) == rebalancing_count.end()){
            rebalancing_count.emplace(k,0);
        }

        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {

            auto is_empty = this->signal.find(k)== this->signal.end();
            if(is_empty){
                this->signal.emplace(k,allocation);
            }else{
                allocation = allocation;
            }

        }
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

        free_cash =max(0, (this->cash_at_hand - short_positions));


//rebalancing date
        if (rebalancing_count.find(k) == rebalancing_count.end()){
            rebalancing_count.emplace(k,0);
        }

        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {

            rebalancing_count.find(k)->second++;

            if(bond_at_hand < ((1-portfolio_alloc)*wealth)){
                double bond_purchase = min(((1-portfolio_alloc)*wealth)-bond_at_hand, free_cash);
                bond_at_hand += bond_purchase;
                this->cash_at_hand -= bond_purchase;
                free_cash -= bond_purchase;
            }else {
                double excess_bond = bond_at_hand - ((1-portfolio_alloc)*wealth);
                bond_at_hand -= excess_bond;
                this->cash_at_hand += excess_bond;
                free_cash = max(0, (this->cash_at_hand - short_positions));
            }


            auto T = free_cash * this->signal.find(k)->second*portfolio_alloc  / (quoted_price_);
            auto C = current_position;

            if(T!=0){
                order_num = concat(1, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

                order_num = concat(2, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);
            }


            if(abs(C+allocation)<=abs(C)) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::market);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((1-limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

                order_num = concat(4, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);


            }

        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else {

            auto T = free_cash * this->signal.find(k)->second  *portfolio_alloc  / (quoted_price_);
            auto C = current_position;

            if(T!=0){
                order_num = concat(1, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

                order_num = concat(2, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);
            }


        }
    }
    return result_;
}


std::map<int,order>
funds::ag_momentum_demand() {
//Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    //dummy threshold TODO: link this variable with actual value
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i

    int t = int(clock->current_time());
    double limit_order_proportion =0.5;




//compute a sum of signals, see McFadden choice function...........................(1)
    for (auto &[k, v] : quotes) {
        if (t > window_size) {
            auto quoted_price_ = v.get_midprice();
            int window_size_MA1 = 200;
            int window_size_MA2 = 500;

            vector<double> hist_prices_MA1 = v.get_price_range(window_size_MA1 - 1);
            vector<double> hist_prices_MA2 = v.get_price_range(window_size_MA2 - 1);
            double MA1 = std::accumulate(hist_prices_MA1.begin(), hist_prices_MA1.end(), 0.);
            double MA2 = std::accumulate(hist_prices_MA2.begin(), hist_prices_MA2.end(), 0.);

            double trend1 = (MA1 + (quoted_price_)) / window_size_MA1;
            double trend2 = (MA2 + (quoted_price_)) / window_size_MA2;
            if (MA2> 0.) { phi =  ((trend1 / trend2) + log(trend1 / trend2))-1; }
            else { phi = 0.; }//stock signal
            sum_of_signals = sum_of_signals + abs(phi);
        }
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for (auto &[k, v] : quotes) {
        double current_position = 0;
        auto i = stocks_at_hand.find(k);
        if (stocks_at_hand.end() != i) { current_position = i->second; }

        auto quoted_price_ = v.get_midprice();
        int window_size_MA1 = 200;
        int window_size_MA2 = 500;

        vector<double> hist_prices_MA1 = v.get_price_range(window_size_MA1 - 1);
        vector<double> hist_prices_MA2 = v.get_price_range(window_size_MA2 - 1);
        double MA1 = std::accumulate(hist_prices_MA1.begin(), hist_prices_MA1.end(), 0.);
        double MA2 = std::accumulate(hist_prices_MA2.begin(), hist_prices_MA2.end(), 0.);

        double trend1 = (MA1 + (quoted_price_)) / window_size_MA1;
        double trend2 = (MA2 + (quoted_price_)) / window_size_MA2;
        if (MA2 > 0.) { phi = ((trend1 / trend2) + log(trend1 / trend2))-1; }
        else { phi = 0.; }//stock signal


//stock wealth allocation
        double allocation = (phi/sum_of_signals);// = phi / sum_of_signals;
        const auto &[time, bid,ask] = v.get_price();


        if (rebalancing_count.find(k) == rebalancing_count.end()){
            rebalancing_count.emplace(k,0);
        }
        
        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {
            auto is_empty = this->signal.find(k)== this->signal.end();
            if(is_empty){
                this->signal.emplace(k,allocation);
            }else{
                allocation = allocation;
            }

        }


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

        free_cash =max(0, (this->cash_at_hand - short_positions));



//rebalancing date
        if (rebalancing_count.find(k) == rebalancing_count.end()){
            rebalancing_count.emplace(k,0);
        }
        
        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {

            rebalancing_count.find(k)->second++;

            if(bond_at_hand < ((1-portfolio_alloc)*wealth)){
                double bond_purchase = min(((1-portfolio_alloc)*wealth)-bond_at_hand, free_cash);
                bond_at_hand += bond_purchase;
                this->cash_at_hand -= bond_purchase;
                free_cash -= bond_purchase;
            }else {
                double excess_bond = bond_at_hand - ((1-portfolio_alloc)*wealth);
                bond_at_hand -= excess_bond;
                this->cash_at_hand += excess_bond;
                free_cash = max(0, (this->cash_at_hand - short_positions));
            }


            auto T = free_cash * this->signal.find(k)->second * portfolio_alloc  / (quoted_price_);
            auto C = current_position;

            if(T!=0){
            order_num = concat(1, order_num);
            market_order.set_id(this->get_identifier());
            market_order.set_order_type(order::market);
            market_order.set_ordered_asset(k);
            market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
            market_order.set_status(order::active);
            result_.emplace(order_num, market_order);

            order_num = concat(2, order_num);
            market_order.set_id(this->get_identifier());
            market_order.set_order_type(order::limit);
            market_order.set_ordered_asset(k);
            market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
            market_order.set_status(order::active);
            result_.emplace(order_num, market_order);
            }


            if(abs(C+allocation)<=abs(C)) {
                order_num = concat(3, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::market);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((1-limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);

                order_num = concat(4, order_num);
                limit_order.set_id(this->get_identifier());
                limit_order.set_order_type(order::limit);
                limit_order.set_ordered_asset(k);
                limit_order.set_order_size((limit_order_proportion)*-C, quoted_price_);
                limit_order.set_status(order::active);
                result_.emplace(order_num, limit_order);


            }
        }

            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
        else {
//
            auto T = free_cash * this->signal.find(k)->second *portfolio_alloc  / (quoted_price_);

            if(T!=0){
                order_num = concat(1, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

                order_num = concat(2, order_num);
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::limit);
                market_order.set_ordered_asset(k);
                market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);
            }

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
    auto t = clock->current_time();






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
        auto allocation = (phi/sum_of_signals);// = phi / sum_of_signals; //stock wealth allocation




//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        // std::cout<<"value "<<value<<std::endl;
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

        free_cash =max(0, (this->cash_at_hand - short_positions));
        investable_wealth = free_cash + long_positions;

        auto disposable_income = max(0,(free_cash - ((1-portfolio_alloc)*investable_wealth)));

//rebalancing date
        if (rebalancing_count.find(k) == rebalancing_count.end()){
            rebalancing_count.emplace(k,0);
        }
        
        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {

            rebalancing_count.find(k)->second++;

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


                double remaining_demand = min(abs(target_position*quoted_price_),
                                              abs(disposable_income*allocation*portfolio_alloc));

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


                double remaining_demand = min(abs(target_position*quoted_price_),
                                              abs(disposable_income*allocation*portfolio_alloc));

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

        }else{
            excess_demand_market = disposable_income*allocation*portfolio_alloc/quoted_price_;

            order_num = concat(order_num,0);
            market_order.set_id(this->get_identifier());
            market_order.set_order_type(order::market);
            std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
            market_order.set_ordered_asset(k);
            market_order.set_order_size(excess_demand_market, quoted_price_);
            market_order.set_status(order::active);
            result_.emplace(order_num, market_order);

            auto target_position = investable_wealth * allocation * portfolio_alloc / (quoted_price_);

            if ((target_position > current_position) && (current_position <= 0)) {
                if (target_position < 0) {
                    auto demand = min(abs(target_position - current_position), abs(current_position));
                    order_num = concat(order_num, 0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(abs(demand), quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                } else {
                    auto demand = min(abs(target_position - current_position), abs(current_position));
                    order_num = concat(order_num, 0);
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(abs(current_position), quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }


        }
    }
    return result_;
}

//
//std::map<int, order>
//funds::value_demand(){
//    //Params definition
//    map<int, order> result_;
//    auto quotes = this->stocks_on_market;
//    double Beta = 1;//dummy threshold TODO: link this variable with actual value threshold
//    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
//    double sum_of_signals;
//    double phi; //trading signal for stock i
//    auto t = clock->current_time();
//
//
//
//
//
//
////compute a sum of exponents of signals, see McFadden choice function...........................(1)
//    for (auto &[k, v] : quotes) {
//        const auto &[time, bid, ask] = v.get_price();
//        double value = v.get_value(t);
//        double quoted_price_ = 0;//v.get_midprice();
//        if(ask<=value){
//            auto fg = v;
//            quoted_price_ = ask;
//        }else
//        if(bid>=value){
//            quoted_price_ = bid;
//        }else{
//            quoted_price_ = v.get_midprice();
//        }
//
//
//        phi = -((quoted_price_)/value + log((quoted_price_)/value) -1);//stock signal
//        sum_of_signals = sum_of_signals + abs(phi);
//    }
//
//
//
////compute the ratio of each stock's signal relative to sum described in (1)
//    for (auto &[k, v] : quotes) {
//        //the following line searches whether a trader already holds inventory of the stock in subject
//        // and stores the value in j
//        auto i = stocks_at_hand.find(k);
//        double current_position = 0;
//        if (stocks_at_hand.end() != i) { current_position = i->second; }
//        const auto &[time, bid, ask] = v.get_price();
//        auto value = v.get_value(t);
//        double quoted_price_ = 0;//v.get_midprice();
//        if(ask<=value){
//            quoted_price_ = ask;
//        }else
//        if(bid>=value){
//            quoted_price_ = bid;
//        }else{
//            quoted_price_ = v.get_midprice();
//        }
//
//
//        phi = -((quoted_price_)/value + log((quoted_price_)/value) -1);
//        double  allocation = phi/sum_of_signals;;
//
//
//
//
////investment in the subject stock
//       //        if(allocation > 0){
////            quoted_price_ = ask;
////        }else{
////            quoted_price_ = bid;
////      }
//
//        double excess_demand_market;
//        double excess_demand_limit;
//        order market_order;
//        order limit_order;
//
//
//        int order_num;
//        order_num = concat(t, k);
//
//        double free_cash = 0;
//        double short_positions = 0;
//        double long_positions = 0;
//        double investable_wealth = 0;
//
//        for (auto &[i, x]:this->stocks_at_hand) {
//            if (x < 0) {
//                short_positions += -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            } else {
//                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            }
//        }
//
//        free_cash = max(0, (this->cash_at_hand - short_positions));
//
//
//
////rebalancing date
//        if (rebalancing_count.find(k) == rebalancing_count.end()){
//            rebalancing_count.emplace(k,0);
//        }
//
//        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {
//
//            rebalancing_count.find(k)->second++;
//
//            if(bond_at_hand < ((1-portfolio_alloc)*wealth)){
//                double bond_purchase = min(((1-portfolio_alloc)*wealth)-bond_at_hand, free_cash);
//                bond_at_hand += bond_purchase;
//                this->cash_at_hand -= bond_purchase;
//                free_cash -= bond_purchase;
//            }else {
//                double excess_bond = bond_at_hand - ((1-portfolio_alloc)*wealth);
//                bond_at_hand -= excess_bond;
//                this->cash_at_hand += excess_bond;
//                free_cash = max(0, (this->cash_at_hand - short_positions));
//            }
//
//
//            auto T = free_cash * allocation*portfolio_alloc  / (quoted_price_);
//            auto C = current_position;
//
//            if(T!=0) {
//                order_num = concat(1, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size(T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//            }
//
//            if(abs(C)>0) {
//                order_num = concat(2, order_num);
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                limit_order.set_order_size(-C, value);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
//
//        }
//
//            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
//        else {
//
//            auto T = free_cash * allocation*portfolio_alloc   / (quoted_price_);
//            auto C = current_position;
//
//            if(T!=0) {
//                order_num = concat(1, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size(T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//            }
//
//            if(abs(C)>0) {
//                order_num = concat(2, order_num);
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
//                limit_order.set_ordered_asset(k);
//                limit_order.set_order_size(-C, value);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
//        }
//    }
//    return result_;
//}
//std::map<int,order>
//funds::growth_demand(){
//    //Params definition
//    map<int,order> result_;
//    auto quotes = this->stocks_on_market;
//    //double Beta = 14;//dummy threshold TODO: link this variable with actual value threshold
//    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
//    double sum_of_signals;
//    double phi; //trading signal for stock i
//    double limit_order_proportion =0.;
//    auto t = clock->current_time();
//
//
//
//
//
////compute a sum of exponents of signals, see McFadden choice function...........................(1)
//    for(auto &[k, v] : quotes){
//        const auto &[time, bid,ask] = v.get_price();
//        auto quoted_price_ = (bid + ask)/2.;
//        double V = v.get_value(t);
//        double V_o = v.get_value(t-252);
//        phi = (V/V_o)+log(V/V_o) - 1;;//stock signal
//        sum_of_signals = sum_of_signals + abs(phi);
//    }
//
//
//
////compute the ratio of each stock's signal relative to sum described in (1)
//    for(auto &[k, v] : quotes) {
//        //the following line searches whether a trader already holds inventory of the stock in subject
//        // and stores the value in j
//        auto i = stocks_at_hand.find(k);
//        double current_position = 0;
//        if (stocks_at_hand.end() != i) { current_position = i->second; }
//        const auto &[time, bid, ask] = v.get_price();
//        auto quoted_price_ = (bid + ask) / 2.;
//        double V = v.get_value(t);
//        double V_o = v.get_value(t-252);
//        phi = (V/V_o)+log(V/V_o) - 1; //stock signal
//        auto allocation = (phi/sum_of_signals);// = phi / sum_of_signals; //stock wealth allocation
//
//       //        if(allocation > 0){
////            quoted_price_ = ask;
////        }else{
////            quoted_price_ = bid;
////      }
////investment in the subject stock
//        double excess_demand_market;
//        double excess_demand_limit;
//        order market_order;
//        order limit_order;
//
//
//        int order_num;
//        order_num = concat(t,k);
//        double free_cash = 0;
//        double short_positions = 0;
//        double long_positions = 0;
//        double investable_wealth = 0;
//
//        for(auto &[i,x]:this->stocks_at_hand){
//            if(x<0){
//                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            } else{
//                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            }
//        }
//
//        free_cash =max(0, (this->cash_at_hand - short_positions));
//
//
////rebalancing date
//        if (rebalancing_count.find(k) == rebalancing_count.end()){
//            rebalancing_count.emplace(k,0);
//        }
//
//        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {
//
//            rebalancing_count.find(k)->second++;
//
//            if(bond_at_hand < ((1-portfolio_alloc)*wealth)){
//                double bond_purchase = min(((1-portfolio_alloc)*wealth)-bond_at_hand, free_cash);
//                bond_at_hand += bond_purchase;
//                this->cash_at_hand -= bond_purchase;
//                free_cash -= bond_purchase;
//            }else {
//                double excess_bond = bond_at_hand - ((1-portfolio_alloc)*wealth);
//                bond_at_hand -= excess_bond;
//                this->cash_at_hand += excess_bond;
//                free_cash = max(0, (this->cash_at_hand - short_positions));
//            }
//
//
//            auto T = free_cash * allocation*portfolio_alloc  / (quoted_price_);
//            auto C = current_position;
//
//            if(T!=0){
//                order_num = concat(1, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::market);
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//
//                order_num = concat(2, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//            }
//
//
//            if(abs(C+allocation)<=abs(C)) {
//                order_num = concat(3, order_num);
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::market);
//                limit_order.set_ordered_asset(k);
//                limit_order.set_order_size((1-limit_order_proportion)*-C, quoted_price_);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//                order_num = concat(4, order_num);
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                limit_order.set_order_size((limit_order_proportion)*-C, quoted_price_);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//
//            }
//
//        }
//
//            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
//        else {
//
//            auto T = free_cash * allocation*portfolio_alloc  / (quoted_price_);
//            auto C = current_position;
//
//            if(T!=0){
//                order_num = concat(1, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::market);
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size((1-limit_order_proportion)*T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//
//                order_num = concat(2, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size((limit_order_proportion)*T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//            }
//
//
//            if(abs(C+allocation)<=abs(C)) {
//                order_num = concat(3, order_num);
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::market);
//                limit_order.set_ordered_asset(k);
//                limit_order.set_order_size((1-limit_order_proportion)*-C, quoted_price_);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//                order_num = concat(4, order_num);
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                limit_order.set_order_size((limit_order_proportion)*-C, quoted_price_);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//
//
//            }
//
//        }
//    }
//    return result_;
//}
//
//
//std::map<int,order>
//funds::momentum_demand() {
//
////Params definition
//    map<int, order> result_;
//    auto quotes = this->stocks_on_market;
//    //double Beta = 1;//dummy threshold TODO: link this variable with actual value
//    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
//    double sum_of_signals;
//    double phi; //trading signal for stock i
//    auto t = clock->current_time();
//    double limit_prop = 0.2;
//
//
//
//
////compute a sum of signals, see McFadden choice function...........................(1)
//    for (auto &[k, v] : quotes) {
//        if (t > window_size) {
//            auto quoted_price_ = (get<1>(v.get_price()) + get<2>(v.get_price()))/2.;
//            int window_size_MA1 = 250;
//            int window_size_MA2 = 500;
//
//            vector<double> hist_prices_MA1 = v.get_price_range(window_size_MA1 - 1);
//            vector<double> hist_prices_MA2 = v.get_price_range(window_size_MA2 - 1);
//            double MA1 = std::accumulate(hist_prices_MA1.begin(), hist_prices_MA1.end(), 0.);
//            double MA2 = std::accumulate(hist_prices_MA2.begin(), hist_prices_MA2.end(), 0.);
//
//            double trend1 = (MA1 + (quoted_price_)) / window_size_MA1;
//            double trend2 = (MA2 + (quoted_price_)) / window_size_MA2;
//            if (MA2> 0.) { phi =  ((trend1 / trend2) + log(trend1 / trend2))-1; }
//            else { phi = 0.; }//stock signal
//            sum_of_signals = sum_of_signals + abs(phi);
//        }
//    }
//
//
//
////compute the ratio of each stock's signal relative to sum described in (1)
//    for (auto &[k, v] : quotes) {
//        double current_position = 0;
//        auto i = stocks_at_hand.find(k);
//        if (stocks_at_hand.end() != i) { current_position = i->second; }
//
//        auto quoted_price_ = (get<1>(v.get_price()) + get<2>(v.get_price())) / 2.;
//        int window_size_MA1 = 50;
//        int window_size_MA2 = 200;
//
//        vector<double> hist_prices_MA1 = v.get_price_range(window_size_MA1 - 1);
//        vector<double> hist_prices_MA2 = v.get_price_range(window_size_MA2 - 1);
//        double MA1 = std::accumulate(hist_prices_MA1.begin(), hist_prices_MA1.end(), 0.);
//        double MA2 = std::accumulate(hist_prices_MA2.begin(), hist_prices_MA2.end(), 0.);
//
//        double trend1 = (MA1 + (quoted_price_)) / window_size_MA1;
//        double trend2 = (MA2 + (quoted_price_)) / window_size_MA2;
//        if (MA2 > 0.) { phi = ((trend1 / trend2) + log(trend1 / trend2))-1; }
//        else { phi = 0.; }//stock signal
//
//
////stock wealth allocation
//        double allocation = (phi/sum_of_signals);// = phi / sum_of_signals;
//        const auto &[time, bid,ask] = v.get_price();
//       //        if(allocation > 0){
////            quoted_price_ = ask;
////        }else{
////            quoted_price_ = bid;
////      }
//
//        if (rebalancing_count.find(k) == rebalancing_count.end()){
//            rebalancing_count.emplace(k,0);
//        }
//
//        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {
//            auto is_empty = this->signal.find(k)== this->signal.end();
//            if(is_empty){
//                this->signal.emplace(k,allocation);
//            }else{
//                allocation = allocation;
//            }
//
//        }
//
////        if(sum_of_signals>1){
////            allocation = phi/sum_of_signals;
////        }else{
////            allocation = (phi);
////        }
//
//
////investment in the subject stock
//        double excess_demand_market;
//        order market_order;
//        order limit_order;
//        auto order_num = concat(t,k);
//
//        double free_cash = 0;
//        double short_positions = 0;
//        double long_positions = 0;
//        double investable_wealth = 0;
//
//        for(auto &[i,x]:this->stocks_at_hand){
//            if(x<0){
//                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            } else{
//                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            }
//        }
//
//        free_cash =max(0, (this->cash_at_hand - short_positions));
//
//
//
////rebalancing date
//        if (rebalancing_count.find(k) == rebalancing_count.end()){
//            rebalancing_count.emplace(k,0);
//        }
//
//        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {
//
//            rebalancing_count.find(k)->second++;
//
//            if(bond_at_hand < ((1-portfolio_alloc)*wealth)){
//                double bond_purchase = min(((1-portfolio_alloc)*wealth)-bond_at_hand, free_cash);
//                bond_at_hand += bond_purchase;
//                this->cash_at_hand -= bond_purchase;
//                free_cash -= bond_purchase;
//            }else {
//                double excess_bond = bond_at_hand - ((1-portfolio_alloc)*wealth);
//                bond_at_hand -= excess_bond;
//                this->cash_at_hand += excess_bond;
//                free_cash = max(0, (this->cash_at_hand - short_positions));
//            }
//
//
//            auto T = free_cash * this->signal.find(k)->second * portfolio_alloc  / (quoted_price_);
//            auto C = current_position;
//
//            if(T!=0){
//                order_num = concat(1, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size(T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);}
//
//
//            if(abs(C+allocation)<=abs(C)) {
//                order_num = concat(2, order_num);
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
//                limit_order.set_ordered_asset(k);
//                limit_order.set_order_size(-C, quoted_price_);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
//        }
//
//            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
//        else {
////            if (this->stocks_at_hand.find(k)!= this->stocks_at_hand.end()) {
////                allocation = (this->stocks_at_hand.find(k)->second*this->stocks_on_market.find(k)->second.get_midprice())/(short_positions+long_positions);
////            } else{
////                allocation = 0.;
////            }
//
//            auto T = free_cash * allocation *portfolio_alloc  / (quoted_price_);
//
//            if(T!=0){
//                order_num = concat(1, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size(T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);}
////////
////
//        }
//    }
//    return result_;
//}
//
//
//std::map<int,order>
//funds::noise_demand(){
//    //Params definition
//    map<int,order> result_;
//    auto quotes = this->stocks_on_market;
//    //dummy threshold TODO: link this variable with actual value
//    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
//    double sum_of_signals;
//    double phi; //trading signal for stock i
//    auto t = clock->current_time();
//
//
//
////Ornstein Uhlenbeck params
//    MatrixXd dX = generateWhiteNoise1(6,3,10000);
//
//    //noise with a half life of 6 years, to match empirical evidence
//    double mean_reversion_rate = 1 - pow(0.5, 1./(6*252.0));
//
//
//    //Ornstein Uhlenbeck noise generation
//    double  X_t ;
//    int i =0;
//    int c =0;
//    map<int, vector<double>> noise;
//
//
//
////compute a sum of exponents of signals, see McFadden choice function...........................(1)
//    for(auto &[k, v] : quotes){
//        const auto &[time, bid, ask] = v.get_price();
//        auto value = v.get_value(t);
//
//        double quoted_price_ =0.;
//        if(ask<=value){
//            quoted_price_ = ask;
//        }else
//        if(bid>=value){
//            quoted_price_ = bid;
//        }else{
//            quoted_price_ = v.get_midprice();
//        }
//
//
//        double var;
//        double mean;
//
//        if (noise1.find(k) == noise1.end()){
//            noise1.emplace(k,value/4);
//        }
//double sigma = 0.2;//sqrt(log((0.04/(pow(value,2)))+1));
//        double mu = v.get_value(t);// log(value)-0.5*pow(sigma,2);
//
//        var = (0.5*pow(sigma,2)) * (1/mean_reversion_rate)*(1-exp(-2*mean_reversion_rate*1));
//        mean = noise1.find(k)->second*(exp(-mean_reversion_rate * 1))+ 0*(1-exp(-mean_reversion_rate * 1));
//
//        X_t =  mean + sqrt(var) * dX(i%3,t);
//        phi = -(((quoted_price_)/(value + X_t))+log((quoted_price_)/(value + X_t))-1);//stock signal
//        sum_of_signals = sum_of_signals + abs(phi);
//
//        noise1.find(k)->second = X_t;
//        i++;
//
//    }
//
//
//    c=0;
//
////compute the ratio of each stock's signal relative to sum described in (1)
//    for(auto &[k, v] : quotes){
//        //the following line searches whether a trader already holds inventory of the stock in subject
//        // and stores the value in j
//        auto it = stocks_at_hand.find(k);
//        double current_position = 0;
//        if (stocks_at_hand.end() != it){current_position = it->second;}
//        const auto &[time, bid, ask] = v.get_price();
//        auto value = v.get_value(t);
//
//        double quoted_price_ =0.;
//        if(ask<=value){
//            quoted_price_ = ask;
//        }else
//        if(bid>=value){
//            quoted_price_ = bid;
//        }else{
//            quoted_price_ = v.get_midprice();
//        }
//
//        double var;
//        double mean;
//        if (noise1.find(k) == noise1.end()){
//            noise1.emplace(k,value/4);
//        }
//double sigma = 0.2;//sqrt(log((0.04/(pow(value,2)))+1));
//        double mu = v.get_value(t);// log(value)-0.5*pow(sigma,2);
//
//        var = (0.5*pow(sigma,2)) * (1/mean_reversion_rate)*(1-exp(-2*mean_reversion_rate*1));
//        mean = noise1.find(k)->second*(exp(-mean_reversion_rate * 1))+ 0*(1-exp(-mean_reversion_rate * 1));
//
//        X_t =  mean + sqrt(var) * dX(i%3,t);
//        c++;
//        phi = -(((quoted_price_)/(value + X_t))+log((quoted_price_)/(value + X_t))-1);
//        auto allocation = (phi/sum_of_signals);// = phi;//stock wealth allocation
//
//        noise1.find(k)->second = X_t;
//
////investment in the subject stock
//       //        if(allocation > 0){
////            quoted_price_ = ask;
////        }else{
////            quoted_price_ = bid;
////      }
//
//        double excess_demand_market;
//        double excess_demand_limit;
//        order market_order;
//        order limit_order;
//
//
//        int order_num;
//        order_num = concat(t, k);
//
//        double free_cash = 0;
//        double short_positions = 0;
//        double long_positions = 0;
//        double investable_wealth = 0;
//
//        for (auto &[i, x]:this->stocks_at_hand) {
//            if (x < 0) {
//                short_positions += -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            } else {
//                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            }
//        }
//
//        free_cash = max(0, (this->cash_at_hand - short_positions));
//
//
////rebalancing date
//        if (rebalancing_count.find(k) == rebalancing_count.end()){
//            rebalancing_count.emplace(k,0);
//        }
//
//        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {
//
//            rebalancing_count.find(k)->second++;
//
//            if(bond_at_hand < ((1-portfolio_alloc)*wealth)){
//                double bond_purchase = min(((1-portfolio_alloc)*wealth)-bond_at_hand, free_cash);
//                bond_at_hand += bond_purchase;
//                this->cash_at_hand -= bond_purchase;
//                free_cash -= bond_purchase;
//            }else {
//                double excess_bond = bond_at_hand - ((1-portfolio_alloc)*wealth);
//                bond_at_hand -= excess_bond;
//                this->cash_at_hand += excess_bond;
//                free_cash = max(0, (this->cash_at_hand - short_positions));
//            }
//
//
//            auto T = free_cash * allocation*portfolio_alloc  / (quoted_price_);
//            auto C = current_position;
//
//            if(T!=0) {
//                order_num = concat(1, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size(T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//            }
//
//            if(abs(C)>0) {
//                order_num = concat(2, order_num);
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                limit_order.set_ordered_asset(k);
//                limit_order.set_order_size(-C, value + X_t);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
//
//        }
//
//            //when time is not a rebalancing date, adjust positions using excess cash (e.g. from dividends)
//        else {
//
//            auto T = free_cash * allocation * portfolio_alloc  / (quoted_price_);
//            auto C = current_position;
//
//            if(T!=0) {
//                order_num = concat(1, order_num);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size(T, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//            }
//
//            if(abs(C)>0) {
//                order_num = concat(2, order_num);
//                limit_order.set_id(this->get_identifier());
//                limit_order.set_order_type(order::limit);
//                std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
//                limit_order.set_ordered_asset(k);
//                limit_order.set_order_size(-C, value + X_t);
//                limit_order.set_status(order::active);
//                result_.emplace(order_num, limit_order);
//            }
//        }
//    }
//    return result_;
//}
//
//std::map<int, order> funds::index_demand(){
//    //Params definition
//    map<int,order> result_;
//    auto quotes = this->stocks_on_market;
//    double portfolio_alloc = 1;// wealth allocated to stock potfolio
//    double sum_of_signals;
//    double phi; //trading signal for stock i
//    auto t = clock->current_time();
//
//
//
//
//
//
////compute a sum of exponents of signals, see McFadden choice function...........................(1)
//    for(auto &[k, v] : quotes){
//        const auto &[time, bid,ask] = v.get_price();
//        auto quoted_price_ = (bid + ask)/2.;
//        auto shares_outstanding = v.get_shares_outstanding();
//        phi = shares_outstanding * quoted_price_;//stock signal
//        sum_of_signals = sum_of_signals + abs(phi);
//    }
//
//
//
////compute the ratio of each stock's signal relative to sum described in (1)
//    for(auto &[k, v] : quotes) {
//        //the following line searches whether a trader already holds inventory of the stock in subject
//        // and stores the value in j
//        auto i = stocks_at_hand.find(k);
//        double current_position = 0;
//        if (stocks_at_hand.end() != i) { current_position = i->second; }
//        const auto &[time, bid, ask] = v.get_price();
//        auto shares_outstanding = v.get_shares_outstanding();
//        auto quoted_price_ = (bid + ask) / 2.;
//
//
//        phi = shares_outstanding * quoted_price_;
//        auto allocation = (phi/sum_of_signals);// = phi / sum_of_signals; //stock wealth allocation
//
//
//
//
////investment in the subject stock
//        double excess_demand_market;
//        double excess_demand_limit;
//        order market_order;
//        order limit_order;
//
//
//        int order_num;
//        order_num = concat(t,k);
//        //this->reb_period = 1;
//
//        double free_cash = 0;
//        double short_positions = 0;
//        double long_positions = 0;
//        double investable_wealth = 0;
//
//        for(auto &[i,x]:this->stocks_at_hand){
//            if(x<0){
//                short_positions+= -1 * (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            } else{
//                long_positions += (x) * this->stocks_on_market.find(i)->second.get_midprice();
//            }
//        }
//
//        free_cash =max(0, (this->cash_at_hand - short_positions));
//        investable_wealth = free_cash + long_positions;
//
//        auto disposable_income = max(0,(free_cash - ((1-portfolio_alloc)*investable_wealth)));
//
////rebalancing date
//        if (rebalancing_count.find(k) == rebalancing_count.end()){
//            rebalancing_count.emplace(k,0);
//        }
//
//        if ( rebalancing_count.find(k)->second < floor(t/reb_period)) {
//
//            rebalancing_count.find(k)->second++;
//
//            auto target_position = investable_wealth * allocation * portfolio_alloc/(quoted_price_);
//            excess_demand_market = target_position - current_position;
//
//
//            //check if there are resources to make the desired trade
//            if(0 <= target_position &&  target_position <= current_position){
//                order_num = concat(order_num,1);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size(excess_demand_market, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//
//            }else if((target_position< 0) && (current_position > 0)){
//                //if target is a negative position but we currently have a positive number of stocks,
//                // we first sell the stocks we have then short sell the remainder
//                order_num = concat(order_num,1);
//                market_order.set_id(this->get_identifier());
//                market_order.set_order_type(order::limit);
//                std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
//                market_order.set_ordered_asset(k);
//                market_order.set_order_size(-current_position, quoted_price_);
//                market_order.set_status(order::active);
//                result_.emplace(order_num, market_order);
//
//
//                double remaining_demand = min(abs(target_position*quoted_price_),
//                                              abs(disposable_income*allocation*portfolio_alloc));
//
//                if(remaining_demand>0){
//                    auto demand = (-remaining_demand/(quoted_price_));
//                    order_num = concat(order_num,1);
//                    market_order.set_id(this->get_identifier());
//                    market_order.set_order_type(order::limit);
//                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
//                    market_order.set_ordered_asset(k);
//                    market_order.set_order_size(demand, quoted_price_);
//                    market_order.set_status(order::active);
//                    result_.emplace(order_num, market_order);
//                }
//            }else {
//
//                double remaining_demand = min(abs(target_position*quoted_price_),
//                                              abs(disposable_income*allocation*portfolio_alloc));
//                if(remaining_demand>0){
//                    auto sign = excess_demand_market/abs(excess_demand_market);
//                    auto demand = sign * remaining_demand/quoted_price_;
//                    order_num = concat(order_num,1);
//                    market_order.set_id(this->get_identifier());
//                    market_order.set_order_type(order::limit);
//                    std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
//                    market_order.set_ordered_asset(k);
//                    market_order.set_order_size(demand, quoted_price_);
//                    market_order.set_status(order::active);
//                    result_.emplace(order_num, market_order);
//
//                }
//            }
//
//        }else{
//
//            excess_demand_market = (disposable_income * allocation * portfolio_alloc / (quoted_price_));
//
//            order_num = concat(order_num,1);
//            market_order.set_id(this->get_identifier());
//            market_order.set_order_type(order::limit);
//            std::cout<<this->fund_philosophy<<"threw market order"<<investable_wealth<<std::endl;
//            market_order.set_ordered_asset(k);
//            market_order.set_order_size(excess_demand_market, quoted_price_);
//            market_order.set_status(order::active);
//            result_.emplace(order_num, market_order);
//
//            auto target_position = investable_wealth * allocation * portfolio_alloc / (quoted_price_);
//
//            if ((target_position > current_position) && (current_position <= 0)) {
//                if (target_position < 0) {
//                    auto demand = min(abs(target_position - current_position), abs(current_position));
//                    order_num = concat(order_num, 0);
//                    market_order.set_id(this->get_identifier());
//                    market_order.set_order_type(order::limit);
//                    std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
//                    market_order.set_ordered_asset(k);
//                    market_order.set_order_size(abs(demand), quoted_price_);
//                    market_order.set_status(order::active);
//                    result_.emplace(order_num, market_order);
//                } else {
//                    auto demand = min(abs(target_position - current_position), abs(current_position));
//                    order_num = concat(order_num, 0);
//                    market_order.set_id(this->get_identifier());
//                    market_order.set_order_type(order::limit);
//                    std::cout << this->fund_philosophy << "threw market order" << investable_wealth << std::endl;
//                    market_order.set_ordered_asset(k);
//                    market_order.set_order_size(abs(current_position), quoted_price_);
//                    market_order.set_status(order::active);
//                    result_.emplace(order_num, market_order);
//                }
//            }
//
//        }
//    }
//    return result_;
//}
//



