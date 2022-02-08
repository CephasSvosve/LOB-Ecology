//
// Created by Cephas Svosve on 22/12/2021.
//

#include "funds.h"

funds::funds() {

}




map<int, order>
        funds::invest(){
    map<int, order> result;
    switch (this->fund_philosophy) {

        case 0:
            result= this->value_demand();
            break;
//        case 1:
//           std::cout<<"growth"<<std::endl;
//            result= this->growth_demand();
//            break;
        case 1:
            result= this->noise_demand();
            break;
        case 2:
            if(this->clock.current_time()>1000 /*window_size*/){
            result= this->momentum_demand();}
            break;

 //       case 6:
//            break;//TODO ADD INDEX DEMAND FUNCTION
        default:
            break;

    }
    return result;

}


void funds::trade_strategy(trading_strategy tradingStrategy){
    this->fund_philosophy = tradingStrategy;
}

double funds::compute_earnings(company &stock, int &time, int &reb_period){


    //earnings

    auto D_o = stock.get_dividend(0);//earnings at time 0
    auto D_t = stock.get_dividend(time);//recent earnings
    //auto t_reb = time - reb_period*floor(time/reb_period);//number of ticks past rebalancing day
    double intrinsic_value;
    if(time > 0){
        double growth = D_t/D_o;
        double g_rate = pow(growth,(1./time))-1;//earnings historical growth rate
        double r_free_rate = pow(1.02,1./252.)-1;
        double capm = this->get_CAPM().find(stock.get_identifier())->second;
//        std::cout<<"time tick "<<capm <<std::endl;
        if(isnan(capm)){
            intrinsic_value =D_t/0.02;
        }else{
        intrinsic_value = D_t*(1+g_rate)/(capm - g_rate);
        }
        }
    else{
        intrinsic_value = D_o /0.02;
    }

    return intrinsic_value;
}


std::map<int,order>
funds::value_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    double Beta = 9;//dummy threshold TODO: link this variable with actual earnings threshold
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto tick = clock.current_time();






//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid,ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto earnings = v.get_earnings(tick);
        phi = -(log((quoted_price_)/earnings)-log(Beta));//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes){
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double j = 0;
        if (stocks_at_hand.end() != i){j = i->second;}
        const auto &[time, bid,ask] = v.get_price();
        auto earnings = v.get_earnings(tick);
        auto quoted_price_ = (bid+ask)/2.;


        phi = -(log((quoted_price_)/earnings)-log(Beta));
        auto allocation = phi/sum_of_signals; //stock wealth allocation




//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


       // std::cout<<"earnings "<<earnings<<std::endl;
        int order_num;
        order_num = (10*tick) +k;
        //this->reb_period = 1;


//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if(tick == this->reb_period*floor(tick/this->reb_period)){

            excess_demand_market = (this->wealth * allocation * portfolio_alloc/(quoted_price_))-j;
            if (excess_demand_market < 0) {
                market_order.set_id(this->get_identifier());
                // std::cout<<"value identifier "<< this->get_identifier()<<std::endl;
                std::cout << tick << "value identifier cash " << this->cash_at_hand << " *j* " << j << " excessdem "
                          << excess_demand_market << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if (excess_demand_market > 0){
                bool liquid = (this->cash_at_hand > (excess_demand_market*quoted_price_));

                if(liquid){
                    market_order.set_id(this->get_identifier());
                    // std::cout<<"value identifier "<< this->get_identifier()<<std::endl;
                    std::cout << tick << "value identifier cash " << this->cash_at_hand << " *j* " << j << " excessdem "
                              << excess_demand_market << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
                    market_order.set_order_type(order::market);
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(excess_demand_market, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
                else{

                    auto demand = (this->cash_at_hand * allocation * portfolio_alloc/(quoted_price_));
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }


            order_num = (10 * (tick + 1)) + k;
            excess_demand_limit = -this->stocks_at_hand.find(k)->second;
            limit_order.set_id(this->get_identifier());
            limit_order.set_order_type(order::limit);
            limit_order.set_ordered_asset(k);
            auto target_price = Beta * earnings;
            limit_order.set_order_size(excess_demand_limit, target_price);
            limit_order.set_status(order::active);
            result_.emplace(order_num, limit_order);

        }else

        {
            auto disposable_income = this->cash_at_hand - (1-portfolio_alloc)*this->wealth;
            if(disposable_income > 0 && this->cash_at_hand > 0){
            excess_demand_market = (disposable_income * allocation * portfolio_alloc/(quoted_price_))-j;
            }
            order_num = (10 * (tick + 1)) + k;
            std::cout<<tick<<"value identifier cash "<<this->cash_at_hand<<" *j* "<<j <<" excessdem "<<excess_demand_market
                    << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
            market_order.set_id(this->get_identifier());
            market_order.set_order_type(order::market);
            market_order.set_ordered_asset(k);
            market_order.set_order_size(excess_demand_market,quoted_price_);
            market_order.set_status(order::active);
            result_.emplace(order_num, market_order);


            order_num = (10*(tick+1)) + k;
            excess_demand_limit = -this->stocks_at_hand.find(k)->second;
            limit_order.set_id(this->get_identifier());
            limit_order.set_order_type(order::limit);
            limit_order.set_ordered_asset(k);
            auto target_price = Beta * earnings;
            limit_order.set_order_size(excess_demand_limit, target_price);
            limit_order.set_status(order::active);
            result_.emplace(order_num, limit_order);
        }


    }
    return result_;
}


std::map<int,order>
funds::growth_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    double Beta = 9;//dummy threshold TODO: link this variable with actual earnings threshold
    double portfolio_alloc = 0.6;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto tick = clock.current_time();




//compute a sum of exponents of signals, see McFadden choice function...........................(1)
    for(auto &[k, v] : quotes){
        const auto &[time, bid,ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto earnings = v.get_earnings(tick);
        phi = (log((quoted_price_)/earnings)-log(Beta));//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes){
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double j = 0;
        if (stocks_at_hand.end() != i){j = i->second;}
        const auto &[time, bid, ask] = v.get_price();
        auto quoted_price_ = (bid + ask)/2.;
        auto earnings = v.get_earnings(tick);

        phi =(log((quoted_price_)/earnings)-log(Beta));
        auto allocation = phi/sum_of_signals; //stock wealth allocation




//investment in the subject stock
        double excess_demand_market;
        double excess_demand_limit;
        order market_order;
        order limit_order;


        // std::cout<<"earnings "<<earnings<<std::endl;
        int order_num;
        order_num = (10*tick) +k;
        //this->reb_period = 1;


//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if(tick == this->reb_period*floor(tick/this->reb_period)){

            excess_demand_market = (this->wealth * allocation * portfolio_alloc/(quoted_price_))-j;
            if (excess_demand_market < 0) {
                market_order.set_id(this->get_identifier());
                // std::cout<<"value identifier "<< this->get_identifier()<<std::endl;
                std::cout << tick << "growth identifier cash " << this->cash_at_hand << " *j* " << j << " excessdem "
                          << excess_demand_market << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if (excess_demand_market > 0){
                bool liquid = (this->cash_at_hand > (excess_demand_market*quoted_price_));

                if(liquid){
                    market_order.set_id(this->get_identifier());
                    // std::cout<<"value identifier "<< this->get_identifier()<<std::endl;
                    std::cout << tick << "growth identifier cash " << this->cash_at_hand << " *j* " << j << " excessdem "
                              << excess_demand_market << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
                    market_order.set_order_type(order::market);
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(excess_demand_market, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
                else if(this->cash_at_hand > 0){

                    auto demand = (this->cash_at_hand * allocation * portfolio_alloc/(quoted_price_));
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }


            order_num = (10 * (tick + 1)) + k;
            excess_demand_limit = -this->stocks_at_hand.find(k)->second;
            limit_order.set_id(this->get_identifier());
            limit_order.set_order_type(order::limit);
            limit_order.set_ordered_asset(k);
            auto target_price = Beta * earnings;
            limit_order.set_order_size(excess_demand_limit, target_price);
            limit_order.set_status(order::active);
            result_.emplace(order_num, limit_order);

        }else

        {
            auto disposable_income = this->cash_at_hand - (1-portfolio_alloc)*this->wealth;
            if(disposable_income > 0 && this->cash_at_hand > 0){
                excess_demand_market = (disposable_income * allocation * portfolio_alloc/(quoted_price_))-j;
            }
            order_num = (10 * (tick + 1)) + k;
            std::cout<<tick<<"growth identifier cash "<<this->cash_at_hand<<" *j* "<<j <<" excessdem "<<excess_demand_market
                    << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
            market_order.set_id(this->get_identifier());
            market_order.set_order_type(order::market);
            market_order.set_ordered_asset(k);
            market_order.set_order_size(excess_demand_market,quoted_price_);
            market_order.set_status(order::active);
            result_.emplace(order_num, market_order);


            order_num = (10*(tick+1)) + k;
            excess_demand_limit = -this->stocks_at_hand.find(k)->second;
            limit_order.set_id(this->get_identifier());
            limit_order.set_order_type(order::limit);
            limit_order.set_ordered_asset(k);
            auto target_price = Beta * earnings;
            limit_order.set_order_size(excess_demand_limit, target_price);
            limit_order.set_status(order::active);
            result_.emplace(order_num, limit_order);
        }


    }
    return result_;
}

std::map<int,order>
funds::momentum_demand() {

//Params definition
    map<int, order> result_;
    auto quotes = this->stocks_on_market;
    double Beta = 9;//dummy threshold TODO: link this variable with actual earnings
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
            auto MA1 = current_balance(v.get_price_range(window_size_MA1 - 1));
            auto MA2 = current_balance(v.get_price_range(window_size_MA2 - 1));
            double trend1 = (MA1 + (quoted_price_)) / window_size_MA1;
            double trend2 = (MA2 + (quoted_price_)) / window_size_MA2;
            if (MA2> 0.) { phi = ((trend1 / trend2) - 1); }
            else { phi = 0.; }//stock signal
            sum_of_signals = sum_of_signals + abs(phi);
        }
    }



//compute the ratio of each stock's signal relative to sum described in (1)
    for (auto &[k, v] : quotes) {
        double j = 0;
        auto i = stocks_at_hand.find(k);
        if (stocks_at_hand.end() != i) { j = i->second; }

            auto quoted_price_ = (get<1>(v.get_price())+get<2>(v.get_price()))/2.;
            int window_size_MA1 = 50;
            int window_size_MA2 = 200;
            auto MA1 = current_balance(v.get_price_range(window_size_MA1 - 1));
            auto MA2 = current_balance(v.get_price_range(window_size_MA2 - 1));
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
            auto order_num = (10*t) + k;


//......we compute a market order which represents an instantaneous trader's entry into a position given
            // the day is a rebalancing date or the fund has received new investment wealth inflow.

            if(t == this->reb_period*floor(t/this->reb_period) && t > window_size){

                excess_demand_market = (this->wealth * allocation * portfolio_alloc/(quoted_price_))-j;
                if (excess_demand_market < 0) {
                    market_order.set_id(this->get_identifier());
                    // std::cout<<"value identifier "<< this->get_identifier()<<std::endl;
                    std::cout << t << "trend identifier cash " << this->cash_at_hand << " *j* " << j << " excessdem "
                              << excess_demand_market << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
                    market_order.set_order_type(order::market);
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(excess_demand_market, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);

                }else if (excess_demand_market > 0){
                    bool liquid = (this->cash_at_hand > (excess_demand_market*quoted_price_));

                    if(liquid){
                        market_order.set_id(this->get_identifier());
                        // std::cout<<"value identifier "<< this->get_identifier()<<std::endl;
                        std::cout << t << "trend identifier cash " << this->cash_at_hand << " *j* " << j << " excessdem "
                                  << excess_demand_market << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
                        market_order.set_order_type(order::market);
                        market_order.set_ordered_asset(k);
                        market_order.set_order_size(excess_demand_market, quoted_price_);
                        market_order.set_status(order::active);
                        result_.emplace(order_num, market_order);
                    }
                    else if(this->cash_at_hand > 0){

                        auto demand = (this->cash_at_hand * allocation * portfolio_alloc/(quoted_price_));
                        market_order.set_id(this->get_identifier());
                        market_order.set_order_type(order::market);
                        market_order.set_ordered_asset(k);
                        market_order.set_order_size(demand, quoted_price_);
                        market_order.set_status(order::active);
                        result_.emplace(order_num, market_order);
                    }
                }


            }if( t > window_size)

            {
                auto disposable_income = this->cash_at_hand - (1-portfolio_alloc)*this->wealth;
                if(disposable_income > 0 && this->cash_at_hand > 0){
                    excess_demand_market = (disposable_income * allocation * portfolio_alloc/(quoted_price_))-j;
                }
                order_num = (10 * (t + 1)) + k;
                std::cout<<t<<"trend identifier cash "<<this->cash_at_hand<<" *j* "<<j <<" excessdem "<<excess_demand_market
                        << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
                market_order.set_id(this->get_identifier());
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market,quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }


        }
        return result_;
    }


MatrixXd funds::generateWhiteNoise1(int seed,int rows, int columns){
    MatrixXd randoms(rows, columns);
    VectorXd a;
    int i =0;

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
funds::noise_demand(){
    //Params definition
    map<int,order> result_;
    auto quotes = this->stocks_on_market;
    double Beta = 9;//dummy threshold TODO: link this variable with actual earnings
    double portfolio_alloc = 1;// wealth allocated to stock potfolio
    double sum_of_signals;
    double phi; //trading signal for stock i
    auto t = clock.current_time();



//Ornstein Uhlenbeck params
    MatrixXd dX = generateWhiteNoise1(6,3,1000);

    //noise with a half life of 6 years, to match empirical evedence
    double mean_reversion_rate = 1 - pow(0.5, 1/(6.0*252.0));
    double sigma = .12;

    //Ornstein Uhlenbeck noise generation
    double  X_t ;
    int i =0;
    int c =0;
    map<int, vector<double>> noise;


//    for(auto &[k, v] : quotes) {
//
//        if (quotes.size() == dX.rows()) {
//            noise.emplace(k, dX(c % 3));
//
//            c++;
//        }
//    }

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
        //std::cout<<"noise *"<<X_t<<std::endl;
        phi = -(log((quoted_price_)/earnings)-log(Beta)) + X_t;//stock signal
        sum_of_signals = sum_of_signals + abs(phi);
        i++;

    }


c=0;
//compute the ratio of each stock's signal relative to sum described in (1)
    for(auto &[k, v] : quotes){
        //the following line searches whether a trader already holds inventory of the stock in subject
        // and stores the value in j
        auto i = stocks_at_hand.find(k);
        double j = 0;
        if (stocks_at_hand.end() != i){j = i->second;}
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
        auto order_num = (10*t) + k;



//......we compute a market order which represents an instantaneous trader's entry into a position given
        // the day is a rebalancing date or the fund has received new investment wealth inflow.
        if(t == this->reb_period*floor(t/this->reb_period)){

            excess_demand_market = (this->wealth * allocation * portfolio_alloc/(quoted_price_))-j;
            if (excess_demand_market < 0) {
                market_order.set_id(this->get_identifier());
                std::cout << t << "noise identifier cash " << this->cash_at_hand << " *j* " << j << " excessdem "
                          << excess_demand_market << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
                market_order.set_order_type(order::market);
                market_order.set_ordered_asset(k);
                market_order.set_order_size(excess_demand_market, quoted_price_);
                market_order.set_status(order::active);
                result_.emplace(order_num, market_order);

            }else if (excess_demand_market > 0){
                bool liquid = (this->cash_at_hand > (excess_demand_market*quoted_price_));

                if(liquid){
                    market_order.set_id(this->get_identifier());
                    // std::cout<<"value identifier "<< this->get_identifier()<<std::endl;
                    std::cout << t << "noise identifier cash " << this->cash_at_hand << " *j* " << j << " excessdem "
                              << excess_demand_market << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
                    market_order.set_order_type(order::market);
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(excess_demand_market, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
                else{

                    auto demand = (this->cash_at_hand * allocation * portfolio_alloc/(quoted_price_));
                    market_order.set_id(this->get_identifier());
                    market_order.set_order_type(order::market);
                    market_order.set_ordered_asset(k);
                    market_order.set_order_size(demand, quoted_price_);
                    market_order.set_status(order::active);
                    result_.emplace(order_num, market_order);
                }
            }


            order_num = (10 * (t + 1)) + k;
            excess_demand_limit = -this->stocks_at_hand.find(k)->second;
            limit_order.set_id(this->get_identifier());
            limit_order.set_order_type(order::limit);
            limit_order.set_ordered_asset(k);
            auto target_price = Beta * earnings;
            limit_order.set_order_size(excess_demand_limit, target_price);
            limit_order.set_status(order::active);
            result_.emplace(order_num, limit_order);

        }else

        {
            auto disposable_income = this->cash_at_hand - (1-portfolio_alloc)*this->wealth;
            if(disposable_income > 0 && this->cash_at_hand > 0){
                excess_demand_market = (disposable_income * allocation * portfolio_alloc/(quoted_price_))-j;
            }
            order_num = (10 * (t + 1)) + k;
            std::cout<<t<<"noise identifier cash "<<this->cash_at_hand<<" *j* "<<j <<" excessdem "<<excess_demand_market
                    << " stocksathand "<<this->stocks_at_hand.find(k)->second<<std::endl;
            market_order.set_id(this->get_identifier());
            market_order.set_order_type(order::market);
            market_order.set_ordered_asset(k);
            market_order.set_order_size(excess_demand_market,quoted_price_);
            market_order.set_status(order::active);
            result_.emplace(order_num, market_order);


            order_num = (10*(t+1)) + k;
            excess_demand_limit = -this->stocks_at_hand.find(k)->second;
            limit_order.set_id(this->get_identifier());
            limit_order.set_order_type(order::limit);
            limit_order.set_ordered_asset(k);
            auto target_price = Beta * earnings;
            limit_order.set_order_size(excess_demand_limit, target_price);
            limit_order.set_status(order::active);
            result_.emplace(order_num, limit_order);
        }


    }
    return result_;
}



//std::map<int,order>
//funds::mm_demand(){
//    //Params definition
//    map<int,order> result_;
//    auto quotes = this->stocks_on_market;
//    auto tick = clock.current_time();
//
////compute the ratio of each stock's signal relative to sum described in (1)
//    for(auto &[k, v] : quotes){
//
//        const auto &[time, bid, ask] = v.get_price();
//        auto quoted_price_ = (bid + ask)/2.;
//
////investment in the subject stock
//        double excess_demand_market;
//        double excess_demand_limit;
//        order market_order;
//        order limit_order;
//
//
//        // std::cout<<"earnings "<<earnings<<std::endl;
//        int order_num;
//        order_num = (10*tick) +k;
//        //this->reb_period = 1;
//
//
////......we compute a market order which represents an instantaneous trader's entry into a position given
//        // the day is a rebalancing date or the fund has received new investment wealth inflow.
//
//            excess_demand_market = 1'000;
//            market_order.set_id(this->get_identifier());
//            market_order.set_order_type(order::market);
//            market_order.set_ordered_asset(k);
//            market_order.set_order_size(excess_demand_market,quoted_price_);
//            market_order.set_status(order::active);
//            result_.emplace(order_num, market_order);
//            order_num = (10*(tick+1)) + k;
//
//            excess_demand_limit = -this->stocks_at_hand.find(k)->second;
//            limit_order.set_id(this->get_identifier());
//            limit_order.set_order_type(order::market);
//            limit_order.set_ordered_asset(k);
//            limit_order.set_order_size(excess_demand_limit, quoted_price_);
//            limit_order.set_status(order::active);
//            result_.emplace(order_num, limit_order);
//
//
//
//
//    }
//    return result_;
//}