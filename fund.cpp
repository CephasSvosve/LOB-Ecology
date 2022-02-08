//
// Created by Cephas Svosve on 20/12/2021.
//
#include "fund.h"
#include <math.h>



//.......create unique identifier for each instance
fund::fund(){

}


//metric getters
map<int, double> fund::get_sum_returns(){
    return this->sum_returns;}

map<int, double> fund::get_sum_squared_returns(){
    return this->sum_squared_returns;
}

map<int, double> fund::get_sum_market_by_stock_returns(){
    return this->sum_market_by_stock_returns;
}

map<int, double> fund::get_beta(){
    return this->beta;
}

map<int, double> fund::get_CAPM(){
    return this->CAPM;
}


double fund::get_sum_mkt_return(){
    return this->sum_mkt_return;
}


double fund::get_sum_squared_mkt_return(){
    return this->sum_squared_mkt_return;
}


//metric setters

double fund::get_mkt_return(int t){
    double mkt_price=0;
    double mkt_prev_price=0;
    double mkt_return_;
    double tot_mkt_share;
    for(auto &[k,v]:this->stocks_on_market){
        tot_mkt_share += v.get_shares_outstanding();
    }

    for(auto &[k,v]:this->stocks_on_market){
        double mkt_share = v.get_shares_outstanding()/tot_mkt_share;
        double prev_quote = v.get_hist_price().find(t-2)->second;

        mkt_prev_price += prev_quote *mkt_share;
        mkt_price += ((get<1>(v.get_price())+get<2>(v.get_price()))/2.)*mkt_share;
    }
    mkt_return_ = (mkt_price-mkt_prev_price)/mkt_prev_price;


    this->mkt_return += mkt_return_;

    return mkt_return_;
}


void fund::set_sum_mkt_return(int t){

    double mkt = this->get_mkt_return(t);
    this->sum_mkt_return += mkt;

}


void fund::set_sum_squared_mkt_return(int t){
    this->sum_squared_mkt_return += pow(this->mkt_return,2.);

}


void fund::set_sum_returns(int t){
    for(auto &[k,v] : this->stocks_on_market) {
        double curr_price = v.get_hist_price().find(t-1)->second;
        double prev_price = v.get_hist_price().find(t-2)->second;
        double returns = (curr_price-prev_price)/prev_price;
        auto i = this->sum_returns.find(k);
        if(i != this->sum_returns.end()){
        this->sum_returns.find(k)->second += returns;
        }
        else{
            this->sum_returns.emplace(k,returns);
        }
    }
}


void fund::set_sum_squared_returns(int t){
    for(auto &[k,v] : this->stocks_on_market) {
        double curr_price = v.get_hist_price().find(t-1)->second;
        double prev_price = v.get_hist_price().find(t-2)->second;
        double returns = (curr_price-prev_price)/prev_price;
        auto i = this->sum_squared_returns.find(k);
        if(i != this->sum_squared_returns.end()){
            this->sum_squared_returns.find(k)->second += (returns*returns);
        }
        else{
            this->sum_squared_returns.emplace(k,(returns*returns));
        }
    }

}


void fund::set_sum_market_by_stock_returns(int t){

    for(auto &[k,v]:this->stocks_on_market){
        double curr_price = v.get_hist_price().find(t-1)->second;
        double prev_price = v.get_hist_price().find(t-2)->second;
        double returns = (curr_price-prev_price)/prev_price;
        double mkt_by_stock_returns = this->get_mkt_return(t)*returns;

        auto i = this->sum_market_by_stock_returns.find(k);

        if(i != this->sum_market_by_stock_returns.end()) {
            this->sum_market_by_stock_returns.find(k)->second += mkt_by_stock_returns;
        }else{
            this->sum_market_by_stock_returns.emplace(k,mkt_by_stock_returns);
        }
    }

}


void fund::set_beta(int t){
    for(auto &[k,v]:stocks_on_market){


        double sum_x = this->get_sum_returns().find(k)->second;
        double sum_y = this->get_sum_mkt_return();
        double sum_xx= this->get_sum_squared_returns().find(k)->second;
        double sum_yy= this->get_sum_squared_mkt_return();
        double sum_xy= this->get_sum_market_by_stock_returns().find(k)->second;


        double cov = ((sum_xy/t) - (sum_x*sum_y)/(t*t))*(t/(t-1));
        double var = ((sum_yy) - (sum_y*sum_y)/(t*t))*(t/(t-1));

        double beta_ = cov/var;
        std::cout<<"BBeta "<<beta_ <<std::endl;

        auto i = this->beta.find(k);
        if(i!=this->beta.end()){
            this->beta.find(k)->second = beta_;
        }else
        {
            this->beta.emplace(k,beta_);
        }
    }


}


void fund::set_CAPM(int t){
    double mkt_return_ = this->get_sum_mkt_return()/t;
    std::cout<< " marketret "<<mkt_return_<<std::endl;
    double risk_free_return = pow(1.02,1./252.)-1;
    for(auto &[k,v]:stocks_on_market){
        double beta_ = this->beta.find(k)->second;
        std::cout<< " beta "<<beta_<<std::endl;
        double CAPM_ = risk_free_return + beta_*(mkt_return_- risk_free_return);
        auto i = this->CAPM.find(k);

        if(i != this->CAPM.end()){
            this->CAPM.find(k)->second = CAPM_;
        }else{
            this->CAPM.emplace(k,CAPM_);
        }
    }

}


void fund::set_reb_period(int rebalancing_period){
    this->reb_period = rebalancing_period;
}

void fund::initial_wealth(double cash, double stocks){

    this->cash_at_hand =  cash;
    for(auto &[k,v]:this->stocks_on_market){
        this->stocks_at_hand.emplace(k,stocks);
    }
    balance_cf(0);
}


//......this is a method to sum up all contents of a map (mp) variable
double
fund::current_balance(map<int, double> mp){
    double result;
        for(auto i = mp.begin(); i!= mp.end(); i++ ){
            result = result + i->second;}
                return result;
}




void
fund::double_entry(map<int, double> &debit
        , map<int, double> &credit
            , map<int, double> amount){

                for(auto i : amount){
                     debit.find(i.first)->second  += i.second;
                        credit.find(i.first)->second -= i.second;
    }
}


void
fund::balance_cf(int t){
//calculate wealth just before executing orders to account for
    // wealth distribution changes due to interest and dividends;
    //......update wealth
        auto aq = this->wealth;
    this->wealth = this->cash_at_hand;
    for(auto &[x,v] : this->stocks_at_hand){
        auto quote = (get<1>(this->stocks_on_market.find(x)->second.get_price())+
                      get<2>(this->stocks_on_market.find(x)->second.get_price()))/2.;
        auto value_of_stocks  = quote * v;
        double total = value_of_stocks;
        this->wealth += total;
    }
    auto ar = this->wealth;
        //......we start by adding interest to the cash at hand, assuming all cash is invested in a govt bond
            double interest_rate = pow(1.02,1/252)-1;
                this->interest_on_cash = interest_rate * this->cash_at_hand;
                    this->cash_at_hand += interest_on_cash;
                        this->wealth += interest_on_cash;

    auto as = this->wealth;
//......we compute dividends income as the sum of dividends received from all stocks at hand
double div_ ;
    for(auto &[i,v] : this->stocks_at_hand){
        auto dividend = this->stocks_on_market.find(i)->second.get_dividend(0);//TODO correct dividend time to t
            div_ += dividend * v;}

    if(this->dividends_received.size()<10){
        this->dividends_received.emplace(this->clock.current_time(),div_);
    }else{
        this->dividends_received.erase(this->dividends_received.begin());
        this->dividends_received.emplace(this->clock.current_time(),div_);
    }

                    this->cash_at_hand += div_;
                        this->wealth += div_;

                        auto at = this->wealth;


}



void
fund::balance_bd(int t,vector<order> &executed_orders){
//we update records from the balance carried forward to account
    // for wealth distribution changes due to trade executions



//......update cash at hand
    for(auto &x : executed_orders){
        auto cost  = x.get_proposed_price() * x.get_order_size();
        this->cash_at_hand -= cost;
    }

//.....update stocks at hand (*NB negative number of stocks at hand represent short positions)
            for(auto &x : executed_orders){
                auto y = this->stocks_at_hand.find(x.get_ordered_asset());
                    if (y != this->stocks_at_hand.end()){
                        this->stocks_at_hand.find(x.get_ordered_asset())->second += x.get_order_size();
                            }else {
                        this->stocks_at_hand.emplace(x.get_ordered_asset(),x.get_order_size());
            }
        }




//.....we the update short-stock record
    for(auto &x : stocks_at_hand) {
        if(x.second< 0.){
            auto y = borrowed_stocks.find(x.first);
                if (y != borrowed_stocks.end()){
                    y->second = x.second;
                        }else{
                            borrowed_stocks.emplace(x);
            }
        }
}



//......update equities, the monetary value of stocks at hand
    for(auto x : stocks_at_hand) {
        auto y = equities.find(x.first);
            //to get price we look for stock on the market data and find its price
                auto quote = get<1>(stocks_on_market.find(x.first)->second.get_price());
                    auto equities_value = quote * x.second;
                    if ( y!= equities.end()){
                        y->second = equities_value;
                            }else {
                                equities.emplace(x.first,equities_value);
                    }
}



//......update loans received, monetary value of short positions
    for(auto &x : borrowed_stocks){
        auto y = loans_received.find(x.first);
            double quote = get<1>(stocks_on_market.find(x.first)->second.get_price());
            if(y!= loans_received.end()){
                loans_received.find(x.first)->second = quote * x.second;
                        } else{
                            loans_received.emplace(x.first, quote * x.second);
                }
    }



//TODO look into transaction cost more carefully
    //map<int, price>  transaction_cost;





//update metrics
//if(stocks_on_market.find(1)->second.get_hist_price().size() > 2){
//    this->set_sum_mkt_return(t);
//    this->set_sum_squared_mkt_return(t);
//    this->set_sum_returns(t);
//    this->set_sum_squared_returns(t);
//    this->set_sum_market_by_stock_returns(t);
//    this->set_beta(t);
//    this->set_CAPM(t);
//    }
}

void
fund::balance(){
//this is a variant of the balance_bd method but specifically for the market maker alone
// who balances after taking care of the executed orders


//.....we the update short-stock record
    for(auto &x : stocks_at_hand) {
        if(x.second< 0.){
            auto y = borrowed_stocks.find(x.first);
            if (y != borrowed_stocks.end()){
                y->second = x.second;
            }else{
                borrowed_stocks.emplace(x);
            }
        }
    }



//......update equities, the monetary value of stocks at hand
    for(auto &x : stocks_at_hand) {
        auto y = equities.find(x.first);
        //to get price we look for stock on the market data and find its price
        auto quote = get<1>(stocks_on_market.find(x.first)->second.get_price());
        auto equities_value = quote * x.second;
        if ( y!= equities.end()){
            equities.find(x.first)->second = equities_value;
        }else {
            equities.emplace(x.first,equities_value);
        }
    }



//......update loans received, monetary value of short positions
    for(auto &x : borrowed_stocks){
        auto y = loans_received.find(x.first);
        double quote = get<1>(stocks_on_market.find(x.first)->second.get_price());
        if(y!= loans_received.end()){
            loans_received.find(x.first)->second = quote * x.second;
        } else{
            loans_received.emplace(x.first, quote * x.second);
        }
    }



//TODO look into transaction cost more carefully
    //map<int, price>  transaction_cost;



}


//method for accessing stocks data
void
fund::get_inform(map<int, company> &assets){
    for(auto &[i,v] : assets){
        this->stocks_on_market.emplace(i,v);
    }

}

