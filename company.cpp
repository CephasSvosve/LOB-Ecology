
#include "company.h"
#include <math.h>
#include <iostream>
#include <utility>


company::company(){
    this->set_identifier(create_identity(company_));
    std::cout<<"id "<<this->identifier;
}
template <typename T, typename H>
void company::stack_mp(std::map<T,H> &mp,int maxsize, T indx, H val){

    if(mp.size()<maxsize){
        mp.emplace(indx,val);
    }else{
        mp.erase(mp.begin());
        mp.emplace(indx,val);
    }

}

//returns the most recent paid-out dividend
double
company::get_dividend(int t){
    auto dividend_ = (this->dividends[floor(t/63)]);
    return dividend_;
}



//returns the most recent reported earnings
double
company::get_earnings(int t) {
    auto earnings_ = this->earnings.at(floor(t));
    return earnings_;
}


//returns the most recent reported free_cash_flow
double
company::get_free_cash_flow(){
    auto free_cash_flow_ = this->free_cash_flow[floor(clock.current_time())];
    return free_cash_flow_;
}
int
company::get_identifier() const{
    return this->identifier;
}


vector<double>
company::get_hist_price(){
    return this->historical_price;
}

int
company::get_shares_outstanding() const{
    return this->shares_outstanding;
}


map<int, double>
company::get_market_capital(){
    return this->market_capital;
}


//updates the market price of a stock
void
company::set_price(int t, double bid, double ask){

    double clearing_price = (bid + ask) / 2.;
    this->bid_price = bid;
    this->ask_price = ask;
    this->market_price = make_tuple(t, clearing_price);

    if(historical_price.size()<250){
        historical_price.push_back(clearing_price);
    }else{
        historical_price.erase(historical_price.begin());
        historical_price.push_back(clearing_price);
    }
    //stack_mp(historical_price,250,t,clearing_price);

}

void
company::set_identifier(int identifier_){
    this->identifier = identifier_;
}
void
company::set_shares_outstanding(int outstanding_shares){
    this->shares_outstanding = outstanding_shares;
}


void
company::set_market_capital(){
    double market_cap = get<1>(this->get_price()) * this->get_shares_outstanding();

    if(market_capital.size()<250){
        market_capital.emplace(int(this->clock.current_time()),market_cap);
    }else{
        market_capital.erase(market_capital.begin());
        market_capital.emplace(int(this->clock.current_time()),market_cap);
    }

}

//returns current price of a stock
tuple<int, double,double>
company::get_price(){
    tuple<int, double,double> result;
    get<0>(result) = get<0>(this->market_price);
    get<1>(result) = this->bid_price;
    get<2>(result) = this->ask_price;

    return result;

}

double
company::get_midprice(){
    double result = 0;
    double a = this->bid_price;
    double b = this->ask_price;
    result = (a+b)/2.;

    return result;

}

// a method overload which returns price at a given time in history
vector<double>
company::get_price_range(int range){
    vector<double> result;
    for(int i = 0; i < range; i++){
        if(i<this->historical_price.size()) {
            result.push_back(this->historical_price.rbegin()[i]);
        } else{
            break;
        }
    }
    return result;
}

void company::set_earnings_process(vector<double> ep){
    this->earnings = std::move(ep);
}


void company::set_dividends_process(vector<double> dp){
    this->dividends = dp;
}

void company::set_free_cash_flow_process(vector<double> fcfp){
    this->free_cash_flow = fcfp;
}


void company::set_clock(market_watch &timer){
    this->clock = timer;
}
