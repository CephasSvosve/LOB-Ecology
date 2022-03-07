//
// Created by Cephas Svosve on 20/12/2021.
//

#ifndef UNTITLED30_COMPANY_H
#define UNTITLED30_COMPANY_H
#include "market_watch.h"
#include <map>
#include <vector>
#include "identity.h"


using namespace std;

class company{
private:
    vector<double> earnings;
    vector<double> dividends;
    vector<double> free_cash_flow;
    tuple<int, double> market_price;//int-market_watch
    vector<double> historical_price;//int-market_watch
    double bid_price;
    double ask_price;

    //properties of a company
    int identifier= 0;
    int shares_outstanding=0; //we assume this value is constant throughout simulations
    map<int, double> market_capital;//int-market_watch


public:
company();

market_watch clock;
template <typename T, typename H>
void stack_mp(std::map<T,H> &mp,int maxsize, T indx, H val);


//setters
void
set_price(int t, double bid, double ask);

void
set_identifier(int identifier);

void
set_shares_outstanding(int outstanding_shares);

void
set_market_capital();

void
set_earnings_process(vector<double> ep);

void
set_dividends_process(vector<double> dp);

void
set_free_cash_flow_process(vector<double> dp);




//getters
double get_dividend(int t);

double get_earnings(int t);

double get_free_cash_flow();

[[nodiscard]] int
get_identifier() const;

[[nodiscard]] int
get_shares_outstanding() const;

map<int, double> get_market_capital();

tuple<int, double, double> get_price();

double get_midprice();

vector<double> get_price_range(int range);

vector<double> get_hist_price();

 void set_clock(market_watch &timer);

};
#endif //UNTITLED30_COMPANY_H
