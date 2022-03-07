//
// Created by Cephas Svosve on 20/12/2021.
//

#ifndef UNTITLED30_FUND_H
#define UNTITLED30_FUND_H

#include "company.h"
#include "order.h"
#include <iostream>
#include "market_watch.h"

class fund
        : public company
        {
public:
    fund();
    int reb_period;

    //Balance Sheet *NB stocks are identified by an integer variable on the map variables
    //Assets
    map<int, company>  stocks_on_market;//int-key,|v



    double  wealth;
    double  cash_at_hand;
    double  interest_on_cash;
    double target_price_ = 0;
    map<int, double>  debtors;
    map<int, double>  dividends_received;
    map<int, double>  stocks_at_hand;
    map<int, double>  equities;//derived from number_of_shares

    //Liabilities
    double interest_on_loans;
    map<int, double>  borrowed_stocks;
    map<int, double>  loans_received;//int-key|v
    map<int, double>  transaction_cost;
    tuple<int, double> market_price ;//int-market_watch
    map<int, double> market_capital;//int-market_watch

    //metrics
    map<int, double>  sum_returns;//sum of returns
    map<int, double>  sum_squared_returns;//sum of squared returns
    map<int, double>  sum_market_by_stock_returns;//sum of (return*market return)

    //derived metrics
    map<int, double>  beta;
    map<int, double>  CAPM;

    //overall market's metrics
    double mkt_return=0;
    double sum_mkt_return=0;
    double sum_squared_mkt_return=0;


static void
   double_entry(map<int, double> &debit
                , map<int, double> &credit
                    , const map<int, double>& amount);


//metric getters

    map<int, double> get_sum_returns();
    map<int, double> get_sum_squared_returns();
    map<int, double> get_sum_market_by_stock_returns();
    map<int, double> get_beta();
    map<int, double> get_CAPM();

    double get_sum_mkt_return();
    double get_sum_squared_mkt_return();

//metric setters
    double get_mkt_return(int t);
    void set_sum_mkt_return(int t);
    void set_sum_squared_mkt_return(int t);

    void set_sum_returns(int t);
    void set_sum_squared_returns(int t);
    void set_sum_market_by_stock_returns(int t);
    void set_beta(int t);
    void set_CAPM(int t);




void set_reb_period(int rebalancing_period);
void initial_wealth(double cash, double stocks);

 //balancesheet updaters
double
current_balance(const map<int, double>& mp);

void
balance_cf(int t);

void
balance_bd(int t, vector<order> &executed_orders);

void
balance();

void get_inform(map<int, company> &assets);



};
#endif //UNTITLED30_FUND_H
