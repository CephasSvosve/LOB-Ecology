


#include <iostream>
#include "company.h"
#include "funds.h"
#include "price.h"
#include "order.h"
#include "price_setter.h"
#include "market_watch.h"
#include <fstream>



int main() {
market_watch watch(market_watch::one_second);


//set overal wealth on the market

double market_wealth = 1.68e10;
double shares_outstanding = 0.5*market_wealth/9;
double market_cash = 0.5*market_wealth;

//Define the price setter
price_setter market;
market.set_clock(&watch);



//create assets and register stocks to the market
market.register_assets(3);
market.generate_fundamentals();


//generate initial quotes
market.set_initial_quotes();

//set price setter's wealth
market.initial_wealth(0.5*market_cash, 0.5*shares_outstanding);



//initialize traders
funds
    aggressive_value
//            ,averse_value
                ,aggressive_noisy
//                    , averse_noisy
//                        ,aggressive_momentum
//                            , averse_momentum
//                                ,aggressive_growth
//                                    , averse_growth
                ;


vector<funds*> participants = {
                                &aggressive_value
//                                    ,&averse_value
                                         ,&aggressive_noisy
//                                            , &averse_noisy
//                                                ,&aggressive_momentum
//                                                        ,&averse_momentum
//                                                             ,&aggressive_growth
//                                                                , &averse_growth
                                            };



//number of participants
double investor_population = participants.size();
//inform traders about assets on the market
    for(auto &p : participants){
        p->set_clock(&watch);
        p->get_inform(market.tradeable_assets());
    }


    aggressive_value.trade_strategy(funds::aggressive_fundamental_value);
    aggressive_value.set_reb_period(63);
    aggressive_value.initial_wealth((1/(investor_population))*0.5*market_cash
                                    ,(1/(investor_population))*0.5*shares_outstanding);


//    averse_value.trade_strategy(funds::fundamental_value);
//    averse_value.set_reb_period(63);
//    averse_value.initial_wealth((1./(investor_population))*0.5*market_cash
//            ,(1./(investor_population))*0.5*shares_outstanding);

//    aggressive_growth.trade_strategy(funds::aggressive_growth);
//    aggressive_growth.set_reb_period(63);
//    aggressive_growth.initial_wealth((1./(investor_population))*0.5*market_cash
//            ,(1./(investor_population))*0.5*shares_outstanding);

//    averse_growth.trade_strategy(funds::growth);
//    averse_growth.set_reb_period(254);
//    averse_growth.initial_wealth((1./(investor_population))*0.5*market_cash
//            ,(1./(investor_population))*0.5*shares_outstanding);

////
    aggressive_noisy.trade_strategy(funds::aggressive_noise);
    aggressive_noisy.set_reb_period(1);
    aggressive_noisy.initial_wealth((1/(investor_population))*0.5*market_cash
            ,(1/(investor_population))*0.5*shares_outstanding);


//    averse_noisy.trade_strategy(funds::noise);
//    averse_noisy.set_reb_period(1);
//    averse_noisy.initial_wealth((1./(investor_population))*0.5*market_cash
//            ,(1./(investor_population))*0.5*shares_outstanding);

////
//    aggressive_momentum.trade_strategy(funds::aggressive_momentum_investment);
//    aggressive_momentum.set_reb_period(21);
//    aggressive_momentum.initial_wealth((1./(investor_population))*0.5*market_cash
//            ,(1./(investor_population))*0.5*shares_outstanding);

//    averse_momentum.trade_strategy(funds::momentum_investment);
//    averse_momentum.set_reb_period(43);
//    averse_momentum.initial_wealth((1./(investor_population))*0.5*market_cash
//            ,(1./(investor_population))*0.5*shares_outstanding);
//


//register traders on the market
market.register_traders(participants);


//market simulation point
std::ofstream myfile;

myfile.open ("NewMarketEcology.csv");

watch.start();
int day_count = 0;
while((&watch)->current_time() < 10000) {

    //adjust traders' balances for interest and dividends


    if (day_count < int((&watch)->current_time())) {
        for (auto &p : participants) {
            p->balance_cf((&watch)->current_time());
        }
    }



//receive orders from traders

    market.receive_orders();

//clear market



//std::cout<< watch.current_time() << ". ";



    if (day_count < int((&watch)->current_time())) {
        for (auto &[k, v] : market.tradeable_assets()) {
            myfile << k << "," << (v.get_midprice()) << ",";

            myfile << (&watch)->current_time() << "," << v.get_dividend((&watch)->current_time()) << ",";
            myfile << (&watch)->current_time() << "," << v.get_value((&watch)->current_time()) << ",";
            if (k == 1) {
                double st = 0;

                for (auto &[x, y]:market.trading_institutions) {
                    myfile << y->wealth << ",";
                    st += y->stocks_at_hand.find(1)->second;
                    std::cout << y->get_identifier() << " " << y->wealth << ",";
                }
                myfile << market.wealth << ",";
                myfile << (v.get_midprice()) << ",";
            }

            //myfile<<(get<1>(v.get_price())+get<2>(v.get_price()))/2<<",";
//        std::cout << v.get_identifier() << ". "<< (get<1>(v.get_price())+get<2>(v.get_price()))/2<<" ";
        }
        std::cout << std::endl;

        myfile << market.trading_institutions.find(5)->second->dX(0,int((&watch)->current_time()))<<",";
        myfile << market.trading_institutions.find(5)->second->dX(1,int((&watch)->current_time()))<<",";
        myfile << market.trading_institutions.find(5)->second->dX(2,int((&watch)->current_time()))<<",";
        myfile << std::endl;
        day_count++;
    }
    market.clear(watch);
    std::cout << "market's watch time " << (&watch)->current_time() << std::endl;
}

//    (&watch)->tick();






myfile.close();
    return 0;
}
