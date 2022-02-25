#include <iostream>
#include "company.h"
#include "funds.h"
#include "price.h"
#include "order.h"
#include "price_setter.h"
#include "market_watch.h"
#include <fstream>


int main() {

//set overal wealth on the market

double market_wealth = 1.68e10;
double shares_outstanding = 0.8*market_wealth/300;
double market_cash = 0.2*market_wealth;

//Define the price setter
price_setter market;



//create assets and register stocks to the market
market.register_assets(3);
market.generate_fundamentals();


//generate initial quotes
market.set_initial_quotes();

//set price setter's wealth
market.initial_wealth(0, shares_outstanding);



//initialize traders
funds  aggressive_value
            , averse_value
                ,aggressive_noisy
                    , averse_noisy
                        ,aggressive_momentum
                            , averse_momentum
                                ,aggressive_growth
                                    , averse_growth
                ;


vector<funds*> participants = {&aggressive_value
                                    , &averse_value
                                        , &aggressive_noisy
                                            , &averse_noisy
                                                , &aggressive_momentum
                                                    , &averse_momentum
                                                        , &aggressive_growth
                                                             , &averse_growth
                                            };



//number of participants
int investor_population = participants.size();
//inform traders about assets on the market
    for(auto &p : participants){
        p->get_inform(market.tradeable_assets());
    }


    aggressive_value.trade_strategy(funds::aggressive_fundamental_value);
    aggressive_value.set_reb_period(63);
    aggressive_value.initial_wealth((1./investor_population)*market_cash,0.);

    averse_value.trade_strategy(funds::fundamental_value);
    averse_value.set_reb_period(63);
    averse_value.initial_wealth((1./investor_population)*market_cash,0.);

    aggressive_growth.trade_strategy(funds::aggressive_growth);
    aggressive_growth.set_reb_period(63);
    aggressive_growth.initial_wealth((1./investor_population)*market_cash,0.);

    averse_growth.trade_strategy(funds::growth);
    averse_growth.set_reb_period(63);
    averse_growth.initial_wealth((1./investor_population)*market_cash,0.);


    aggressive_noisy.trade_strategy(funds::aggressive_noise);
    aggressive_noisy.set_reb_period(1);
    aggressive_noisy.initial_wealth((1./investor_population)*market_cash,0.);

    averse_noisy.trade_strategy(funds::noise);
    averse_noisy.set_reb_period(1);
    averse_noisy.initial_wealth((1./investor_population)*market_cash,0.);


    aggressive_momentum.trade_strategy(funds::aggressive_momentum_investment);
    aggressive_momentum.set_reb_period(21);
    aggressive_momentum.initial_wealth((1./investor_population)*market_cash,0.);

    averse_momentum.trade_strategy(funds::momentum_investment);
    averse_momentum.set_reb_period(21);
    averse_momentum.initial_wealth((1./investor_population)*market_cash,0.);




//register traders on the market
market.register_traders(participants);


//market simulation point
std::ofstream myfile;

myfile.open ("NewMarketEcology.csv");


while(watch.current_time() < 10000){



int g = watch.current_time();

int kj = watch.current_time();
//receive orders from traders
    market.receive_orders();

//clear market



//std::cout<< watch.current_time() << ". ";

for(auto &[k,v] : market.tradeable_assets()){
        if(k==2){
            double st=0;
                myfile<<watch.current_time()<<","<<(get<1>(v.get_price())+get<2>(v.get_price()))/2<<","<<v.get_dividend(watch.current_time())<<","<<market.wealth<<",";
                for(auto &[x,y]:market.trading_institutions){
                    myfile << y->stocks_at_hand.find(2)->second<<",";
                     st += y->stocks_at_hand.find(2)->second;
                    std::cout <<y->get_identifier() <<" "<<y->stocks_at_hand.find(2)->second<<",";

                }
            myfile<< market.stocks_at_hand.find(2)->second<<",";
                myfile<< st + market.stocks_at_hand.find(2)->second<<",";
            }

        //myfile<<(get<1>(v.get_price())+get<2>(v.get_price()))/2<<",";
//        std::cout << v.get_identifier() << ". "<< (get<1>(v.get_price())+get<2>(v.get_price()))/2<<" ";
}
//    std::cout<<std::endl;
    myfile << std::endl;
    market.clear();
    watch.tick();

    //update traders balances before orders
    for(auto &p : participants){
        p->balance_cf(watch.current_time());
    }

}

myfile.close();

    return 0;
}
