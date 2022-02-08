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

double market_wealth = 10'000;
//Define the price setter
price_setter market;



//create assets and register stocks to the market
market.register_assets(3);
market.generate_fundamentals();


//generate initial quotes
market.set_initial_quotes();

//set price setter's wealth
market.initial_wealth(market_wealth, 5'000);

//initialize traders
funds value_investor,  noisy_investors, momentum_investors;
vector<funds*> participants = {&value_investor
                                    , &noisy_investors
                                        , &momentum_investors
                                            };




//inform traders about assets on the market
    for(auto &p : participants){
        p->get_inform(market.tradeable_assets());
    }


    value_investor.trade_strategy(funds::fundamental_value);
    value_investor.set_reb_period(1);
    value_investor.initial_wealth(0.25*market_wealth,500);

//    growth_investor.trade_strategy(funds::growth);
//    growth_investor.set_reb_period(63);
//    growth_investor.initial_wealth(0.25*market_wealth,500'000);

    noisy_investors.trade_strategy(funds::noise);
    noisy_investors.set_reb_period(63);
    noisy_investors.initial_wealth(0.25*market_wealth,500);

    momentum_investors.trade_strategy(funds::momentum_investment);
    momentum_investors.set_reb_period(200);
    momentum_investors.initial_wealth(0.25*market_wealth,500);




//register traders on the market
market.register_traders(participants);


//market simulation point
std::ofstream myfile;

myfile.open ("NewMarketEcology.csv");


do{


    //update traders balances before orders
    for(auto &p : participants){
        p->balance_cf(watch.current_time());


}


     //receive orders from traders
    market.receive_orders();

    //clear market
   market.clear();
//   std::cout<< watch.current_time() << ". ";

    for(auto &[k,v] : market.tradeable_assets()){
        if(k==2){
            if(!market.bid.find(k)->second.empty()){
            myfile<<watch.current_time()<<","<<v.get_earnings(watch.current_time())<<","<<market.bid.find(k)->second.rbegin()->get_id()<<","<<market.bid.find(k)->second.rbegin()->get_proposed_price();
        }else{
                myfile<<watch.current_time()<<","<<v.get_earnings(watch.current_time())<<",";
            } }
        //myfile<<(get<1>(v.get_price())+get<2>(v.get_price()))/2<<",";
//        std::cout << v.get_identifier() << ". "<< (get<1>(v.get_price())+get<2>(v.get_price()))/2<<" ";
    }
//    std::cout<<std::endl;
    myfile << std::endl;

    watch.tick();

}while(watch.current_time() < 1000);

myfile.close();

    return 0;
}
