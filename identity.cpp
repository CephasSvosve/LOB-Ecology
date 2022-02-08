//
// Created by Cephas Svosve on 20/12/2021.
//
#include "identity.h"
#include <iostream>

int create_identity(identity_type a){
    int i;
    try {
        switch(a){

            case 0:
                i = companies.size();
                companies.emplace(i, company_);
                break;
            case 1:
                i = assets.size();
                assets.emplace(i,asset);
                break;
            case 2:
                i = traders.size();
                traders.emplace(i,trader);
                break;
            case 3:
                i = fund_traders.size();
                fund_traders.emplace(i,fund_trader);
                break;
            case 4:
                i = pricers.size();
                pricers.emplace(i,pricer);
                break;
            case 5:
                i = trade_strategies.size();
                trade_strategies.emplace(i, strategy);
                break;
            default:
                throw "failed to create identity";

        }
        return i;
    }catch (...) {
        std::cout<<"failed to create identity"<<std::endl;
    }

}
