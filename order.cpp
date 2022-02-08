//
// Created by Cephas Svosve on 4/1/2022.
//

#include "order.h"

//getters
int order::get_id(){ return this->order_id;}
bool order::get_status(){return this->order_status;}
int order::get_ordered_asset(){return this->asset_id;}
bool order::get_order_type(){return this->order_type;}
double order::get_order_size(){return this->order_size;}
double order::get_proposed_price(){return this->proposed_price;}


//setters
void order::set_id(int id){
    this->order_id = id;
}

void order::set_status(ORDER_STATUS st){
    this->order_status = st;
}

void order::set_ordered_asset(int id){
    this->asset_id = id;
}

void order::set_order_type(ORDER_TYPE st){
    this->order_type = st;
}

void order::set_order_size(double number_of_stocks, double proposed_price_){
    this->order_size = number_of_stocks;
    this->proposed_price = proposed_price_;
}

//define addition of two orders
order order::operator+ (const order& p){
        this->order_size = this->order_size + p.order_size;
        return *this;
}


order order::operator+= (const order& p){
        this->order_size += p.order_size;
        return *this;
}
