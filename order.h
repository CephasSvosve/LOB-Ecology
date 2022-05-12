//
// Created by Cephas Svosve on 4/1/2022.
//

#ifndef UNTITLED32_ORDER_H
#define UNTITLED32_ORDER_H
#include "price.h"

class order {
private:

    int order_id=0;//customer identity
    bool order_status=false;
    int asset_id=0;
    bool order_type=false;// can be limit or market
    double order_size=0;//number of stocks required
    double proposed_price=0;



public:
order(){};

enum ORDER_TYPE{limit = 1, market = 0};
enum ORDER_STATUS{active = 1, filled = 0};

//getters
int get_id();
bool get_status();
int get_ordered_asset();
bool get_order_type();
double get_order_size();
double get_proposed_price();



//setters
void set_id(int id);
void set_status(ORDER_STATUS);
void set_ordered_asset(int id);
void set_order_type(ORDER_TYPE);
void set_order_size(double number_of_stocks, double limit_price_);

//define post-addition of two prices
order operator+ (const order& p);

order operator+= (const order& p);

};


#endif //UNTITLED32_ORDER_H
