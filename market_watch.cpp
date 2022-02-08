//
// Created by Cephas Svosve on 20/12/2021.
//
#include "market_watch.h"
#include <iostream>

void market_watch::tick(){
    this->lower++;
    this->upper++;
}

int market_watch::current_time(){
    return this->lower;
}