//
// Created by Cephas Svosve on 20/12/2021.
//

#ifndef UNTITLED30_TIME_POINT_H
#define UNTITLED30_TIME_POINT_H
class market_watch{
public:
    int lower=0;
    int upper;

    void tick();
    int current_time();
};

static market_watch watch;
#endif //UNTITLED30_TIME_POINT_H
