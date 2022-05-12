//
// Created by Cephas Svosve on 20/12/2021.
//
#include "market_watch.h"
#include <iostream>

void market_watch::start() {
    start_time = std::chrono::high_resolution_clock::now();
}


double
market_watch::current_time(){
    auto stop = std::chrono::high_resolution_clock::now();

    switch (day_length) {
        case 0:
            return duration_cast<std::chrono::microseconds>(stop - start_time).count();
        case 1:
            return duration_cast<std::chrono::milliseconds>(stop - start_time).count();
        case 2:
            return duration_cast<std::chrono::seconds>(stop - start_time).count();
        case 3:
            return duration_cast<std::chrono::minutes>(stop - start_time).count();

            // return (end_time-start_time)/std::chrono::miiseconds(1);


    }

}

void
market_watch::reset(){
    start_time = std::chrono::high_resolution_clock::now();
}

time_point<std::chrono::steady_clock, std::chrono::duration<long long int, std::ratio<1, 1000000000>>>
market_watch::get_creation_time() {
    return creation_time;
}
