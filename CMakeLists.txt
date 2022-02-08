cmake_minimum_required(VERSION 3.17)
project(untitled32)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)
find_package (Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS ALL)
find_package(PythonLibs)

add_executable(untitled32 main.cpp company.cpp company.h currencies.h market_watch.h market_watch.cpp price.h price.cpp identity.h identity.cpp value_investor.h value_investor.cpp fund.cpp fund.h momentum_investor.cpp momentum_investor.h noisy_investor.cpp noisy_investor.h funds.cpp funds.h order.cpp order.h stochastic_math.cpp stochastic_math.h price_setter.cpp price_setter.h wiener.h)
target_link_libraries(untitled32 ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
target_link_libraries(untitled32 Eigen3::Eigen)