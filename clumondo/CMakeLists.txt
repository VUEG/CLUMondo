cmake_minimum_required(VERSION 3.6)

set(SOURCE_FILES main.cpp clumondomodules.cpp globals.cpp util.cpp
        ../include/CLUMondoConfig.h ../include/clumondomodules.h ../include/globals.h ../include/util.h)
add_executable(clumondo ${SOURCE_FILES})
