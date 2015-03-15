#!/bin/bash

# PATH TO HEADER FILES
CXXFLAGS="-I ../lib/netevo/ -I ../lib/include -I /Users/Tom/Development/Library/include -I ."

g++ $CXXFLAGS ../lib/netevo/gml.cc ../lib/netevo/simulate.cc ../lib/netevo/system.cc dynamic_nets.cc -o dynNet -O3 
