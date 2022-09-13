#ifndef __INTERVAL_H_INCLUDED__ 
#define __INTERVAL_H_INCLUDED__

#include <iostream>
#include <vector>

using namespace std;

class Interval
{
public:
    //the length of interval;
    int len;
    //the start and end position on the graph;
    int start;
    //the end point is the starting position of the last node in interval;
    int end;
    //the index of matches in the interval
    vector<int> matchIndex;
    //index of parent and children intervals;
    ////vector<int> parentIntIndex;
    ////vector<int> childIntIndex;
    Interval(int _len, int _start, int _end);
};


#endif