#ifndef __NAIVESDP_H_INCLUDED__ //   if naiveSDP.h hasn't been included yet
#define __NAIVESDP_H_INCLUDED__ //   #define this so the compiler knows it has been included

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <utility>

#include "graph.h"
#include "interval.h"

using namespace std;

class Node;
class Interval;

class Match
{
public:
    int len;
    //len after possible 'cut' in SDP traceback
    int len_cut;
    pair<int, int> position;
    int cost;
    int intervalID;
    Match *preMatch;
    Match *nextMatch;
    Match ();
    Match(int _len, pair<int, int> _position);
};

//return min value of x and y
template <typename T>
T findMin(T x, T y);

//return max value of x and y
template <typename T>
T findMax(T x, T y);

//return match1 diag - match2 diag: the difference of diagnal
int diagDiff(Match &match1, Match &match2);

//extend intervals
//void extendIntervals(vector<Node> &graph, vector<Interval> &intervals, int k);

//get intervals in between junction nodes in topo order;
void getIntervals(vector<Node> &graph, vector<Interval> &intervals, int k);

//find matches on each interval;
void findMatchesHelp(int n, int m, int k, vector<Interval> &intervals, 
                        vector<Match> &matches, string &seq2, 
                        vector<int> &ifMatch, vector<Node> &graph);


//find matches of two sequences given k
void findMatches(string &seq2, vector<Match> &matches,
                 vector<int> &ifMatch, int n, int m, int k, 
                 vector<Node> &graph, vector<Interval> &intervals);

//output Matches SDP cost in terminal
void outputMatches(int n, int m, 
               vector<int> &ifMatch, vector<Match> &matches);

//trace back SDP results and add them to sdpMatches
void traceSDP(vector<Match> &matches, vector<Match> &sdpMatches,
                vector<int> &graphEnd, vector<int> &graphStart, vector<Node> &graph,
                vector<Interval> &intervals, int m, int n);

//merge overlapping matches
void mergeMatch(vector<Match> &matches, vector<int> &ifMatch, int n, int m, vector<Interval> &intervals);

//store interval ID to nodes;
void storeIntervalID(vector<Node> &graph, vector<Interval> &intervals);

//help function for getValidIntervals();
void getValidIntervalsHelp(int curIntervalID, vector<Interval> &intervals,
                            vector<int> &validIntervals, vector<Node> &graph);

//New: help function for getValidIntervals(): faster;
void NewgetValidIntervalsHelp(int curIntervalID, vector<Interval> &intervals,
                            vector<int> &validIntervals, vector<Node> &graph);

//get valid interval IDs for a given match;
vector<int> getValidIntervals(int curIntervalID, vector<Interval> &intervals, vector<Node> &graph);

//BFS on the graph and store depth in each node: for calculating distance between matches;
void bfsGraph(vector<Node> &graph, int pre, int now);

//get nodes with no preNode/no nextNode;
void getStartEnd(vector<Node> &graph, vector<int> &graphStart, vector<int> &graphEnd);

//run naiveSDP on given sequences, O(n^2) time;
void naiveSDP(string &seq2, vector<Match> &sdpMatches, int k, 
    vector<Node> &graph, int &cost);

//modified scoring function: for longest increasing subseq;
void traceSDP_1(vector<Match> &matches, vector<Match> &sdpMatches,
                vector<int> &graphEnd, vector<int> &graphStart, vector<Node> &graph, 
                vector<Interval> &intervals, int m, int n);

void naiveSDP_1(string &seq2, vector<Match> &sdpMatches, int k, vector<Node> &graph, int &cost);

#endif