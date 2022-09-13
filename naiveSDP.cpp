#include <iostream>
#include <string>
#include <assert.h>
#include <math.h>
#include <vector>
#include <array>
#include <utility>
#include <stack>
#include <list>

#include "naiveSDP.h"
#include "graph.h"
#include "interval.h"

using namespace std;

//constructor
Match::Match() {};
Match::Match(int _len, pair<int, int> _position) : len(_len),
                                                position(_position) {}

//return min value of x and y
template <typename T>
T findMin(T x, T y) {
    return x < y ? x : y;
}

//return max value of x and y
template <typename T>
T findMax(T x, T y) {
    return x > y ? x : y;
}

//return match1 diag - match2 diag: the difference of diagnal;
int diagDiff(Match &match1, Match &match2) {
    return -(match1.position.first - match2.position.first) +
           (match1.position.second - match2.position.second);
}

//merge overlapping matches
void mergeMatch(vector<Match> &matches, vector<int> &ifMatch, int n, int m, vector<Interval> &intervals) {
    //store the number of merges so far
    //int no_merge = 0;
    vector<Interval>::iterator it_int = intervals.begin();
    for (; it_int != intervals.end(); ++it_int) {
        if ((*it_int).matchIndex.size() > 0) {
            int match_start = (*it_int).matchIndex[0];
            int match_end = (*it_int).matchIndex[(*it_int).matchIndex.size() - 1];
            //naive: for each match in current interval, loop throght all previous matches
            for (int i = match_start ; i <= match_end; ++i) {
                //test
                cout << "have " << matches[i].position.first << " "  << matches[i].position.second << 
                        " len " << matches[i].len << 
                        " interval no " << it_int - intervals.begin() + 1 << endl << endl;
                for (int j = match_start; j < i; ++j){
                    if (diagDiff(matches[i], matches[j]) == 0 && 
                        (matches[i].position.first - matches[j].position.first <= matches[j].len)) {
                        matches[j].len = matches[i].len + (matches[i].position.first - matches[j].position.first);
                        ifMatch[matches[i].position.first * m + matches[i].position.second] = 0;
                        //test
                        cout << "delete " << matches[i].position.first << " " << matches[i].position.second << endl << endl;
                        matches.erase(matches.begin() + i);
                        (*it_int).matchIndex.erase((*it_int).matchIndex.begin() + (*it_int).matchIndex.size() - 1);
                        i = i - 1;
                        match_end = match_end - 1;
                        break;                
                    }
                }
            }
        }
    }
}

//get intervals in between junction nodes in topo order;
void getIntervals(vector<Node> &graph, vector<Interval> &intervals, int k) {
    int intStart, intEnd;
    stack<int> endPoints;
    ////stack<vector<int>> preIntervals;
    intStart = 0;
    vector<Node>::iterator it_graph = graph.begin();
    for(; it_graph != graph.end(); ++it_graph) {
        //conditions can be simplified;
        if ((it_graph->nextNodeVec).size() != 1) {
            intEnd = it_graph - graph.begin();
            intervals.push_back(Interval(intEnd - intStart + 1, intStart, intEnd));
            intStart = intEnd + 1;
        }
        else if ((it_graph->nextNodeVec)[0] != it_graph - graph.begin() + 1) {
            if (endPoints.size() != 0 && (it_graph->nextNodeVec)[0] == endPoints.top()) {
                intEnd = it_graph - graph.begin();
                intervals.push_back(Interval(intEnd - intStart + 1, intStart, intEnd));
                intStart = intEnd + 1;
            }
            else {
                endPoints.push((it_graph->nextNodeVec)[0]);
                intEnd = it_graph - graph.begin();
                intervals.push_back(Interval(intEnd - intStart + 1, intStart, intEnd));
                intStart = intEnd + 1;
            }
        }
        else {
            if (endPoints.size() != 0 && (it_graph->nextNodeVec)[0] == endPoints.top()) {
                endPoints.pop();
                intEnd = it_graph - graph.begin();
                intervals.push_back(Interval(intEnd - intStart + 1, intStart, intEnd));
                intStart = intEnd + 1;
            }
        }
    }
}

//find matches on each interval;
void findMatchesHelp(int n, int m, int k, vector<Interval> &intervals, 
                        vector<Match> &matches, string &seq2, 
                        vector<int> &ifMatch, vector<Node> &graph) {
    //TODO: Search from nodes only, otherwise notes that each end of an interval contribute k bases to the seq1;
    //int offsetBase = 0;
    vector<Interval>::iterator it_int = intervals.begin();
    int start, end;
    for (; it_int != intervals.end(); ++it_int) {
        start = it_int->start;
        end = it_int->end;
        //note that end is an index, m is length;
        //start and end could be the same;
        for (int i = start; i <= end; ++i) {
            for (int j = 0; j < m - k + 1; ++j) {
                for (int l = 0; l < k; ++l) {
                    if ((graph[i].base)[l] != seq2[j + l])
                        break;
                    if (l == k - 1) {
                        ////ifMatch[i * m + j] = 1;
                        Match curMatch(k, make_pair(i, j));
                        //intervalID starting from 0;
                        curMatch.intervalID = it_int - intervals.begin();
                        matches.push_back(curMatch);
                        //store match's index to current interval;
                        (*it_int).matchIndex.push_back(matches.size() - 1);
                    }
                }
            }
        }
    }
}

//find matches of two sequences given k;
void findMatches(string &seq2, vector<Match> &matches,
                 vector<int> &ifMatch, int n, int m, int k, 
                 vector<Node> &graph, vector<Interval> &intervals) {
    //Question: if we search by 'intervals', matches span two intervals will
    //not be considered.??? (Set k small???)
    //merge for each interval

    //get interval of each branch in the order of branch's topo order;
    intervals.clear();
    getIntervals(graph, intervals, k);

    //extend intervals
    //extendIntervals(graph, intervals, k);

    //TODO: double check the results, and handle the boundry problem (question above);
    findMatchesHelp(n, m, k, intervals, matches, seq2, ifMatch, graph);

    //merge matches in each interval;
    //mergeMatch(matches, ifMatch, n, m, intervals);
}

//output Matches SDP cost in terminal
void outputMatches(int n, int m, 
               vector<int> &ifMatch, vector<Match> &matches) {
    int match_count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (ifMatch[i * m + j] == 1) {
                cout << matches[match_count].cost;
                //cout << ifMatch[i * m + j];
                match_count += 1;
            }
            else
                cout << "_";
                //cout << ifMatch[i * m + j];
        }
        cout << endl;
    }
}

//trace back SDP results and add them to sdpMatches
//TODO: modify: too complicated, use a seperate function for scoring;
void traceSDP(vector<Match> &matches, vector<Match> &sdpMatches,
                vector<int> &graphEnd, vector<int> &graphStart, vector<Node> &graph, 
                vector<Interval> &intervals, int m, int n) {
    //calculate score for each end point;
    int bestEndPointPreMatch, curEndPointPreMatch;
    vector<int> curValidIntervals;
    //TODO: change with scoring function accordingly;
    int sdpCost = m + n;
    //TODO: change scoring function;
    for (vector<int>::iterator it_graphEnd = graphEnd.begin(); it_graphEnd != graphEnd.end(); ++it_graphEnd) {
        int curEnd = *it_graphEnd;
        int curIntervalID = graph[curEnd].intervalID;
        curValidIntervals = getValidIntervals(curIntervalID, intervals, graph);
        //boolean: current end point has valid match or not;
        bool ifValidMatch = false;
        //loop through valid intervals;
        vector<int>::iterator it_validInt = curValidIntervals.begin();
        for (; it_validInt != curValidIntervals.end(); ++it_validInt) {
            int curValidInt = *it_validInt;
            //if there is a match in curValidInt, set it (the last match) to curEndPointPreMatch;
            int size = intervals[curValidInt].matchIndex.size();
            if (size != 0) {
                curEndPointPreMatch = intervals[curValidInt].matchIndex[size - 1];
                ifValidMatch = true;
                break;
            }
        }
        //use the same scoring here as in the naiveSDP();
        if (ifValidMatch) {
            bfsGraph(graph, matches[curEndPointPreMatch].position.first, curEnd);
            //x axis: graph axis, y axis: seq axis;
            //the x axis distance;
            int xDis = findMax<int>(0, graph[curEnd].depth
                                        - graph[matches[curEndPointPreMatch].position.first].depth 
                                        - matches[curEndPointPreMatch].len + 1);
            //precedes;
            if (xDis > 0 && matches[curEndPointPreMatch].position.second + matches[curEndPointPreMatch].len <= m - 1) {
                //the y axis distance;
                int yDis = m - 1 - matches[curEndPointPreMatch].position.second - matches[curEndPointPreMatch].len + 1;
                if (matches[curEndPointPreMatch].cost + xDis + yDis < sdpCost) {
                    sdpCost = matches[curEndPointPreMatch].cost + xDis + yDis;
                    bestEndPointPreMatch = curEndPointPreMatch;
                }
            }
            //visl or visa;
            else {
                //the y axis distance;
                int yDis = findMax<int>(0, m - 1  
                                            - matches[curEndPointPreMatch].position.second - matches[curEndPointPreMatch].len + 1);
                if (matches[curEndPointPreMatch].cost + xDis + yDis < sdpCost) {
                    sdpCost = matches[curEndPointPreMatch].cost + xDis + yDis;
                    bestEndPointPreMatch = curEndPointPreMatch;
                }        
            }
        }
        else {
            //use the distance from the graphStarts;
            vector<int>::iterator it_graphStart = graphStart.begin();
            for (; it_graphStart != graphStart.end(); ++it_graphStart) {
                bfsGraph(graph, *it_graphStart, curEnd);
                int nodeDistance = graph[curEnd].depth - graph[*it_graphStart].depth;
                int curRawCost = nodeDistance + m - 1;
                if (curRawCost < sdpCost) {
                    sdpCost = curRawCost;
                    //no best prematch: not inherit scores from previous matches;
                    bestEndPointPreMatch = -1;
                }
            }
        }
        //matches.push_back(Match(1, make_pair(n - 1, m - 1)));
    }

    //add the last match;
    if (bestEndPointPreMatch != -1) {
        Match curMatch = matches[bestEndPointPreMatch];
        curMatch.len_cut = curMatch.len;
        sdpMatches.push_back(curMatch);
        Match oldMatch = curMatch;
        while (curMatch.preMatch != NULL) {
            curMatch = *(curMatch.preMatch);
            //when visl or visa happens, store part of match to sdpMatches
            //also, make sure there is no overlapping matches
            //TODO: the match before the fake match, if contact to edges, will be cut unproperly;
            if (curMatch.position.first + curMatch.len > oldMatch.position.first &&
                    diagDiff(oldMatch, curMatch) > 0)
                curMatch.len_cut = oldMatch.position.first - curMatch.position.first;
            else if (curMatch.position.second + curMatch.len > oldMatch.position.second &&
                    diagDiff(oldMatch, curMatch) < 0)
                curMatch.len_cut = oldMatch.position.second - curMatch.position.second;
            else 
                curMatch.len_cut = curMatch.len;
            sdpMatches.push_back(curMatch);
            oldMatch = curMatch;
        }
    }
    //no match selected;
    //else {
        //
    //}
}

//store interval ID to nodes;
//O(n);
//TODO: then do not need store interval info in matches;
void storeIntervalID(vector<Node> &graph, vector<Interval> &intervals) {
    vector<Interval>::iterator it_int = intervals.begin();
    int curStart, curEnd;
    for (; it_int != intervals.end(); ++it_int) {
        curStart = (*it_int).start;
        curEnd = (*it_int).end;
        for (int i = curStart; i <= curEnd; ++i) {
            graph[i].intervalID = it_int - intervals.begin();
        }
    }
}

//help function for getValidIntervals();
void getValidIntervalsHelp(int curIntervalID, vector<Interval> &intervals,
                            vector<int> &validIntervals, vector<Node> &graph) {
    Node curBeginNode = graph[intervals[curIntervalID].start];
    //the idea is, every preNode of curBeginNode comes from a valid interval, store them;
    //call the function recursively;
    int preIntervalID;
    if (curBeginNode.preNodeVec.size() != 0) {
        vector<int>::iterator it_pre = curBeginNode.preNodeVec.begin();
        for (; it_pre != curBeginNode.preNodeVec.end(); ++it_pre) {
            preIntervalID = graph[*it_pre].intervalID;
            if (!checkExist(validIntervals, preIntervalID))
                validIntervals.push_back(preIntervalID);
            getValidIntervalsHelp(preIntervalID, intervals, validIntervals, graph);
        }
    }
}

//New: help function for getValidIntervals(): faster;
//TODO!!!: test NewgetValidIntervalsHelp;
void NewgetValidIntervalsHelp(int curIntervalID, vector<Interval> &intervals,
                            vector<int> &validIntervals, vector<Node> &graph) {
    //store visit information;
    vector<bool> visit(intervals.size(), false);
    //use list as a queue: store node positions in graph;
    list<int> que;
    //initialization
    que.push_back(curIntervalID);
    visit[curIntervalID] = true;
    int preIntervalID;
    //the beginning node of an interval;
    Node curBeginNode;
    //do bfs
    while (!que.empty()) {
        curIntervalID = que.front();
        que.pop_front();
        //get the beginning node of the current interval;
        curBeginNode = graph[intervals[curIntervalID].start];
        vector<int>::iterator it_pre = curBeginNode.preNodeVec.begin();
        for (; it_pre != curBeginNode.preNodeVec.end(); ++it_pre) {
            if (!visit[graph[*it_pre].intervalID]) {
                que.push_back(graph[*it_pre].intervalID);
                visit[graph[*it_pre].intervalID] = true;
                preIntervalID = graph[*it_pre].intervalID;
                validIntervals.push_back(preIntervalID);
            }
        }
    }
}

//get valid interval IDs for a given match;
vector<int> getValidIntervals(int curIntervalID, vector<Interval> &intervals, vector<Node> &graph) {
    vector<int> validIntervals;
    //intervalID: the position in intervals;
    //int curIntervalID = match.intervalID;
    //add current interval;
    validIntervals.push_back(curIntervalID);
    //add intervals by preNode of the first nodes of valid intervals;
    NewgetValidIntervalsHelp(curIntervalID, intervals, validIntervals, graph);
    return validIntervals;
}

//BFS on two nodes: for calculating distance between matches;
//TODO: too slow, think about it;
void bfsGraph(vector<Node> &graph, int pre, int now) {
    //use list as a queue: store node positions in graph;
    list<int> que;
    //store visit information;
    vector<bool> visit(graph.size(), false);
    //initialization
    que.push_back(pre);
    visit[pre] = true;
    //depth starts from 0;
    graph[pre].depth = 0;
    //do BFS and store depth
    int cur;
    bool find = false;
    while (!que.empty()) {
        cur = que.front();
        que.pop_front();
        vector<int>::iterator it = graph[cur].nextNodeVec.begin();
        for (; it != graph[cur].nextNodeVec.end(); ++it) {
            if (!visit[*it]) {
                if (*it == now) {
                    graph[*it].depth = graph[cur].depth + 1;
                    find = true;
                    break;
                }
                que.push_back(*it);
                visit[*it] = true;
                graph[*it].depth = graph[cur].depth + 1;
            }
        }
        if (find) break;
    }
    //if there is no tree struc between pre and now;
    //this is for matches on starting intervals;
    if (!find) graph[now].depth = graph.size();
}

//get nodes with no preNode/no nextNode;
//TODO: still bugs;
void getStartEnd(vector<Node> &graph, vector<int> &graphStart, vector<int> &graphEnd) {
    vector<Node>::iterator it_graph = graph.begin();
    for (; it_graph != graph.end(); ++it_graph) {
        if ((*it_graph).preNodeVec.size() == 0)
            graphStart.push_back(it_graph - graph.begin());
        if ((*it_graph).nextNodeVec.size() == 0)
            graphEnd.push_back(it_graph - graph.begin());
    }
}


void naiveSDP(string &seq2, vector<Match> &sdpMatches, int k, vector<Node> &graph, int &cost) {
    //TODO: learn implementations of DP/SDP;
    //TODO: modify scoring function;

    //seq1: the graph seq and the x axis;
    //seq2: the new seq and the y axis;

    //clear the previous sdp results
    sdpMatches.clear();
    //store all the matches;
    vector<Match> matches;
    //int n = seq1.length();
    int n = graph.size();
    int m = seq2.length();
    //store intervals
    vector<Interval> intervals;
    vector<int> ifMatch(n * m, 0);

    //matches are in (x, y) order;
    findMatches(seq2, matches, ifMatch, n, m, k, graph, intervals);

    //store interval ID to each node;
    storeIntervalID(graph, intervals);

    //get start and end nodes of current graph;
    vector<int> graphStart;
    vector<int> graphEnd;
    getStartEnd(graph, graphStart, graphEnd);

    if (matches.size() != 0) {
        //add fake matches at graphEnd positions: easier to trace back;
        ////ifMatch[(n - 1) * m + (m - 1)] = 1;

        cout << "No. of matches: " << matches.size() << endl;
        //store valid intervals
        vector<int> validIntervals;
        //loop through all matches
        //TODO: find the correct cost for the last and first matches with BFS depth;
        for (int i = 0; i < matches.size(); ++i) {
            //get valid intervals for current match;
            validIntervals.clear();
            validIntervals = getValidIntervals(matches[i].intervalID, intervals, graph);

            //TODO: test more about getValidIntervals();
            //test
            //for (int p = 0; p < validIntervals.size(); p++)
            //    cout << i << " intervals:  " << validIntervals[p] << endl;

            //starting from the current interval;
            int curMinCost = m + n;
            //Raw cost: not using any previous matches;
            int curRawCost;
            //distance between two nodes on the bfs tree;
            int nodeDistance;
            //set curMinCost to be the smallest curRawCost;
            vector<int>::iterator it_graphStart = graphStart.begin();
            for (; it_graphStart != graphStart.end(); ++it_graphStart) {
                //TODO: improve: the way bfs is used here is super slow;
                bfsGraph(graph, *it_graphStart, matches[i].position.first);
                nodeDistance = graph[matches[i].position.first].depth 
                                - graph[*it_graphStart].depth;
                curRawCost = nodeDistance + matches[i].position.second;
                curMinCost = findMin<int>(curMinCost, curRawCost);
            }

            //pointer of the selected previous match;
            Match *preMatch = NULL;

            //TODO: debug;
            //loop through all the valid intervals;
            vector<int>::iterator it_validInt = validIntervals.begin();
            for (; it_validInt != validIntervals.end(); ++it_validInt) {
                //loop through matches;
                vector<int>::iterator it_matchInd = intervals[*it_validInt].matchIndex.begin();
                for (; it_matchInd != intervals[*it_validInt].matchIndex.end(); ++it_matchInd) {
                    //if valid matches: both x and y coordinates;
                    if (matches[*it_matchInd].position.first < matches[i].position.first && 
                        matches[*it_matchInd].position.second < matches[i].position.second) {
                        bfsGraph(graph, matches[*it_matchInd].position.first, matches[i].position.first);
                        //x axis: graph axis, y axis: seq axis;
                        //the x axis distance;
                        int xDis = findMax<int>(0, graph[matches[i].position.first].depth
                                                    - graph[matches[*it_matchInd].position.first].depth 
                                                    - matches[*it_matchInd].len + 1);
                        //precedes;
                        if (xDis > 0 && matches[*it_matchInd].position.second + matches[*it_matchInd].len <= matches[i].position.second) {
                            //the y axis distance;
                            int yDis = matches[i].position.second - matches[*it_matchInd].position.second - matches[*it_matchInd].len + 1;
                            if (matches[*it_matchInd].cost + xDis + yDis < curMinCost) {
                                curMinCost = matches[*it_matchInd].cost + xDis + yDis;
                                preMatch = &matches[*it_matchInd];
                            }
                        }
                        //visl or visa;
                        else {
                            //the y axis distance;
                            int yDis = findMax<int>(0, matches[i].position.second 
                                                        - matches[*it_matchInd].position.second - matches[*it_matchInd].len + 1);
                            if (matches[*it_matchInd].cost + xDis + yDis < curMinCost) {
                                curMinCost = matches[*it_matchInd].cost + xDis + yDis;
                                preMatch = &matches[*it_matchInd];
                            }
                            
                        }
                    }
                }
            }

            /*
            int curMinCost = m + n;
            //Raw cost: not using any previous matches;
            //TODO: need to be modify in POA (BFS);
            int curRawCost = matches[i].position.first + matches[i].position.second;
            curMinCost = findMin<int>(curMinCost, curRawCost);
            //pointer of the selected previous match
            Match *preMatch = NULL;

            //check if match j is valid for match i: at least share one sequence;
            bool valid = false;
            //double check +1/-1 in cost/distance calculations
            for (int j = 0; j < i; ++j) {
                valid = false;
                //for the artificial match, all previous matches are valid;
                if (i == matches.size() - 1)
                    valid = true;
                else {
                    for (vector<int>::iterator it_id = graph[matches[i].position.first].seq_ID.begin();
                                            it_id != graph[matches[i].position.first].seq_ID.end();
                                            ++it_id) {
                        if (checkExist(graph[matches[j].position.first].seq_ID, *it_id)){
                            valid = true;
                            break;
                        }
                    }
                }
                
                if (valid) {
                    if (matches[j].position.first < matches[i].position.first &&
                        matches[j].position.second < matches[i].position.second) {
                        //precedes
                        if (matches[j].position.first + matches[j].len <= matches[i].position.first &&
                            matches[j].position.second + matches[j].len <= matches[i].position.second) {
                            if (matches[j].cost +
                                    (matches[i].position.first - (matches[j].position.first + matches[j].len - 1) +
                                    (matches[i].position.second - (matches[j].position.second + matches[j].len - 1))) <
                                curMinCost) {
                                curMinCost = matches[j].cost +
                                            (matches[i].position.first - (matches[j].position.first + matches[j].len - 1) +
                                            (matches[i].position.second - (matches[j].position.second + matches[j].len - 1)));
                                preMatch = &matches[j];
                            }
                        }
                        //visl or visa
                        else {
                            if (matches[j].cost + abs(diagDiff(matches[i], matches[j])) < curMinCost) {
                                curMinCost = matches[j].cost + abs(diagDiff(matches[i], matches[j]));
                                preMatch = &matches[j];
                            }
                        }
                    }
                }
            }
            */
            matches[i].cost = curMinCost;
            matches[i].preMatch = preMatch;
        }
    }

    if (matches.size() != 0) {
        traceSDP(matches, sdpMatches, graphEnd, graphStart, graph, intervals, m, n);
        cost = sdpMatches[0].cost;
        outputMatches(n, m, ifMatch, matches);
    }

    //delete the last fake match;
    //note it's added after merging;
    //sdpMatches.erase(sdpMatches.begin());

    //test
    cout << "naive sdpMatches: " << endl;
    for (int i = 0; i < sdpMatches.size(); ++i) {
        cout << sdpMatches[i].position.first << " " << sdpMatches[i].position.second 
            << " cost " << sdpMatches[i].cost 
            << " len " <<  sdpMatches[i].len 
            << " len_cut " << sdpMatches[i].len_cut << endl << endl;
    }
}


void traceSDP_1(vector<Match> &matches, vector<Match> &sdpMatches,
                vector<int> &graphEnd, vector<int> &graphStart, vector<Node> &graph, 
                vector<Interval> &intervals, int m, int n) {
    int bestMatch, bestCost;
    bestCost = 0;
    vector<Match>::iterator it_matches = matches.begin();
    for (; it_matches != matches.end(); ++it_matches) {
        if ((*it_matches).cost > bestCost) {
            bestCost = (*it_matches).cost;
            bestMatch = it_matches - matches.begin();
        }
    }

    Match curMatch = matches[bestMatch];
    curMatch.len_cut = curMatch.len;
    //curMatch.len_cut = curMatch.len;
    sdpMatches.push_back(curMatch);
    Match oldMatch = curMatch;
    while (curMatch.preMatch != NULL) {
        curMatch = *(curMatch.preMatch);
        curMatch.len_cut = curMatch.len;
        sdpMatches.push_back(curMatch);
        oldMatch = curMatch;
    }
}

//longest increasing subseq: count the most matches;
void naiveSDP_1(string &seq2, vector<Match> &sdpMatches, int k, vector<Node> &graph, int &cost) {
    //TODO: learn implementations of DP/SDP;
    //TODO: modify scoring function;

    //seq1: the graph seq and the x axis;
    //seq2: the new seq and the y axis;

    //clear the previous sdp results
    sdpMatches.clear();
    //store all the matches;
    vector<Match> matches;
    //int n = seq1.length();
    int n = graph.size();
    int m = seq2.length();
    //store intervals
    vector<Interval> intervals;
    //TODO: change all the ifMatch: in two cases;
    vector<int> ifMatch(1, 0);

    //matches are in (x, y) order;
    //do not merge in this case;
    findMatches(seq2, matches, ifMatch, n, m, k, graph, intervals);

    //store interval ID to each node;
    storeIntervalID(graph, intervals);

    //get start and end nodes of current graph;
    vector<int> graphStart;
    vector<int> graphEnd;
    getStartEnd(graph, graphStart, graphEnd);

    //scoring matches;
    if (matches.size() != 0) {
        cout << "No. of matches: " << matches.size() << endl;
        //store valid intervals
        vector<int> validIntervals;
        //loop through all matches
        for (int i = 0; i < matches.size(); ++i) {
            //get valid intervals for current match;
            //validIntervals.clear();
            validIntervals = getValidIntervals(matches[i].intervalID, intervals, graph);

            //TODO: test more about getValidIntervals();
            //test
            //cout << i << endl;
            //for (int p = 0; p < validIntervals.size(); p++)
            //    cout << i << " intervals:  " << validIntervals[p] << endl;

            //starting from the current interval;
            int curMinCost = m + n;
            //Raw cost: not using any previous matches;
            //here cost is the number of matches so far;
            //here curRawCost = 1 because its the first match;
            int curRawCost  = 1;
            curMinCost = findMin<int>(curMinCost, curRawCost);

            //pointer of the selected previous match;
            Match *preMatch = NULL;

            //loop through all the valid intervals;
            vector<int>::iterator it_validInt = validIntervals.begin();
            for (; it_validInt != validIntervals.end(); ++it_validInt) {
                //loop through matches;
                vector<int>::iterator it_matchInd = intervals[*it_validInt].matchIndex.begin();
                for (; it_matchInd != intervals[*it_validInt].matchIndex.end(); ++it_matchInd) {
                    //if valid matches: both x and y coordinates;
                    if (matches[*it_matchInd].position.first < matches[i].position.first && 
                        matches[*it_matchInd].position.second < matches[i].position.second) {
                        if (matches[*it_matchInd].cost >= curMinCost) {
                            curMinCost = matches[*it_matchInd].cost + 1;
                            preMatch = &matches[*it_matchInd];
                        }
                    }
                }
            }
            matches[i].cost = curMinCost;
            matches[i].preMatch = preMatch;
        }
    }

    if (matches.size() != 0) {
        traceSDP_1(matches, sdpMatches, graphEnd, graphStart, graph, intervals, m, n);
        cost = sdpMatches[0].cost;
        //outputMatches(n, m, ifMatch, matches);
    }

    //test
    /*
    cout << "naive sdpMatches: " << endl;
    for (int i = 0; i < sdpMatches.size(); ++i) {
        cout << sdpMatches[i].position.first << " " << sdpMatches[i].position.second 
            << " cost " << sdpMatches[i].cost 
            << " len " <<  sdpMatches[i].len 
            << " len_cut " << sdpMatches[i].len_cut << endl << endl;
    }
    */
    cout << "naive sdpMatches score: " << cost << endl;
}