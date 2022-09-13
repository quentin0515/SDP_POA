#ifndef __GRAPH_H_INCLUDED__   // if naiveSDP.h hasn't been included yet
#define __GRAPH_H_INCLUDED__   //   #define this so the compiler knows it has been included

#include <string>
#include <vector>
//#include <stack>

#include "naiveSDP.h"

using namespace std;

class Match;

class Node {
public:
    string base;
    //store positions of the previous nodes in the topo sorted vector: can have multiple;
    vector <int> preNodeVec;
    //store positions of the next nodes in the topo sorted vector: can have multiple;
    vector <int> nextNodeVec;
    //store the sequences (ID);
    vector <int> seq_ID;
    //store the intervalID;
    int intervalID;
    //depth in the BFS tree of the graph;
    int depth;
    Node();
    Node(string _base, int _id);
};

void TopoSortHelp(vector<Node> &graph, bool visit[], int i, vector<int> &Stack);

void TopoSort (vector<Node> &graph);

void printGraph(vector<Node> &graph);

void initialGraph(string &seq, vector<Node> &seqGraph, int k, int id);

//check if an integer x exists in a vector vec;
bool checkExist(vector<int> &vec, int x);

void updateGraph (vector<Node> &newSeqGraph, vector<Node> &graph, vector<Match> &sdpMatches, int k);

//Extract seq given an ID;
void extactSeq(vector<Node> &graph, int seqID, string &seq);

//extract matches of two given seqs from the graph;
void extractMatch(vector<Node> &graph, int firstID, int secondID, 
                    int k, vector<Match> &matchesTwo);
                    
//find next node of given seq ID;
int seqNextPos(vector<Node> &graph, int seqID, int seqCurPos);

//find starting node of given seq ID;
int seqStartPos(vector<Node> &graph, int seqID);

//extract and print any pairwise alignment;
void printPairAlign(vector<Node> &graph, int firstID, int secondID, int k);

#endif
