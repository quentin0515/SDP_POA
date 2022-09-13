#include <string>
#include <vector>
//#include <stack>
#include <assert.h>
#include <list>

#include "graph.h"
#include "naiveSDP.h"

using namespace std;

//TODO: modify cpp files with classes;

Node::Node() = default;
Node::Node(string _base, int _id) : base(_base) {
    seq_ID.push_back(_id);
}

void TopoSortHelp(vector<Node> &graph, bool visit[], int i, vector<int> &Stack) {
    visit[i] = true;
    vector<int> children = graph[i].nextNodeVec;
    for (int j = 0; j < children.size(); ++j) {
        if (!visit[children[j]])
            TopoSortHelp(graph, visit, children[j], Stack);
    }
    Stack.push_back(i);
}

//TODO: change TopoSort (and other functions) to be functions of graph class maybe;
void TopoSort(vector<Node> &graph){
    //use a vector implemented Stack to store visited nodes;
    vector<int> Stack;
    //store visited information;
    bool *visit = new bool[graph.size()];
    for (int i = 0; i < graph.size(); ++i)
        visit[i] = false;

    //starting from the first node of the previous graph, recursively sort nodes;
    for (int i = 0; i < graph.size(); ++i) {
        if (!visit[i]) {
            TopoSortHelp(graph, visit, i, Stack);
        }
    }

    //pop out the stored vector Stack, this is the topoorder;
    //test
    //cout << endl;
    //for (int i = 0; i < Stack.size(); ++i) {
    //    cout << " " << Stack[Stack.size() - 1 - i];
    //}
    cout << endl;

    int curPosition;
    vector<int> curNextNodeVec, curPreNodeVec;

    //TODO: graph programming techniques;
    //clear all nextNodeVec first;
    for (vector<Node>::iterator it = graph.begin(); it != graph.end(); ++it)
        it->nextNodeVec.clear();
    //change all nextNodeVec;
    //O(n^2), n: the no. of nodes;
    for (vector<int>::reverse_iterator rit_Stack = Stack.rbegin(); rit_Stack != Stack.rend(); ++rit_Stack) {
        curPosition = *rit_Stack;
        curPreNodeVec = graph[curPosition].preNodeVec;
        //change nextNodeVec by curPreNodeVec;
        for (vector<int>::iterator it_Pre = curPreNodeVec.begin(); it_Pre != curPreNodeVec.end(); ++it_Pre) {
            graph[*it_Pre].nextNodeVec.push_back(rit_Stack - Stack.rbegin());
        }
    }

    //update graph to TopoSorted
    Node curNode;
    vector<Node> topoSortGraph;
    while (!Stack.empty()) {
        curPosition = Stack.back();
        curNode = graph[curPosition];
        topoSortGraph.push_back(curNode);
        Stack.pop_back();
    }
    graph = topoSortGraph;
    delete[] visit;
    
    //clear all preNodeVec first;
    for (vector<Node>::iterator it = graph.begin(); it != graph.end(); ++it)
        it->preNodeVec.clear();
    //change all preNodeVec;
    //O(n^2), n: the no. of nodes;
    for (vector<Node>::iterator it_graph = graph.begin(); it_graph != graph.end(); ++it_graph) {
        curNextNodeVec = it_graph->nextNodeVec;
        //change nextNodeVec by curPreNodeVec;
        for (vector<int>::iterator it_Next = curNextNodeVec.begin(); it_Next != curNextNodeVec.end(); ++it_Next) {
            graph[*it_Next].preNodeVec.push_back(it_graph - graph.begin());
        }
    }
}

void printGraph(vector<Node> &graph) {
    int count = 0;
    cout << endl << "Graph: " << endl;
    cout << "nextNode info: ";
    for (vector<Node>::iterator it = graph.begin(); it != graph.end(); ++it) {
        cout << count << ")" << "  ";
        for (vector<int>::iterator it_in = (*it).nextNodeVec.begin(); it_in != (*it).nextNodeVec.end(); ++it_in)
            cout << *it_in << " ";
        count++;
    }
    cout << endl;

    count = 0;
    cout << "preNode info: ";
    for (vector<Node>::iterator it = graph.begin(); it != graph.end(); ++it) {
        cout << count << ")" << "  ";
        for (vector<int>::iterator it_in = (*it).preNodeVec.begin(); it_in != (*it).preNodeVec.end(); ++it_in)
            cout << *it_in << " ";
        count++;
    }
    cout << endl;

    count = 0;
    cout << "ID info: ";
    for (vector<Node>::iterator it = graph.begin(); it != graph.end(); ++it) {
        cout << count << ")" << "  ";
        for (vector<int>::iterator it_in = (*it).seq_ID.begin(); it_in != (*it).seq_ID.end(); ++it_in)
            cout << *it_in << " ";
        count++;
    }
    cout << endl;
}
    

void initialGraph(string &seq, vector<Node> &seqGraph, int k, int id) {
    assert(seq.size() >= k);
    seqGraph.clear();
    //Question: use new or not? how to delete? add destructor in class Node???
    //Node *curNode = new Node(seq[0]);
    seqGraph.push_back(Node(seq.substr(0, k), id));
    //build a line graph by adding edges (de bruijn like graph);
    for (int i = 1; i < seq.size() - (k - 1); ++i) {
        seqGraph.push_back(Node(seq.substr(i, k), id));
        seqGraph[i-1].nextNodeVec.push_back(i);
        seqGraph[i].preNodeVec.push_back(i-1);
    }
}

//check if int x exists in vec;
//TODO: improve time;
bool checkExist(vector<int> &vec, int x) {
    bool exist = false;
    for (vector<int>::iterator it = vec.begin(); it != vec.end(); ++it) {
        if (x == *it) {
            exist = true;
            break;
        }
    }
    return exist;
}

void updateGraph (vector<Node> &newSeqGraph, vector<Node> &graph, vector<Match> &sdpMatches, int k) {
	//sdpMatches: may not be a whole match in a case of being visl or visa;
    //if a node matched, merge;
    //if not matched, push_back a new node to graph;
    //note that all_matched = true only when k consecutive bases matched;
    //in other words, node ends at current position matched;
    bool all_matched = false;
    //number of remained matched nodes for current match;
    int remain_no_match_node;
    //the position (in graph) of current matched node;
    int matched_pos;
    //the previous node position on the graph;
    int preNode;

    if ((*(sdpMatches.rbegin())).position.second == 0) {
        //all_matched is true at position 0;
        all_matched = true;
        remain_no_match_node = findMin<int>((*(sdpMatches.rbegin())).len_cut, (*(sdpMatches.rbegin())).len - k + 1);
        matched_pos = (*(sdpMatches.rbegin())).position.first;
    }

    //if first k positions matched;
    if (all_matched) {
        //add sequence id to matched node in graph
        graph[matched_pos].seq_ID.push_back(newSeqGraph[0].seq_ID[0]);
        preNode = matched_pos;
        --remain_no_match_node;
    }
    else {
        //Note the size of graph changes;
        graph.push_back(Node(newSeqGraph[0].base, newSeqGraph[0].seq_ID[0]));
        preNode = graph.size() - 1;
    }

    //loop through all nodes on newSeqGraph;
    for (int i = 1; i < newSeqGraph.size(); ++i) {
        //all_matched = true;
        if (all_matched) {
            if (remain_no_match_node > 0) {
                ++matched_pos;
                //Question: minimize time for exist()???Is it too slow
                if (!checkExist(graph[matched_pos].preNodeVec, preNode)) 
                    graph[matched_pos].preNodeVec.push_back(preNode);
                if (!checkExist(graph[preNode].nextNodeVec, matched_pos))
                    graph[preNode].nextNodeVec.push_back(matched_pos);
                graph[matched_pos].seq_ID.push_back(newSeqGraph[i].seq_ID[0]);
                preNode = matched_pos;
                --remain_no_match_node;
            }
            else {
                all_matched = false;
                for (vector<Match>::reverse_iterator rit = sdpMatches.rbegin(); rit != sdpMatches.rend(); ++rit) {
                    //newSeqGraph seq2, is on the y axis;
                    //TODO: add search method;
                    if (i == (*rit).position.second) {
                        all_matched = true;
                        remain_no_match_node = findMin<int>((*rit).len_cut, (*rit).len - k + 1);
                        matched_pos = (*rit).position.first;
                        break;
                    }
                }
                //a new match start;
                if (all_matched) {
                    if (!checkExist(graph[matched_pos].preNodeVec, preNode)) 
                        graph[matched_pos].preNodeVec.push_back(preNode);
                    if (!checkExist(graph[preNode].nextNodeVec, matched_pos))
                        graph[preNode].nextNodeVec.push_back(matched_pos);
                    graph[matched_pos].seq_ID.push_back(newSeqGraph[i].seq_ID[0]);
                    preNode = matched_pos;
                    --remain_no_match_node;
                }
                else {
                    graph.push_back(Node(newSeqGraph[i].base, newSeqGraph[i].seq_ID[0]));
                    graph[graph.size() - 1].preNodeVec.push_back(preNode);
                    graph[preNode].nextNodeVec.push_back(graph.size() - 1);
                    preNode = graph.size() - 1;
                }
            }
        }
        //previous node is not matched;
        else {
            //all_matched = false;
            for (vector<Match>::reverse_iterator rit = sdpMatches.rbegin(); rit != sdpMatches.rend(); ++rit) {
                if (i == (*rit).position.second) {
                    all_matched = true;
                    remain_no_match_node = findMin<int>((*rit).len_cut, (*rit).len - k + 1);
                    matched_pos = (*rit).position.first;
                    break;
                }
            }
            //a new match start;
            if (all_matched) {
                if (!checkExist(graph[matched_pos].preNodeVec, preNode)) 
                    graph[matched_pos].preNodeVec.push_back(preNode);
                if (!checkExist(graph[preNode].nextNodeVec, matched_pos))
                    graph[preNode].nextNodeVec.push_back(matched_pos);
                graph[matched_pos].seq_ID.push_back(newSeqGraph[i].seq_ID[0]);
                preNode = matched_pos;
                --remain_no_match_node;
            }
            else {
                graph.push_back(Node(newSeqGraph[i].base, newSeqGraph[i].seq_ID[0]));
                graph[graph.size() - 1].preNodeVec.push_back(preNode);
                graph[preNode].nextNodeVec.push_back(graph.size() - 1);
                preNode = graph.size() - 1;
            }
        }
    }
}

//Extract seq given an ID;
void extactSeq(vector<Node> &graph, int seqID, string &seq) {
    int seqStart = seqStartPos(graph, seqID);
    assert(seqStart >= 0);
    int curPos = seqStart;
    //seq.push_back((graph[curPos].base)[0]);
    while (seqNextPos(graph, seqID, curPos) >= 0) {
        seq.push_back((graph[curPos].base)[0]);
        curPos = seqNextPos(graph, seqID, curPos);
    }
    seq += graph[curPos].base;
    /*
    vector<Node>::iterator it = graph.begin();
    for (; it != graph.end(); ++it) {
        if ((it->nextNodeVec).size() > 0)
            graphSeq.push_back((it->base)[0]);
        else
            graphSeq += it->base;
    }
    */
}

//extract matches of two given seqs from the graph;
//TODO: way to be faster: combine functions together;
void extractMatch(vector<Node> &graph, int firstID, int secondID, 
                    int k, vector<Match> &matchesTwo) {
    //get matches for given seqs;
    int firstStart = seqStartPos(graph, firstID);
    int secondStart = seqStartPos(graph, secondID);
    vector<int> matchPosFirst, matchPosSecond;

    //extract match positions for first seq;
    assert(firstStart >= 0);
    int curPos = firstStart;
    //position on the first seq;
    int firstPos = 0;
    while (curPos >= 0) {
        if (checkExist(graph[curPos].seq_ID, secondID))
            matchPosFirst.push_back(firstPos);
        curPos = seqNextPos(graph, firstID, curPos);
        ++firstPos;
    }

    //extract match positions for second seq;
    assert(secondStart >= 0);
    curPos = secondStart;
    //position on the first seq;
    int secondPos = 0;
    while (curPos >= 0) {
        if (checkExist(graph[curPos].seq_ID, firstID))
            matchPosSecond.push_back(secondPos);
        curPos = seqNextPos(graph, secondID, curPos);
        ++secondPos;
    }

    //construct the matches for the two seqs;
    assert(matchPosFirst.size() == matchPosSecond.size());
    //TODO: iterator???
    for (int i = 0; i < matchPosFirst.size(); ++i) {
        matchesTwo.push_back(Match(k, make_pair(matchPosFirst[i], matchPosSecond[i])));
    }
}

//find next node of given seq ID;
int seqNextPos(vector<Node> &graph, int seqID, int seqCurPos) {
    int seqNext = -1;
    vector<int>::iterator it_next = graph[seqCurPos].nextNodeVec.begin();
    for (; it_next != graph[seqCurPos].nextNodeVec.end(); ++it_next) {
        if (checkExist(graph[*it_next].seq_ID, seqID)) {
            seqNext = *it_next;
            break;
        }
    }
    //return -1 if not found;
    return seqNext;
}

//find starting node of given seq ID;
int seqStartPos(vector<Node> &graph, int seqID) {
    //search for the start of seqID seq on the graph: bfs;
    int seqStart = -1;
    bool findStart = false;
    vector<int> graphStart;
    vector<int> graphEnd;
    getStartEnd(graph, graphStart, graphEnd);
    //use list as a queue: store node positions in graph;
    list<int> que;
    //store visit information;
    vector<bool> visit(graph.size(), false);
    //initialization
    vector<int>::iterator it_graphStart = graphStart.begin();
    for (; it_graphStart != graphStart.end(); ++it_graphStart) {
        if (checkExist(graph[*it_graphStart].seq_ID, seqID)) {
            seqStart = *it_graphStart;
            return seqStart;
        }
        que.push_back(*it_graphStart);
        visit[*it_graphStart] = true;
    }

    int curNode;
    while (!que.empty()) {
        curNode = que.front();
        que.pop_front();
        vector<int>::iterator it_next = graph[curNode].nextNodeVec.begin();
        for (; it_next != graph[curNode].nextNodeVec.end(); ++it_next) {
            if (!visit[*it_next]) {
                if (checkExist(graph[*it_next].seq_ID, seqID)) {
                    seqStart = *it_next;
                    return seqStart;
                }
                que.push_back(*it_next);
                visit[*it_next] = true;
            }
        }
    }

    //return -1 if not found;
    return seqStart;
}

//prepare print of insetion/deletion of firstSeq, or match;
void preparePrint(int option, int start, int end, 
                    string &seq1, string &seq2, 
                    string &str1, string &str2, string &strBtw) {
    assert(option >= 0 && option <= 2);
    //insetion on seq1;
    if (option == 1) {
        str1 += seq1.substr(start, end - start + 1);
        string tem2(end - start + 1, '_');
        str2 += tem2;
        string temBtw(end - start + 1, ' ');
        strBtw += temBtw;
    }
    //insetion on seq2;
    else if (option == 2) {
        str2 += seq2.substr(start, end - start + 1);
        string tem1(end - start + 1, '_');
        str1 += tem1;
        string temBtw(end - start + 1, ' ');
        strBtw += temBtw;
    }
    //match;
    else {
        //start end are on seq1;
        str1 += seq1.substr(start, end - start + 1);
        str2 += seq1.substr(start, end - start + 1);;
        string temBtw(end - start + 1, '|');
        strBtw += temBtw;
    }
    
}


//extract and print any pairwise alignment;
void printPairAlign(vector<Node> &graph, int firstID, int secondID, int k) {
    //extract seq;
    string firstSeq, secondSeq;
    extactSeq(graph, firstID, firstSeq);
    extactSeq(graph, secondID, secondSeq);

    //extract matches for given two seqs;
    vector<Match> matchesTwo;
    extractMatch(graph, firstID, secondID, k, matchesTwo);

    //preapare printing seqs of alignment according to matches;
    string str1, str2, strBtw;
    //TODO: refer to LRA printing method and improve;
    int curFirstPos = 0;
    int curSecondPos = 0;
    Match curMatch, nextMatch;
    vector<Match>::iterator it_matchesTwo = matchesTwo.begin();
    for (; it_matchesTwo != matchesTwo.end(); ++it_matchesTwo) {
        //idea: each iteration print to the end of curMatch;
        curMatch = *it_matchesTwo;
        if (it_matchesTwo + 1 != matchesTwo.end()) {
            nextMatch = *(it_matchesTwo + 1);
            //precedes;
            if (curMatch.position.first + curMatch.len <= nextMatch.position.first &&
                curMatch.position.second + curMatch.len <= nextMatch.position.second) {
                preparePrint(1, curFirstPos, curMatch.position.first - 1, firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "1st: " << firstSeq.substr(curFirstPos, curMatch.position.first - curFirstPos) << endl;
                preparePrint(2, curSecondPos, curMatch.position.second - 1, firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "2ns: " << secondSeq.substr(curSecondPos, curMatch.position.second - curSecondPos) << endl;
                preparePrint(0, curMatch.position.first, curMatch.position.first + curMatch.len - 1, 
                                firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "Match: " << firstSeq.substr(curMatch.position.first, curMatch.len) << endl;
                curFirstPos = curMatch.position.first + curMatch.len;
                curSecondPos = curMatch.position.second + curMatch.len;
            }
            //visl;
            else if (diagDiff(nextMatch, curMatch) >= 0) {
                preparePrint(1, curFirstPos, curMatch.position.first - 1, firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "1st: " << firstSeq.substr(curFirstPos, curMatch.position.first - curFirstPos) << endl;
                preparePrint(2, curSecondPos, curMatch.position.second - 1, firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "2ns: " << secondSeq.substr(curSecondPos, curMatch.position.second - curSecondPos) << endl;
                preparePrint(0, curMatch.position.first, nextMatch.position.first - 1, 
                                firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "Match: " 
                //     << firstSeq.substr(curMatch.position.first, nextMatch.position.first - curMatch.position.first) 
                //     << endl;
                curFirstPos = curMatch.position.first + nextMatch.position.first - curMatch.position.first;
                curSecondPos = curMatch.position.second + nextMatch.position.first - curMatch.position.first;
            }
            //visa
            else if (diagDiff(nextMatch, curMatch) < 0) {
                preparePrint(1, curFirstPos, curMatch.position.first - 1, firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "1st: " << firstSeq.substr(curFirstPos, curMatch.position.first - curFirstPos) << endl;
                preparePrint(2, curSecondPos, curMatch.position.second - 1, firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "2ns: " << secondSeq.substr(curSecondPos, curMatch.position.second - curSecondPos) << endl;
                preparePrint(0, curMatch.position.first, 
                                curMatch.position.first + nextMatch.position.second - curMatch.position.second - 1, 
                                firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "Match: " 
                //     << firstSeq.substr(curMatch.position.first, nextMatch.position.second - curMatch.position.second) 
                //     << endl;
                curFirstPos = curMatch.position.first + nextMatch.position.second - curMatch.position.second;
                curSecondPos = curMatch.position.second + nextMatch.position.second - curMatch.position.second;
            }
        }
        //curMatch is the last match
        else {
            preparePrint(1, curFirstPos, curMatch.position.first - 1, firstSeq, secondSeq, str1, str2, strBtw);
            //cout << "1st: " << firstSeq.substr(curFirstPos, curMatch.position.first - curFirstPos) << endl;
            preparePrint(2, curSecondPos, curMatch.position.second - 1, firstSeq, secondSeq, str1, str2, strBtw);
            //cout << "2ns: " << secondSeq.substr(curSecondPos, curMatch.position.second - curSecondPos) << endl;
            preparePrint(0, curMatch.position.first, curMatch.position.first + curMatch.len - 1, 
                            firstSeq, secondSeq, str1, str2, strBtw);
            //cout << "Match: " << firstSeq.substr(curMatch.position.first, curMatch.len) << endl;
            curFirstPos = curMatch.position.first + curMatch.len;
            curSecondPos = curMatch.position.second + curMatch.len;
            if (curFirstPos < firstSeq.size())
                preparePrint(1, curFirstPos, firstSeq.size() - 1, firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "1st: " << firstSeq.substr(curFirstPos, firstSeq.size() - curFirstPos) << endl;
            if (curSecondPos < secondSeq.size())
                preparePrint(2, curSecondPos, secondSeq.size() - 1, firstSeq, secondSeq, str1, str2, strBtw);
                //cout << "2ns: " << secondSeq.substr(curSecondPos, secondSeq.size() - curSecondPos) << endl;
            cout << "<<<<<<<<<<<<<<<<<" << endl;
        }
    }

    //no of char each line;
    int noCharLine = 200;
    //remained length of printing;
    //int sizePrint = seq1.size();
    int curPos = 0;
    //print;
    assert(str1.size() == str2.size() && str1.size() == strBtw.size());
    while (curPos + noCharLine <= str1.size()) {
        cout << str1.substr(curPos, noCharLine) << endl;
        cout << strBtw.substr(curPos, noCharLine) << endl;
        cout << str2.substr(curPos, noCharLine) << endl;
        curPos += noCharLine;
    }
    cout << str1.substr(curPos, str1.size() - curPos) << endl;
    cout << strBtw.substr(curPos, strBtw.size() - curPos) << endl;
    cout << str2.substr(curPos, str2.size() - curPos) << endl;
    cout << "END<<<<<<<<<<<<<<" << endl;

}
