//test file: /home/cmb-16/mjc/quentin/testksw/racon/build/bin/test/w4000/h1.reads.to_cons.bam
// /home/cmb-16/mjc/quentin/testksw/racon/build/bin/test/w2000e02/h1.reads.to_cons.bam
//https://github.com/quentin0515/sdp_poa.git

#include <iostream>
#include <string>
#include <assert.h>
#include <math.h>
#include <vector>

#include "./naiveSDP.h"
#include "./graph.h"
#include "alignedRead.h"
#include "htslib/hts.h"
#include "htslib/sam.h"


using namespace std;

//TODO: change updateGraph according to len_cut;
//TODO: make sure fused nodes have correct seq ID;
//write function to get seq given two ID;

//TODO: clean code;
//1. separate modes: including ifMatch changes;
//2. add function bfs: to be used;

int main(int argc, char *argv[])
{
    //string seq2 = "gttagagcggcccaatcttaaaatgatattaacacgtatgcgagtcagtctaatccgttttaaagcgacttcaggtctaggccc";
    string seq4 = "gccgcccttccctacctcaccccgttgttacaacctatcccgccacacccttttccggtgaaagatggaatacacgcgggggaacttacgtccatgggag";
    string seq5 = "agccgatcgaatcgctgataaaaacacgccattggttccgatcagtcaacccgctgcacgagcgaatatgcagtagttctaactataggtctacgttgta";
    string seq3 = "gtttcagcgtcccata";
    string seq2 = "gccgtctcctacc";
    string seq1 = "agccgagtcatac";

    //k-mer length;
    int k = 3;
    //SDP final cost;
    int cost;
    //current sequence ID
    int id = 0;
    //store graph: the vector order is the topological order
    vector<Node> graph;
    //store graph of a new seq: a line structure
    vector<Node> newSeqGraph;
    //sdpMatches: in a reversed order: the last selected match is at the begining;
    vector<Match> sdpMatches;

    ++id;
    initialGraph(seq1, newSeqGraph, k, id);
    //note the preNode (nextNode) of the first (last) node is undefined;
    graph = newSeqGraph;

    //naiveSDP(seq2, sdpMatches, k, graph, cost);
    naiveSDP_1(seq2, sdpMatches, k, graph, cost);

    ////if (sdpMatches.size() > 0) {
    ////    cout << "SDP cost: " << cost << endl;
    ////}

    //given sdp results: a vector of selected matches, update POA graph
    ++id;
    initialGraph(seq2, newSeqGraph, k, id);

    //test
    cout << "original graph:" << endl;
    printGraph(graph);
    updateGraph(newSeqGraph, graph, sdpMatches, k);
    //TODO:test updateGraph();

    //test
    cout << "after updateGraph:" << endl;
    printGraph(graph);

    //output a vector: graph; vector order is the topo order, each node stores the graph structure;
    TopoSort(graph);
    //test
    cout << "after TopoSort:" << endl;
    printGraph(graph);

    //SDP again with graph and seq3;
    naiveSDP_1(seq3, sdpMatches, k, graph, cost);
    ++id;
    initialGraph(seq3, newSeqGraph, k, id);

    updateGraph(newSeqGraph, graph, sdpMatches, k);
    TopoSort (graph);
    ////if (sdpMatches.size() > 0) {
    ////    cout << "SDP cost: " << cost << endl;
    ////}

    //test
    cout << "after TopoSort:" << endl;
    printGraph(graph);


    /////////////////////////////////////////////////////////////////////////////////
    
    cout << endl 
         << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl 
         << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl
         << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl
         << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl
         << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl
         << endl;
    //test on real reads;
    k = 15;
    const char* testFile = "/home/cmb-16/mjc/quentin/testksw/racon/build/bin/test/w4000/h1.reads.to_cons.bam";
    //the starting reads ID;
    int startRead = 200;
    int countRead = 0;
    htsFile *htsfp;
	bam_hdr_t *samHeader;

    htsfp = hts_open(testFile, "r");
    samHeader = sam_hdr_read(htsfp);
	bam1_t *alignedRead = bam_init1();

    //loop while there is a new read and countRead < startRead;
    //find a random position in bam file to perform SDP_POA;
    while (sam_read1(htsfp, samHeader, alignedRead) >= 0 and countRead < startRead) {
        ++countRead;
    }

    //test
    cout << endl << "break countRead: " << countRead << endl;

    //the number of sequences to align;
    int noReadAlign = 5;
    countRead = 0;

    AlignedRead* curRead = new AlignedRead;
	int res = sam_read1(htsfp, samHeader, alignedRead);
    //TODO: why here?;
	while (res >= 0 and alignedRead->core.flag & 2304) {
		res = sam_read1(htsfp, samHeader, alignedRead);
	}

    if (res < 0) {
		cerr << "NO reads" << endl;
		exit(0);
	}
	
    IndexBamRead(alignedRead, *curRead);
    ++countRead;
    //index starting from 1;
    curRead -> index = countRead;

    //extract seq stored in read;
    string temp1(curRead -> seq);
    //test
    seq1 = temp1.substr(curRead->alignStart, 
            (curRead->alignEnd) - (curRead->alignStart) + 1);
    //test
    cout << "length: " << (curRead->alignEnd) - (curRead->alignStart) + 1 << endl;
    //seq1 = temp1;
    //cout << endl << seq1[0] << endl;
    initialGraph(seq1, newSeqGraph, k, curRead -> index);
    graph = newSeqGraph;
    //test
    //cout << endl << newSeqGraph[0].base << endl;

    while (countRead < noReadAlign) {
        AlignedRead* curRead = new AlignedRead;
        //TODO: double check get read from bam files;
	    int res = sam_read1(htsfp, samHeader, alignedRead);
	    while (res >= 0 and alignedRead->core.flag & 2304) {
		    res = sam_read1(htsfp, samHeader, alignedRead);
	    }
        if (res < 0) {
		    cerr << "NO reads" << endl;
		    exit(0);
	    }
        IndexBamRead(alignedRead, *curRead);
        ++countRead;
        curRead -> index = countRead;

        //extract seq stored in read;
        string temp2(curRead -> seq);
        //test
        seq2 = temp2.substr(curRead->alignStart, 
                (curRead->alignEnd) - (curRead->alignStart) + 1);
        //test
        cout << "length: " << (curRead->alignEnd) - (curRead->alignStart) + 1 << endl;
        //seq2 = temp2;

        naiveSDP_1(seq2, sdpMatches, k, graph, cost);

        initialGraph(seq2, newSeqGraph, k, curRead -> index);

        updateGraph(newSeqGraph, graph, sdpMatches, k);

        TopoSort(graph);
    }
    
    printPairAlign(graph, 2, 3, k);

}