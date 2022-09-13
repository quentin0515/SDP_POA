#ifndef __ALIGNEDREAD_H_INCLUDED__ //   if naiveSDP.h hasn't been included yet
#define __ALIGNEDREAD_H_INCLUDED__ //   #define this so the compiler knows it has been included

#include "htslib/hts.h"
#include "htslib/sam.h"
#include <string>
#include <fstream>
#include <vector>

using namespace std;

class AlignedRead {
public:
	int index;
	string chrom;
	//alignment starting position on read
	int alignStart, alignEnd;
	//read starting position on reference
	int startPos;
	int endPos;
	int length;
	char* seq, *seqrc;
	string name;
	int n;
	bam1_t *aln;	

	~AlignedRead() {
		if (seq != NULL) {
			delete seq;
		}
	}

	int LookupRefPositionInRead(int refPos) {
		int curReadPos = alignStart;
		int curRefPos = startPos;
    	int k, l;
		uint32_t *cigar = bam_get_cigar(aln);
		if (refPos < startPos) {
			return -1;
		}
		if (refPos > endPos ) {
			return -1;
		}
    	for (k = 0; k < aln->core.n_cigar; ++k) {
			int opLen = bam_cigar_oplen(cigar[k]);
			int op = bam_cigar_op(cigar[k]);
			if (op == BAM_CSOFT_CLIP) {
				curReadPos = opLen;
				continue;
			}
			if ((op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF) and 
					curRefPos <= refPos and curRefPos + opLen > refPos) {
				return curReadPos + (refPos-curRefPos);
			}
			else if (op == BAM_CDEL and curRefPos <= refPos and curRefPos + opLen > refPos) {
				return curReadPos;
			}
			if (bam_cigar_type(cigar[k]) & 1 ) {
				curReadPos += opLen;
			}
			if (bam_cigar_type(cigar[k]) & 2 ){
				curRefPos += opLen;
			}
		}
		return -1;
	}
};

void IndexBamRead(bam1_t *alignedRead, AlignedRead &read);

#endif