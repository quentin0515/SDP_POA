#include <iostream>
#include <vector>

#include "alignedRead.h"

using namespace std;

void IndexBamRead(bam1_t *alignedRead, AlignedRead &read) {
	//
	// Convert from bam to read.
	//
	read.startPos=alignedRead->core.pos;
	read.endPos=read.startPos;

		 
	read.length = alignedRead->core.l_qseq;
	read.seq=new char[read.length];
	read.name = string(bam_get_qname(alignedRead));


	uint8_t *q = bam_get_seq(alignedRead);
	for (int i=0; i < read.length; i++) {
		read.seq[i]=seq_nt16_str[bam_seqi(q,i)];	
	}
	uint32_t *cigarPtr = bam_get_cigar(alignedRead);

	int leftClip = 0;
	read.alignStart = 0;
	if (alignedRead->core.n_cigar > 0 and bam_cigar_op(cigarPtr[0]) == BAM_CSOFT_CLIP) {
		read.alignStart = bam_cigar_oplen(cigarPtr[0]);
	}
	else if (alignedRead->core.n_cigar > 1 and bam_cigar_op(cigarPtr[2]) == BAM_CSOFT_CLIP) {
		read.alignStart = bam_cigar_oplen(cigarPtr[1]);
	}

	read.alignEnd = read.alignStart;
	for (int i = 0; i < alignedRead->core.n_cigar; i++) {
		/* bam_cigar_type returns a bit flag with:
		*   bit 1 set if the cigar operation consumes the query
		*   bit 2 set if the cigar operation consumes the reference
		*/
	    //if consume reference
		if (bam_cigar_type(cigarPtr[i]) & 0x2 ) {
			read.endPos += bam_cigar_oplen(cigarPtr[i]);
		}
	    //if consume query
		if (bam_cigar_type(cigarPtr[i]) & 1 ) {
			read.alignEnd += bam_cigar_oplen(cigarPtr[i]);
		}
	}
	read.aln=bam_dup1(alignedRead);
}