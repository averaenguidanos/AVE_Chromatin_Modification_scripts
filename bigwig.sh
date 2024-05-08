#!/bin/bash

bamCoverage -b h3.bam --normalizeUsing RPKM -o h3.bw
bamCoverage -b h3k36me1.bam --normalizeUsing RPKM -o h3k36me1.bw
bamCoverage -b h3k36me2.bam --normalizeUsing RPKM -o h3k36me2.bw
bamCoverage -b h3k36me3.bam --normalizeUsing RPKM -o h3k36me3.bw
bamCoverage -b h3k37me1.bam --normalizeUsing RPKM -o h3k37me1.bw
