#!/bin/bash

epic2 --treatment h3k36me1.bam --control h3.bam --chromsizes chromsizes > h3k36me1_over_h3.peaks.bed
epic2 --treatment h3k36me3.bam --control h3.bam --chromsizes chromsizes > h3k36me3_over_h3.peaks.bed
epic2 --treatment h3k37me1.bam --control h3.bam --chromsizes chromsizes > h3k37me1_over_h3.peaks.bed
