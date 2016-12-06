#!/bin/bash
./StereoGene cfg=conf.cfg H3K4me1.bed H3K4me3.bed 
Rscript H3K4me1~H3K4me3_report.r
