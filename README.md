network-histogram-code
======================

Code to implement the network histogram (Olhede and Wolfe, arXiv:1312.5306)

Demo: 

1. Save nethist.m and polblogs.mat to the same directory.  

2. Open Matlab and change to this directory.

3. Type or paste in the following sequence of commands:

clear all

load polblogs

A = full(polblogsAdjMat);

help nethist

idx = nethist(A,72);

4. The variable idx now contains the bin indices of each network node.
