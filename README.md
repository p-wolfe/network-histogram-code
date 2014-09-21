network-histogram-code
======================

Code to implement the network histogram (Olhede and Wolfe, arXiv:1312.5306)

Windows executable and Matlab files.  Demos below use the political blogs dataset
of Adamic and Glance (2005), see http://www-personal.umich.edu/~mejn/netdata/.
First 586 blogs are liberal; remaining 638 are conservative.  Nodes with zero degree
have been removed.

Matlab demo:

1. Save nethist.m and polblogs.mat to the same directory. 
2. Open Matlab and change to this directory. 
3. Type or paste in the following sequence of commands: 

clear all  
load polblogs  
A = full(polblogsAdjMat);
help nethist  
idx = nethist(A,72); % or nethist('polblogs.txt',72,'polblogs_out.txt');

The variable idx now contains the bin indices of each network node.

A quick visualization of the result may be obtained by typing

[~,ind] = sort(idx);  
spy(A(ind,ind));  

Standalone Windows demo:

1. Save nethist.exe, nethist.bat, and polblogs.txt to the same directory.
2. Make sure you have MATLAB Compiler Runtime (MCR) version 7.16, 
freely available from http://www.mathworks.com/products/compiler/mcr/
3. Double-click on nethist.bat, OR open up a command window and type:

nethist polblogs.txt 72 polblogs_out.txt
