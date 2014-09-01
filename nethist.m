function idx = nethist(A,h)
%NETHIST Network histogram.
%   IDX = NETHIST(A, H) computes the network histogram of an N-by-N
%   adjacency matrix A, which may be sparse or dense but must be 0-1
%   valued, symmetric, and with zero main diagonal.  Optional bandwidth
%   parameter H specifies the number of nodes in each histogram bin, which
%   is automatically determined if H is not specified.  NETHIST returns an
%   N-by-1 vector IDX containing the bin indices of each network node.
%   
%   Copyright (C) 2013 Sofia C. Olhede and Patrick J. Wolfe (arXiv:1312.5306)
%   NETHIST comes with ABSOLUTELY NO WARRANTY; for details type `TYPE NETHIST'.
%   This is free software, and you are welcome to redistribute it
%   under certain conditions; type `TYPE NETHIST' for details.

%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or (at
%     your option) any later version.
% 
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%     02110-1301 USA.

% Error-check inputs
assert(ismatrix(A),'Input A must be a matrix');
assert(all(size(A)==size(A')),'Matrix input A must be square');
assert(all(size(A)>1),'Matrix input A must be of dimension at least 2 x 2');
assert(~norm(A-A','fro'),'Matrix input A must be symmetric');
assert(norm(A-A.*A,'fro')==0,'All entries of A must be 0 or 1');

% Compute necessary summaries from A
n = size(A,2);
rhoHat = sum(A(:))/(n*(n-1));

% Create random number stream with new 'shuffled' seed based on system clock
rStreamArray = cell(1,1);
[rStreamArray{:}] = RandStream.create('mrg32k3a','NumStreams',1,'Seed','shuffle');
rngSeed = rStreamArray{1}.Seed;
rStream = rStreamArray{1};
display(['Generated seed based on system clock: ' int2str(rngSeed)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick an analysis bandwidth and initialize via regularized spectral clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('h','var'), % use data-driven h if optional input not specified
    c = min(4,sqrt(n)/8);
    [h,~] = oracbwplugin(A,c,'degs',1); % fractional h, prior to rounding
    display(['Determining bandwidth from data; initial estimate ' num2str(h)]);
else
    display(['Determining bandwidth from user input; intially set to ' num2str(h)]);
end
h = max(2,min(n,round(h)));
display(['Final bandwidth ' num2str(h)]);
lastGroupSize = mod(n,h);
while (lastGroupSize==1) && (h>2) % step down hbarHist, to avoid singleton final group
    h = h - 1;
    lastGroupSize = mod(n,h);
    display('NB: Bandwidth reduced to avoid singleton group');
end
display(['Adjacency matrix has ' int2str(n) ' rows/cols'])

% Initialize using regularized spectral clustering based on row similarity
tstart = tic;
regParam = rhoHat/4;

distVec = pdist(A+regParam,'hamming'); % vector of pairwise distances between regularized rows of A
L = 1 - squareform(distVec); % exponential Taylor approximation to L_ij = exp(-||A_i. - A_j.||^2 / 2) for small ||.||
clear distVec;
d = sum(L,2);
OPTS.vo = ones(n,1);
[u,~] = eigs( (d.^(-1/2)*(d.^(-1/2))').*L-sqrt(d)*sqrt(d)'/sqrt(d'*d),1,'LA',OPTS); % 2nd eigenvector of normalized Laplacian
clear L
u = u.*sign(u(1)); %  set 1st coord >= 0 wlog, to fix an arbitrary sign permutation
[~,ind]=sort(u,'ascend'); % sort on this embedding
k = ceil(n/h);

idxInit = zeros(n,1);
for i = 1:k,
    idxInit(ind(((i-1)*h+1):min(n,i*h))) = i;
end
display(['Initial label vector assigned from row-similarity ordering; time ' num2str(toc(tstart))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute histogram from adjacency matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = graphest_fastgreedy(A,h,idxInit,rStream);

% END MAIN FUNCTION NETHIST.M

function [h, estMSqrd] = oracbwplugin(A,c,type,alpha)
%ORACBWPLUGIN - Oracle bandwidth plug-in estimtor for network histograms 
%   [h estMSqrd] = oracbwplugin(A,c,type,alpha) returns a plug-in estimate
%   of the optimal histogram bandwidth (blockmodel community size) as a 
%   function of the following inputs:
%
%  A - input adjacency matrix (must be a simple random graph)
%  c - positive multiplier by which to estimate slope from +/- c sqrt(n)
%  type - Estimate slope from sorted vector (optional; 'degs' [default] or 'eigs')
%  alpha = Holder exponent (optional; defaults to 1)
%
%   Example:
%   load polblogs
%   A = full(polblogsAdjMat);
%   h = oracbwplugin(A,3,'eigs',1); % returns h = 73.5910
%   h = oracbwplugin(A,3,'degs',1); % returns h = 74.1031
%   Copyright 2013 Sofia C. Olhede and Patrick J. Wolfe (arXiv:1312.5306)

showPlots = 1; % Boolean for plotting

% Default inputs
if ~exist('type','var'), type = 'degs'; end
if ~exist('alpha','var'), alpha = 1; end

n = size(A,2);
midPt = round(n/2-c*sqrt(n)):round(n/2+c*sqrt(n));
selfLoops = any(diag(A)~=0); % check for self-loops in A
sampleSize = nchoosek(n+selfLoops,2);
rhoHat = full(sum(A(triu(true(n),1-selfLoops)))) / sampleSize;

% Rank-1 graphon estimate via fhat(x,y) = mult*u(x)*u(y)*pinv(rhoHat);
switch type
    case 'eigs'
        [u,mult] = eigs(A,1);
    case 'degs'
        u = full(sum(A,2));
        mult = (u'*A*u)/(u'*u)^2;
    otherwise
        error('Invalid input ''type''');
end

% Calculate bandwidth
u = sort(u);
uMid = u(midPt(1):midPt(end),1);
p = polyfit((1:length(uMid))',uMid,1); % linear fit to uMid
if alpha~=1, error('Currently only supports alpha = 1'); end
h = full( ( 2^(alpha+1) * alpha * mult^2 * (p(2)+p(1)*length(uMid)/2)^2 * p(1)^2 * pinv(rhoHat) )^(-1/(2*(alpha+1))) );

estMSqrd = full( 2^(alpha+1) * alpha * mult^2 * (p(2)+p(1)*length(uMid)/2)^2 * p(1)^2 * pinv(rhoHat)^2 * (n+1)^2 );
MISEfhatBnd = estMSqrd * ( (2/sqrt(estMSqrd)) * (sampleSize*rhoHat)^(-1/2) + 1/n );
display(['\hat M^2: ' num2str(estMSqrd) '; \hat MISE bound: ' num2str(MISEfhatBnd)]);

% Diagnostic plots
if showPlots==1,
    subplot(1,2,1);
    plot(u);
    title('Graphon projection for bandwidth estimation');
    subplot(1,2,2)
    plot(uMid);
    title('Chosen patch of projection component (adjust using c)');
end


function bestLabelVec = graphest_fastgreedy(A,hbar,inputLabelVec,rStream)
%GRAPHEST_FASTGREEDY implements likelihood-based optimization for nethist.m
%   Copyright 2013 Sofia C. Olhede and Patrick J. Wolfe (arXiv:1312.5306)

absTol = 2.5*10^-4; % absolute tolerance - change in normalized (order-one) likelihood for 3 consec iterations
maxNumRestarts = 500; %;
if size(A,2)<=256, % ~n-choose-2 for n = 256
    allInds = 1; % Consider basing this choice on size (or sparse-vs-dense storage) of A
else
    allInds = 0;
end
if allInds
    numMHsteps = nchoosek(size(A,2),2); %round( nchoosek(size(A,2),2) / 40 );
else
    numMHsteps = 2*10^4; %round(size(A,2).^1.67); %4*10^5; % 10^7 is a typical laptop memory limit
end

% Compute necessary quantities
n = size(A,2);
sampleSize = nchoosek(n,2);

% Initialize cluster assignments in order
smallerLastGroup = ~(mod(n,hbar)==0); % 0 iff h divides n
k = ceil(n/hbar);
display(['Fitting a(n) ' int2str(k) '-group blockmodel']);
equalSizeInds = 1:(k-smallerLastGroup);
orderedLabels = zeros(n,1);
h = zeros(k,1);
orderedClusterInds = zeros(k,hbar);
for a = 1:(k-smallerLastGroup)
    orderedInds = ((a-1)*hbar+1):(a*hbar);
    h(a) = length(orderedInds);
    orderedLabels(orderedInds) = a;
    orderedClusterInds(a,:) = orderedInds;
end
if smallerLastGroup,
    orderedIndsLast = ((k-1)*hbar+1):n;
    h(k) = length(orderedIndsLast);
    orderedClusterInds(k,:) = [orderedIndsLast, zeros(1,hbar-length(orderedIndsLast))];
    display(['Final group of size ' int2str(h(k)) ' smaller than others of size ' int2str(h(1))])
else
    display(['All groups of equal size ' int2str(h(1))])
end
aLeqb = triu(true(k,k),0);
habSqrd = h*h'-diag(h.^2-h.*(h-1)./2);
assert(all(habSqrd(:)>=1),'All clusters must contain at least 2 nodes');
assert(max(orderedClusterInds(k,:))==n,'All nodes must be assigned to a cluster');
assert(sum(h)==n,'Number of cluster assignments must equal number of nodes');

initialLabelVec = inputLabelVec;
initialClusterInds = zeros(k,hbar);
initialClusterCentroids = zeros(k,n);
for a = 1:(k-smallerLastGroup) % update cluster indices
    initialClusterInds(a,:) = find(initialLabelVec==a)';
    initialClusterCentroids(a,:) = sum(A(:,initialClusterInds(a,:)),2)';
end
if smallerLastGroup
    initialClusterInds(k,1:length(find(initialLabelVec==k))) = find(initialLabelVec==k)';
    initialClusterCentroids(k,:) = sum(A(:,initialClusterInds(k,1:length(find(initialLabelVec==k)))),2)';
end
initialACounts = getSampleCounts(A,initialClusterInds);
initialLL = fastNormalizedBMLogLik(initialACounts(aLeqb)./habSqrd(aLeqb),habSqrd(aLeqb),sampleSize);

bestLL = initialLL;
oldNormalizedBestLL = bestLL*2*sampleSize/sum(A(:));
bestLabelVec = initialLabelVec;
bestClusterCentroids = initialClusterCentroids;
bestCount = 0;
consecZeroImprovements = 0;
tolCounter = 0;

tic;
tStartOuter = tic;
for mm = 1:maxNumRestarts
    
    oneTwoVec = 1 + (rand(rStream,numMHsteps,1)>2/3); % 1 wp x; 2 wp 1-x
    iVec = max(1,ceil(rand(rStream,numMHsteps,1)*n));
    jVec = max(1,ceil(rand(rStream,numMHsteps,1)*n));
    kVec = max(1,ceil(rand(rStream,numMHsteps,1)*n)); % random pairs of elements between 1 and n
        
    bestClusterInds = zeros(k,hbar);
    for a = 1:(k-smallerLastGroup) % update cluster indices
        bestClusterInds(a,:) = find(bestLabelVec==a)';
        bestClusterCentroids(a,:) = sum(A(:,bestClusterInds(a,:)),2)';
    end
    if smallerLastGroup
        bestClusterInds(k,1:length(find(bestLabelVec==k))) = find(bestLabelVec==k)';
        bestClusterCentroids(k,:) = sum(A(:,bestClusterInds(k,1:length(find(bestLabelVec==k)))),2)';
    end
    bestACounts = getSampleCounts(A,bestClusterInds);
    bestLL = fastNormalizedBMLogLik(bestACounts(aLeqb)./habSqrd(aLeqb),habSqrd(aLeqb),sampleSize);
    
    currentACounts = bestACounts;
    currentClusterInds = bestClusterInds;
    currentLL = bestLL;
    currentLabelVec = bestLabelVec;
        
    for m = 1:numMHsteps
        
        % Prepare to update quantities for trial clustering
        trialClusterInds = currentClusterInds;
        trialLabelVec = currentLabelVec;
        trialACounts = currentACounts;
        trialLL = currentLL;
        
        % Implement consecutive pairwise swaps to obtain trial clustering
        for mmm = 1:oneTwoVec(m)
            if mmm==1 % Step 1 of 2:
                i = iVec(m); % ideally here i,j are very similar, but in different groups
                j = jVec(m);
                a = trialLabelVec(i); % get group labels of nodes in chosen pair
                b = trialLabelVec(j);
            elseif a~=b % Step 2 of 2: Check that pairwise swap was made in Step 1
                i = jVec(m); % ideally here i,j are very similar, but in different groups
                j = kVec(m);
                a = trialLabelVec(i); % get group labels of nodes in chosen pair
                b = trialLabelVec(j);                    
            end
            if a~=b % Swap and update trial likelihood only if nodes i and j are in different clusters
                
                trialLabelVec(i) = b; % assign to node i the new cluster label previously belonging to node j
                trialLabelVec(j) = a; % assign to node j the new cluster label previously belonging to node i
                
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                habSqrdCola = habSqrd(:,a);
                habSqrdColb = habSqrd(:,b);
                habSqrdEntryab = habSqrd(a,b);
                %
                oldThetaCola = trialACounts(:,a)./habSqrdCola;
                oldThetaColb = trialACounts(:,b)./habSqrdColb;
                oldThetaEntryab = trialACounts(a,b)./habSqrdEntryab;
                %
                oldThetaCola(oldThetaCola<=0) = eps;     % else 0*log(0) evaluates to NaN
                oldThetaCola(oldThetaCola>=1) = 1 - eps; % else (1-1)*log(1-1) evaluates to NaN
                %
                oldThetaColb(oldThetaColb<=0) = eps;     % else 0*log(0) evaluates to NaN
                oldThetaColb(oldThetaColb>=1) = 1 - eps; % else (1-1)*log(1-1) evaluates to NaN
                %
                oldThetaEntryab(oldThetaEntryab<=0) = eps;     % else 0*log(0) evaluates to NaN
                oldThetaEntryab(oldThetaEntryab>=1) = 1 - eps; % else (1-1)*log(1-1) evaluates to NaN
                
                % Begin updating
                trialClusterInds(a,trialClusterInds(a,:)==i) = j; % update that node j has replaced node i
                trialClusterInds(b,trialClusterInds(b,:)==j) = i; % update that node i has replaced node j
                ARowiMinusRowj = full( A(i,:) - A(j,:) ); % "full" needed for speed if A is sparse
                %
                clusterIndMat = trialClusterInds(equalSizeInds,:); % concatenate all possible group indices into a matrix
                sumAijc = sum(ARowiMinusRowj(clusterIndMat),2);
                trialACounts(equalSizeInds,a) = trialACounts(equalSizeInds,a) - sumAijc;
                trialACounts(equalSizeInds,b) = trialACounts(equalSizeInds,b) + sumAijc;
                if smallerLastGroup % take care of last group separately, if unequal size
                    sumAijEnd = sum(ARowiMinusRowj(trialClusterInds(k,trialClusterInds(k,:)>0)));
                    trialACounts(k,a) = trialACounts(k,a) - sumAijEnd;
                    trialACounts(k,b) = trialACounts(k,b) + sumAijEnd;
                end
                %
                % Update the above for special cases (c==a) || (c==b)
                trialACounts(a,a) = trialACounts(a,a) + A(i,j);
                trialACounts(b,b) = trialACounts(b,b) + A(i,j);
                if smallerLastGroup && (b==k)
                    trialACounts(a,b) = trialACounts(a,b) - sum(ARowiMinusRowj(trialClusterInds(b,trialClusterInds(b,:)>0))) - 2*A(i,j);
                else
                    trialACounts(a,b) = trialACounts(a,b) - sum(ARowiMinusRowj(trialClusterInds(b,:))) - 2*A(i,j);
                end
                trialACounts(b,a) = trialACounts(a,b);
                %
                % Normalize and respect symmetry of trialAbar matrix
                trialACounts(a,:) = trialACounts(:,a)';
                trialACounts(b,:) = trialACounts(:,b)';
                
                % Now calculate changed likelihood directly
                thetaCola = trialACounts(:,a)./habSqrdCola;
                thetaColb = trialACounts(:,b)./habSqrdColb;
                thetaEntryab = trialACounts(a,b)./habSqrdEntryab;
                %
                thetaCola(thetaCola<=0) = eps;     % else 0*log(0) evaluates to NaN
                thetaCola(thetaCola>=1) = 1 - eps; % else (1-1)*log(1-1) evaluates to NaN
                %
                thetaColb(thetaColb<=0) = eps;     % else 0*log(0) evaluates to NaN
                thetaColb(thetaColb>=1) = 1 - eps; % else (1-1)*log(1-1) evaluates to NaN
                %
                thetaEntryab(thetaEntryab<=0) = eps;     % else 0*log(0) evaluates to NaN
                thetaEntryab(thetaEntryab>=1) = 1 - eps; % else (1-1)*log(1-1) evaluates to NaN
                %
                % for this to work, we will have had to subtract out terms prior to updating
                deltaNegEnt = sum( habSqrdCola.*( thetaCola.*log(thetaCola) + (1-thetaCola).*log(1-thetaCola) ) ...
                    + habSqrdColb.*( thetaColb.*log(thetaColb) + (1-thetaColb).*log(1-thetaColb) ) )...
                    - habSqrdEntryab*( thetaEntryab.*log(thetaEntryab) + (1-thetaEntryab).*log(1-thetaEntryab) );
                %
                oldDeltaNegEnt = sum( habSqrdCola.*( oldThetaCola.*log(oldThetaCola) + (1-oldThetaCola).*log(1-oldThetaCola) ) ...
                    + habSqrdColb.*( oldThetaColb.*log(oldThetaColb) + (1-oldThetaColb).*log(1-oldThetaColb) ) ) ...
                    - habSqrdEntryab*( oldThetaEntryab.*log(oldThetaEntryab) + (1-oldThetaEntryab).*log(1-oldThetaEntryab) );
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Update log-likelihood - O(k)
                trialLL = trialLL + (deltaNegEnt-oldDeltaNegEnt)/sampleSize;
                
            end % if a~=b
            
        end % for mmm = 1:geomVec(m)
        
        % Metroplis or greedy step; if trial clustering accepted, then update current <-- trial
        if trialLL > currentLL
            currentLabelVec = trialLabelVec;
            currentLL = trialLL;
            currentACounts = trialACounts;
            currentClusterInds = trialClusterInds;
            
        end
        
    end
        
    % Keep track of best clustering overall
    if currentLL > bestLL % replace and save if trialLL is an improvement
        bestLL = currentLL; % update globally best visited likelihood
        bestLabelVec = currentLabelVec;
        bestCount = bestCount + 1;
        
    end % if currentLL > bestLL
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~mod(mm,5)
        tElapsedOuter = toc(tStartOuter);
        normalizedBestLL = bestLL*2*sampleSize/sum(A(:));
        display([num2str(normalizedBestLL) ' LL.  Iter ' int2str(mm) ' of max ' int2str(maxNumRestarts) '; ' int2str(bestCount) ' global improvements; took ' num2str(tElapsedOuter) ' s']);
        tStartOuter = tic;
        if bestCount == 0
            consecZeroImprovements = consecZeroImprovements + 1;
        else
            bestCount = 0;
            consecZeroImprovements = 0;
        end
        if normalizedBestLL - oldNormalizedBestLL < absTol,
            tolCounter = tolCounter + 1;
        else
            tolCounter = 0;
        end
        oldNormalizedBestLL = normalizedBestLL;
        if tolCounter >= 3
            display('3 consecutive likelihood improvements less than specified tolerance; quitting now');
            break
        end
        if allInds == 1
            if consecZeroImprovements == 2
                display('Local optimum likely reached in random-ordered greedy likelihood search; quitting now');
                break
            end
        else
            if consecZeroImprovements == ceil(k*nchoosek(n,2)/numMHsteps)
                display('Local optimum likely reached in random-ordered greedy likelihood search; quitting now');
                break
            end
        end
    end
    
end

% END MAIN FUNCTION GRAPHEST_FASTGREEDY.M

function Xsums = getSampleCounts(X,clusterInds)

numClusters = size(clusterInds,1);
Xsums = zeros(numClusters,numClusters); % initialize
for b = 1:numClusters % sum over strict upper-triangular elements
    for a = 1:(b-1)
        clusterIndsa = clusterInds(a,clusterInds(a,:)>0);
        clusterIndsb = clusterInds(b,clusterInds(b,:)>0);
        Xsums(a,b) = sum(sum(X(clusterIndsa,clusterIndsb),1));
    end
end
Xsums = Xsums + Xsums'; % diagonal still 0
for a = 1:numClusters % relies on A begin symmetric and with no self-loops
    clusterIndsa = clusterInds(a,clusterInds(a,:)>0);
    Xsums(a,a) = sum(sum(X(clusterIndsa,clusterIndsa),1)) ./ 2;
end

function normLogLik = fastNormalizedBMLogLik(thetaVec,habSqrdVec,sampleSize)

thetaVec(thetaVec<=0) = eps;     % else 0*log(0) evaluates to NaN
thetaVec(thetaVec>=1) = 1 - eps; % else (1-1)*log(1-1) evaluates to NaN
negEntVec = thetaVec.*log(thetaVec) + (1-thetaVec).*log(1-thetaVec);
normLogLik = sum(habSqrdVec.*negEntVec) / sampleSize;
