function ci = fcn_edgecomm(ts,k,nMax)
%   ci = fcn_edgecomm(ts,k,nMax)
%
%   Repeat the clustering 250 times and select the 
%   partition that is most similar to the other partitions
%
%   Inputs:
%       ts, 
%           time series, size: (time)x(node) 
%       k, 
%           number of edge communities
%       nMax
%           epochs
%   Outputs:
%       ci,
%          	clustering result
%    

ets = fcn_edgets(ts);
N = size(ts,2);
M = N*(N-1)/2;
ci_all =zeros(M,nMax);
for epoch = 1: nMax
    tic
    disp(['*****',num2str(epoch),'******']);
    ci_all(:,epoch) = kmeans(ets',k,...
        'distance','sqeuclidean',...
        'maxiter',1000);  
    toc
end

% calculate the standardized mutual information
s=size(ci_all);
mi=zeros(s(2));
for i=1:s(1,2)-1
    for j=i+1:s(1,2)
        mi(i,j) = calmi(ci_all(:,i),ci_all(:,j),N);
    end
end

mi = mi+mi';
mi_mean = sum(mi)/(nMax-1);
[~,mi_max_index] = max(mi_mean);
ci = ci_all(:,mi_max_index);
end

