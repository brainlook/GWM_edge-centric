function mi = calmi(u1, u2, wind_size)
%   mi = calmi(u1, u2, wind_size)
%
%   Standardized mutual information between two clustering results
% 
%   Inputs:
%       u1, 
%           clustering result1 
%       u2, 
%           clustering result2
%       wind_size
%           size of window
%   Outputs:
%       mi,
%           standardized mutual information
%           range [0,1]
% 

x = [u1, u2];
n = wind_size;
[xrow, xcol] = size(x);
bin = zeros(xrow,xcol);
pmf = zeros(n, 2);
for i = 1:2
    minx = min(x(:,i));
    maxx = max(x(:,i));
    binwidth = (maxx - minx) / n;
    edges = minx + binwidth*(0:n);
    histcEdges = [-Inf edges(2:end-1) Inf];
    [occur,bin(:,i)] = histc(x(:,i),histcEdges,1); %The histogram distribution of a single vector is computed by means of histogram
    pmf(:,i) = occur(1:n)./xrow;
end
%Calculate the joint probability density of u1 and u2
jointOccur = accumarray(bin,1,[n,n]);  %（xi，yi）The number of pieces of data that fall into the n*n square at the same time is the joint probability density
jointPmf = jointOccur./xrow;
Hx = -(pmf(:,1))'*log2(pmf(:,1)+eps);
Hy = -(pmf(:,2))'*log2(pmf(:,2)+eps);
Hxy = -(jointPmf(:))'*log2(jointPmf(:)+eps);
MI = Hx+Hy-Hxy;
mi = MI/sqrt(Hx*Hy);	
