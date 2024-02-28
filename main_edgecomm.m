% main step of generating edge community and related indexes
addpath(genpath('function'));
%% generate edge community
load ./Data/sample.mat
k = 10; % number of edge community
nMax = 250; % number of clustering epochs
n_gm = 200; % number of GM regions
N = 328; % number of all regions
ci_gwm = fcn_edgecomm(ts,k,nMax);
ci_gm = fcn_edgecomm(ts(:,1:n_gm),k,nMax);

% GWM/GM edge community
mat = zeros(N); mat(triu(ones(N),1) > 0) = ci_gwm; mat_gwm = mat + mat';
mat = zeros(n_gm); mat(triu(ones(n_gm),1) > 0) = ci_gm; mat_gm = mat + mat';

%% calculate the GM-WM triplet
triplet = fcn_GWMtriplet(ci_gwm);
load ./Data/label.mat
wm_lab = lab(n_gm+1:N)-16;   % WM node serial number
triplet_WMnets = grpstats(triplet, wm_lab);  % Find the average value of WM network

%% evulation indexes and delta
[u,v] = find(triu(ones(N),1));
 [~,enorm1] = fcn_node_entropy(ci_gwm,u,v,N);
[u,v] = find(triu(ones(n_gm),1));
 [~,enorm2] = fcn_node_entropy(ci_gm,u,v,n_gm);
enorm_delta = enorm1(1:n_gm)-enorm2;

s_gwm = fcn_profilesim(mat_gwm);
s_gm = fcn_profilesim(mat_gm);
s_delta = s_gwm(1:n_gm,1:n_gm)-s_gm;

enorm_net= cell(8,1);
s_net = cell(8,1);
for ind = 1:8   % number of GM networks
    index_tmp = find(lab==ind);
    s_net{ind} = reshape(s_delta(index_tmp,index_tmp),[],1);
    enorm_net{ind} = enorm_delta(index_tmp);
end


