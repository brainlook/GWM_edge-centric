function triplet = fcn_GWMtriplet(ci)
%   The GWM triplets are calculated from the contiguous communities of the clusters
%   closed loop--The three connecting edges belong to the same communities
%   forked loop-- Both sides of the WM reference brain region belong to the same community,
                % and the other side belongs to a different community
%   L_shape loop--One of the sides connecting the WM reference areas belongs to the same community 
                % as the side connecting the two GM areas, and the other side belongs to a different community
%   diverse loop-- The three sides belong to three different communities
%   Inputs:
%       ci,
%          	clustering result 
%   Outputs:
%       triplet,
%          	GWM triplet,size: (number of WM regions)x4(loop)
%    

n_gm = 200;
n_wm = 128;
n = n_gm + n_wm;
M = n_gm*(n_gm-1)/2;  % Number of GM edges
[u,v] = find(triu(ones(n),1));  % get edges

ci_index = zeros(n,n);
kk=find(triu(ones(n),1));
ci_index(kk)=ci; % locate  the value of the upper triangular matrix

% A certain WM brain region as the reference node
triplet_model = zeros(n_wm,4);
ssum = 0;

tic
for ind=n_gm+1:n
    gm_index = find(v==ind);   % There are 200 GM brain regions connected to it,
                               % 200 gray matter brain regions arranged in order
    gm_index = gm_index(1:n_gm);
    for i=1:n_gm-1
        for j=i+1:n_gm

            wm2gm_1 = ci(gm_index(i));
            wm2gm_2 = ci(gm_index(j));
            gm2gm = ci_index(i,j);

            % Determine the triplet type
            if length(unique([gm2gm,wm2gm_1,wm2gm_2]))==1  % four types of triplets
                triplet_model(ind-n_gm,1)=triplet_model(ind-n_gm,1)+1;
            elseif wm2gm_1==wm2gm_2
                triplet_model(ind-n_gm,2)=triplet_model(ind-n_gm,2)+1;
            elseif wm2gm_1==gm2gm || wm2gm_2==gm2gm    
                triplet_model(ind-n_gm,3)=triplet_model(ind-n_gm,3)+1;
            else 
                triplet_model(ind-n_gm,4)=triplet_model(ind-n_gm,4)+1;
            end

            ssum = ssum+1;
        end
    end

end
toc
triplet=triplet_model./M;


end
