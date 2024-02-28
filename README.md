## code for generating WM functional regions and GWM edge-centric network
Repository for the paper:
Title: fMRI signals in white matter rewire gray matter community organization
Author: Huanyu Xu, Wenjing Hu, Huanxin Wang, Yiwen Gao, Luyao Wang, Jiehui Jiang

There are two main codes:
region_grow.m (use random growth to divide WM functional networks to regions)
main_edgecomm.m (generate the stable edge community and related evalution indexes)

for main_edgecomm.m there are some subfunctions  
calmi.m -- calculate the mutual information and find the stable edge community  
fcn_edgets.m -- calculate the edge time series  
fcn_node_entropy.m --  calculate the community overlap entropy  
fcn_profilesim.m -- calculate the edge community similarity  
fcn_edgecomm.m -- multiple epochs for select the stable edge community  
fcn_GWMtriplet.m -- calculate the GM-WM triplet 

where fcn_edgets.m/fcn_node_entropy.m/fcn_edgecomm.m/are from https://github.com/brain-networks/edge-centric_demo
