function [V_edge, V_ROI] = variability( FC1, FC2)
%==========================================================================
% This code calculates the inter-subject variability of functional 
% connection strength (FCS) of the FC matrix at the edge and ROI level. 
% Note: the variability is calculated as the standard deviation of FCS 
% across all subjects 
%
% Syntax:  [V_edge, V_ROI] = variability( FC1, FC2)
%
% Input:
%        FC1,FC2:
%               1*N cell, FC matrices of N subjects in Section 1 
%               and Section 2, respectively.
%
% Outputs:
%     V_edge:
%        The FCS variability of each edge.
%     V_ROI
%        The FCS variability of each ROI.
%
% Qiushi Wang, BNU, BeiJing, 2021/2/23, aqiubshic@126.com
%==========================================================================

n_sub = length(FC1);
n_edge = length( squareform(FC1{1}) );

Edge1 = zeros(n_sub, n_edge);
Edge2 = zeros(n_sub, n_edge);
for i = 1 : n_sub
    Edge1(i,:) = squareform( FC1{i} );
    Edge2(i,:) = squareform( FC2{i} );
end

V_edge = (( std(Edge1) + std(Edge2) ) / 2);
tmp = squareform(V_edge);
V_ROI = (sum(tmp) / (length(tmp)-1));
