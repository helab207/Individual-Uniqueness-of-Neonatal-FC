function [mask_FC1, mask_FC2, mask] = strong_mask(FC1,FC2)
%==========================================================================
% This code is used to generate a strong edge mask that can be used to 
% extract strong edges in individual FC matrix.
% Note: each strong edge was determined by a one-sample t-test across all
% subjects and multiple comparison correction.
%
% Syntax:  [mask_FC1, mask_FC2, mask] = strong_mask(FC1,FC2)
%
% Input:
%        FC1,FC2:
%               1*N cell, FC matrices of N subjects in Section 1 
%               and Section 2, respectively.
%
% Outputs:
%     mask_FC1, mask_FC2
%        The strong edge mask of Section1 and Section 2, respectively.
%     mask:
%        A strong edge mask derived by mask_FC1 intersecting mask_FC2.
%
% Qiushi Wang, BNU, BeiJing, 2021/2/23, aqiubshic@126.com
%==========================================================================

n_sub = length(FC1);
n_edge = length( squareform(FC1{1}) );

Edge1 = zeros(n_edge, n_sub);
Edge2 = zeros(n_edge, n_sub);
for i = 1 : n_sub
    Edge1(:,i) = squareform( FC1{i} );
    Edge2(:,i) = squareform( FC2{i} );
end

% For Section 1
p_array = zeros(n_edge,1);
for i = 1 : n_edge
    tmp = Edge1(i,:);
    [~, p_array(i)] = ttest(tmp);
end
q = 0.05;
[fdr1,~] = gretna_FDR(p_array,q);
p_array( p_array >= fdr1) = 0;
p_array( p_array ~=0 ) =1;
mask_FC1 = squareform(p_array);

% For Section 2
p_array = zeros(n_edge,1);
for i = 1 : n_edge
    tmp = Edge2(i,:);
    [~, p_array(i)] = ttest(tmp);
end
q = 0.05;
[fdr2,~] = gretna_FDR(p_array,q);
p_array( p_array >= fdr2) = 0;
p_array( p_array ~=0 ) =1;
mask_FC2 = squareform(p_array);

mask = mask_FC1 .* mask_FC2;
