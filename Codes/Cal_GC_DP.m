function [GCmatrix, DPmatrix] = Cal_GC_DP( FC1, FC2 )
%==========================================================================
% This code calculates the group consistency (i.e. GC) and differential 
% power (i.e. DP) values for each edge in the matrix based on the 
% connectivity matrices of all subjects in two sections.
%
% Syntax:  function [GCmatrix, DPmatrix] = GC_DP( FC1, FC2 )
%
% Input:
%        FC1,FC2:
%               1*N cell, FC matrices of N subjects in Section 1 
%               and Section 2, respectively.
%
% Outputs:
%     GCmatrix:
%        The matrix that stores the group consistency value for each edge.
%     DPmatrix:
%        The matrix that stores the differential power value for each edge.
%
% Reference: Finn (2015) -- Nat neurosci. 18(11):1664-1671.
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
Edge1=zscore(Edge1); Edge2=zscore(Edge2);

R_ii = Edge1 .* Edge2;    % R_ii:the product of all corresponding edges between identical subjects in two sections

% Calculation of GC value
GCedge_eachsubject = R_ii;
GCedge = mean( R_ii, 2);
GCmatrix = squareform(GCedge);

% Calculation of DP value
R_ij = cell(1,40); R_ji = cell(1,40);
for i=1:n_sub        % R_ij & R_ji:the product of all corresponding edges between different subjects in two sections
    for j=1:n_sub
        R_ij{i}(:,j) = Edge1(:,i) .* Edge2(:,j);
        R_ji{i}(:,j) = Edge2(:,i) .* Edge1(:,j);
    end
    R_ij{i}(:,i)=[];
    R_ji{i}(:,i)=[];
end

count_ij = cell(1,40); count_ji = cell(1,40);
for i=1:n_sub    % for each edge in each subject, count the times for "R_ii(:,i) < R_ij" and "R_ii(:,i) < R_ji".
    for j=1:n_sub-1
        count_ij{i}(:,j) = R_ii(:,i) < R_ij{i}(:,j);
        count_ji{i}(:,j) = R_ii(:,i) < R_ji{i}(:,j);
    end
end

for i=1:n_sub
    P(:,i) = (sum(count_ij{i},2)+sum(count_ji{i},2)) ./ (2*(n_sub-1));
end

P(P==0)=10^-7;  % since log(0) doesn't make sense, we set the zeroes in 'P' as a minimum of 10^-7. 
                % ref: DP_edge_example.m from Horien et al.              
DPedge_eachsubject = -log(P);
DPedge = sum(DPedge_eachsubject,2);
DPmatrix = squareform(DPedge);
