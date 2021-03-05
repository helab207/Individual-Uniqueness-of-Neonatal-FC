function [IDrate12, IDrate21, coef12, coef21] = ID_Predict( FC1, FC2 )
%==========================================================================
% This code identifies individuals between two sets of matrices based on
% the similarity between the matrices. Similarity is measured by Pearson 
% correlation coefficient.
%
% Syntax:  function [IDrate12, IDrate21, coef12, coef21] = ID_Predict( FC1, FC2 )
%
% Input:
%        FC1,FC2:
%               1*N cell, FC matrices of N subjects in Section 1 
%               and Section 2, respectively.
%
% Outputs:
%        IDrate:
%                The identification accuracy from FC1 to FC2 (IDrate12) or 
%                from FC2 to FC1 (IDrate21).
%        coef:
%                The subjects similarity matrix from FC1 to FC2 (coef12) or 
%                from FC2 to FC1 (coef21).
%
% Reference: Finn (2015) -- Nat neurosci. 18(11): 1664-1671.
%
% Qiushi Wang, BNU, BeiJing, 2021/2/23, aqiubshic@126.com
%==========================================================================

n_sub = length( FC1 );

coef12 = zeros(n_sub, n_sub);
coef21 = zeros(n_sub, n_sub);
for i=1:n_sub
    FC1_all_edge = squareform( FC1{i} )';
    for j=1:length(FC2)
        FC2_all_edge = squareform( FC2{j} )';
        coef12(i,j) = corr(FC1_all_edge, FC2_all_edge);
        coef21(j,i) = corr(FC2_all_edge, FC1_all_edge);
    end
end

count12 = 0; count21 = 0;
for i=1:n_sub
    [~,max12] = max(coef12(i,:));
    [~,max21] = max(coef21(i,:));
    if i==max12
        count12 = count12 + 1;
    end
    if i==max21
        count21 = count21 + 1;
    end
end
IDrate12 = sum(count12) / n_sub;
IDrate21 = sum(count21) / n_sub;
