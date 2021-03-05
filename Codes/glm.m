function [Pvalue, FDRvalue, Effect] = glm( data, age, sex, interaction, mFD )
%==========================================================================
% This code calculates which edges in the data have age, sex, or 
% interaction effects based on a general linear model, using mFD 
% (i.e. mean framewise displacement) as a covariate.
%
% Syntax: [Pvalue, FDRvalue, Effect] = glm( data, age, sex, interaction, mFD )
%
% Input:
%        data:
%              The dependent variable in the general linear model, the
%              values of all the edges of all the subjects are stored in 
%              an M*N marix, with M denotes all the edges of a subject and 
%              N denotes the number of subjects.
%        age:
%              A vector that stores the age of all subjects.
%        sex:
%              A vector that stores the sex of all subjects.
%              Use different numbers for male and female (e.g. 1 and 0).
%        interaction:
%              A vector that stores the interaction between age and sex of 
%              all subjects, obtained by age*sex.
%        mFD:
%              A vector that stores the head motion parameters of all 
%              subjects. For each subject, the mFD value is calculated as 
%              the average FD value of section 1 and section 2.
%
% Outputs:
%     Pvalue:
%          The p value of each edge under each predictor (i.e. age, sex, 
%          and interaction).
%     FDRvalue:
%          The FDR value of each predictor (i.e. age, sex, and interaction).
%     Effect:
%          It stores which edges show significant effect under each 
%          predictor (i.e. age, sex, and interaction) after multiple 
%          comparison correction, with 1 indicating significant effect and 
%          0 indicating no significant effect.
%
% Qiushi Wang, BNU, BeiJing, 2021/2/23, aqiubshic@126.com
%==========================================================================

n_edge = length(data(:,1));

age_p = zeros( n_edge,1); sex_p = zeros( n_edge,1); interaction_p = zeros( n_edge,1);
for i = 1 : n_edge
    edgeX = data(i,:)';
    tstats_linear = regstats(edgeX,[age sex interaction mFD],'linear','tstat');
    age_p(i) = tstats_linear.tstat.pval(2);
    sex_p(i) = tstats_linear.tstat.pval(3);
    interaction_p(i) = tstats_linear.tstat.pval(4);
end
q = 0.05;
[age_fdr,~] = gretna_FDR(age_p, q);
[sex_fdr,~] = gretna_FDR(sex_p, q);
[interaction_fdr,~] = gretna_FDR(interaction_p, q);
%
if ~isempty(age_fdr)
    age_effect = age_p; age_effect(age_effect >= age_fdr) = 0; age_effect(age_effect ~=0) = 1;
else
    age_effect = [];
end
if ~isempty(sex_fdr)
    sex_effect = sex_p; sex_effect(sex_effect >= sex_fdr) = 0; sex_effect(sex_effect ~=0) = 1;
else
    sex_effect = [];
end
if ~isempty(interaction_fdr)
    interaction_effect = interaction_p; interaction_effect(interaction_effect >=interaction_fdr) = 0;
    interaction_effect(interaction_effect ~=0) = 1;
else
    interaction_effect=[];
end
% exclude the edges with interaction effect from the edges with age effect
% as the final edges with ege effect
if ~isempty(age_effect)
ind =  (age_effect + interaction_effect)==2 ;
age_effect(ind) = 0;
end
% exclude the edges with interaction effect from the edges with sex effect
% as the final edges with sex effect
if ~isempty(sex_effect)
ind =  (sex_effect + interaction_effect)==2 ;
sex_effect(ind) = 0;
end

Pvalue{1} = age_p; Pvalue{2} = sex_p; Pvalue{3} = interaction_p;
FDRvalue{1} = age_fdr; FDRvalue{2} = sex_fdr; FDRvalue{3} = interaction_fdr;
Effect{1} = age_effect; Effect{2} = sex_effect; Effect{3} = interaction_effect;
