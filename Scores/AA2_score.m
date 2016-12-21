function [ output_activity ] = AA2_score( AAsequence, Nsequence, skip )
% Dec-09-2016, janghyun.yoo@nih.gov
%
% Get amber codon termination prob. Combine Mao_score and the last 2 AAs
% based on EMBO, 13, 249 (1994)
% 
% Usage
% termination_percentage=AA2_score(AAsequence, Nsequence)
% termination_percentage=AA2_score(AAsequence, Nsequence, correction)
%
% AAsequence: the last two aminoacids (single character) before TAG.
% Sequence: 9 neucleotide sequence NNNTAGNNN
% correction: 0 or 1. If it is 1, correct the score of +4 position 
% based on Mottagui-tabar, Björnsson, Isaksson , EMBO, 14, 152 (1995)

if nargin<3
    skip=0;
end
AA2_multiplier={
'A',	1.445325405;
'R',	0.506802415;
'N',	0.871976041;
'D',	1.998383595;
'C',	1.268346784;
'Q',	1.138013691;
'E',    1.832466138;
'G',	1.64145712;
'H',	0.836551399;
'I',	1.238850347;
'L',	1.268346784;
'K',	1.29771285;
'M',	1.445325405;
'F',    1.788401339;
'P',	1.496623556;
'S',	1.235745459;
'T',	0.655861948;
'W',	1.351473106;
'Y',	1.255705454;
'V',	1.486620416
};

for ii=1:size(AA2_multiplier,1)
    if strcmpi(AA2_multiplier{ii,1},AAsequence(1))
        if fix(skip/2)==1
            output_activity=1.228429299*AA1_score(AAsequence(2),Nsequence, skip); % skip this AA, use average value
        else
            output_activity=AA2_multiplier{ii,2}*AA1_score(AAsequence(2),Nsequence, skip);
        end
        return;
    end
end
h=errordlg(['AA2, unknown AA: ' AAsequence(1)]);
