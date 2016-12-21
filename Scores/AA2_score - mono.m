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
'A',	1.362609982;
'R',	0.477798305;
'N',	0.822073184;
'D',	1.884016863;
'C',	1.19575978;
'Q',	1.072885601;
'E',    1.727594799;
'G',	1.547517154;
'H',	0.788675881;
'I',	1.167951413;
'L',	1.19575978;
'K',	1.223445237;
'M',	1.362609982;
'F',    1.686051811;
'P',	1.410972359;
'S',	1.165024217;
'T',	0.618327219;
'W',	1.274128814;
'Y',	1.183841909;
'V',	1.401541696
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
