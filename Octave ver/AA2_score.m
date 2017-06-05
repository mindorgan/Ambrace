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
'A',	1.251519702;
'R',	0.438844571;
'N',	0.755051556;
'D',	1.730417548;
'C',	1.098272392;
'Q',	0.985415845;
'E',    1.586748194;
'G',	1.421351842;
'H',	0.724377053;
'I',	1.072731174;
'L',	1.098272392;
'K',	1.123700721;
'M',	1.251519702;
'F',    1.548592106;
'P',	1.295939213;
'S',	1.070042624;
'T',	0.567916504;
'W',	1.170252189;
'Y',	1.087326156;
'V',	1.287277408
};
for ii=1:size(AA2_multiplier,1)
    if strcmpi(AA2_multiplier{ii,1},AAsequence(1))
        if fix(skip/2)==1
            output_activity=AA2_multiplier{ii,2}*AA1_score(AAsequence(2),Nsequence, skip);
        else
            output_activity=mean(cell2mat(AA2_multiplier(:,2)))*AA1_score(AAsequence(2),Nsequence, skip); % skip this AA, use average value
        end
        return;
    end
end
h=errordlg(['AA2, unknown AA: ' AAsequence(1)]);
