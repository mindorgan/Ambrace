function [ output_activity ] = AA1_score( AA1, Nsequence, skip )
% Dec-09-2016, janghyun.yoo@nih.gov
%
% Get amber codon termination prob. Combine Mao_score and the last AA score
% based on Zhang, Rydén-Aulin, Isaksson JMB 261, 98 (1996).
% 
% Usage
% termination_percentage=AA1_score(AA1, Nsequence)
% termination_percentage=AA1_score(AA1, Nsequence, correction)
%
% AA1: the last aminoacid (single character) before TAG.
% Sequence: 9 neucleotide sequence NNNTAGNNN

if nargin < 3
    skip=0;
end

AA1_multiplier={
'A',	1.450106682;
'R',	1.720063992;
'N',	2.033714466;
'D',	1.033651878;
'C',	1.794885956;
'Q',	1.284976468;
'E',    1.153846154;
'G',	2.553121546;
'H',	1.392126392;
'I',	1.354643033;
'L',	2.021795508;
'K',	1.271700155;
'M',	1.465565864;
'F',    1.80496674;
'P',	1.800511249;
'S',	1.526060519;
'T',	1.611815239;
'W',	1.921878003;
'Y',	1.360940993;
'V',	0.703150113
};

for ii=1:size(AA1_multiplier,1)
    if strcmpi(AA1_multiplier{ii,1},AA1)
        if rem(skip,2)==0
            output_activity=1.562976048*Mao_score(Nsequence); % 1.5629... is average of AAs.
        else
            output_activity=AA1_multiplier{ii,2}*Mao_score(Nsequence);
        end
        return;
    end
end
h=errordlg(['AA1, unknown AA: ' AA1]); % no AA matching to AA1