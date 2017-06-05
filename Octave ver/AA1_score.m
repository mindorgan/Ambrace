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
    skip=0; %skip the AA1 scoring. Set the average value
end

AA1_multiplier={
'A',	2.141927627;
'R',	2.507009011;
'N',	3.045582214;
'D',	1.271922436;
'C',	2.702253263;
'Q',	1.329914548;
'E',    1.853568119;
'G',	3.061428297;
'H',	1.677073624;
'I',	2.07038582;
'L',	2.348628477;
'K',	1.567719993;
'M',	1.720667828;
'F',    3.016391513;
'P',	2.872761363;
'S',	2.039614316;
'T',	2.442721375;
'W',	2.411112287;
'Y',	1.650281727;
'V',	1.512793914
};

for ii=1:size(AA1_multiplier,1)
    if strcmpi(AA1_multiplier{ii,1},AA1)
        if rem(skip,2)==1
            output_activity=AA1_multiplier{ii,2}*Mao_score(Nsequence);
        else
            output_activity=mean(cell2mat(AA1_multiplier(:,2)))*Mao_score(Nsequence); %mean of AAs.
        end
        return;
    end
    
end
h=errordlg(['AA1, unknown AA: ' AA1]); % no AA matching to AA1