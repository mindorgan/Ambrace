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
'A',	1.81692113;
'R',	1.967057842;
'N',	2.288496088;
'D',	1.238941244;
'C',	2.122647537;
'Q',	1.499777501;
'E',    1.917071227;
'G',	3.121580018;
'H',	1.634406598;
'I',	1.936040048;
'L',	2.168404875;
'K',	1.41980403;
'M',	1.539502609;
'F',    2.415600388;
'P',	2.193923775;
'S',	1.736121872;
'T',	2.057517085;
'W',	2.144225255;
'Y',	1.608616082;
'V',	1.357269819
};

for ii=1:size(AA1_multiplier,1)
    if strcmpi(AA1_multiplier{ii,1},AA1)
        if rem(skip,2)==1
            output_activity=1.909196251*Mao_score(Nsequence); % 1.5629... is average of AAs.
        else
            output_activity=AA1_multiplier{ii,2}*Mao_score(Nsequence);
        end
        return;
    end
end
h=errordlg(['AA1, unknown AA: ' AA1]); % no AA matching to AA1