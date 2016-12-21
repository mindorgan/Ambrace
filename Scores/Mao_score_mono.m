function [ output_activity ] = Mao_score( sequence)
% Dec-09-2016, janghyun.yoo@nih.gov
%
% Get amber codon termination prob. based on the mononeucleotide model
% that was decribed in Computational biology and chemistry 28 (2004)
% by P.-L. Mao.
% 
% Usage
% termination_percentage=Mao_score(sequence)
% termination_percentage=Mao_score(sequence, correction)
%
% Sequence: 9 neucleotide sequence NNNTAGNNN

for ii=1:length(sequence)
    if strcmpi(sequence(ii),'U')
        sequence(ii)='T';
    end
end
if length(sequence) ~= 9
    h = errordlg(['Mao, the sequence should be 9 bp length: ' num2str(length(sequence))]);
    return;
end

if ~strcmpi(sequence(4:6),'TAG')
    h = errordlg(['Mao, Stop codon is not Amber: ' sequence(4:6)]);
    return;
end

    function [vector]=p_decompose(neucleotide)
        if neucleotide=='T' || neucleotide=='t'
            vector=[-3.96;1.61;-1.72;0.53];
        elseif neucleotide=='G' || neucleotide=='g'
            vector=[-0.31;-3.76;0.03;1.41];
        elseif neucleotide=='A' || neucleotide=='a'
            vector=[2.25;0.77;-1.17;0.75];
        elseif neucleotide=='C' || neucleotide=='c'
            vector=[-2.36;0.72;-2.07;-1.76];
        else
            h = errordlg(['Mao, Unknown neucleotide: ' neucleotide]);
            return;
        end
    end

P=zeros(6,4); % Score coefficient for each principal components
P(1,:)=[-0.0373 0.0028 -0.0225 -0.0186]; % for -3 position
P(2,:)=[-0.0153 -0.0016 -0.0522 -0.0526]; % for -2 position
P(3,:)=[0.0961 0.0945 -0.131 -0.0054]; % for -1 position
P(4,:)=[0.2049 -0.0722 0.2106 0.0578]; % for +4
P(5,:)=[-0.1603 0.1712 -0.1132 0.2899];
P(6,:)=[-0.0123 0.0282 0.0042 0.079];

sequence_target=[sequence(1:3) sequence(7:9)];
output_activity=1.7520;
for ii=1:length(sequence_target)
    output_activity=output_activity+P(ii,:)*p_decompose(sequence_target(ii));    
end

%output_activity=-9.43*exp(-0.015*output_activity)+9.51857;
output_activity=0.0675151*exp(0.433047*output_activity)+0.01*1.6814; % 1.6814% is average basal activity
end