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
    if strcmpi(sequence(ii),'T')
        sequence(ii)='U';
    end
end
if length(sequence) ~= 9
    h = errordlg(['Mao, the sequence should be 9 bp length: ' num2str(length(sequence))]);
    return;
end

if ~strcmpi(sequence(4:6),'UAG')
    h = errordlg(['Mao, Stop codon is not Amber: ' sequence(4:6)]);
    return;
end

    dn1=sequence(3:4);
    dn2=sequence(6:7);
    dn3=sequence(7:8);

    di_pp_list={'AA' 'AG' 'AC' 'AU' 'GA' 'GG' 'GC' 'GU' 'CA' 'CG' 'CC' 'CU' 'UA' 'UG' 'UC' 'UU'};
    pp1_values=[1.030 0.746 0.819 0.370 0.748 0.762 0.770 0.812 0.905 0.922 0.950 -0.050 0.983 0.945 0.996 0.954];
    pp2_values=[-0.585 -0.192 0.848 2.158 -0.813 -0.398 0.609 0.651 -1.653 -1.394 -0.412 -0.775 -1.706 -1.343 -0.319 -0.448];
    pp3_values=[0.640 0.733 0.911 0.868 0.396 0.569 0.761 0.864 0.211 0.310 0.463 -0.611 0.105 0.243 0.406 0.430];

    output_activity=-9.5415;
    
    %position (-1,+1)
    for ii=1:numel(di_pp_list)
        if strcmpi(di_pp_list{ii},dn1)
            output_activity=output_activity+[0.09999 0.1408 0.2263]*[pp1_values(ii); pp2_values(ii); pp3_values(ii)];
            break;
        end
    end
    %position (+3,+4)
    for ii=1:numel(di_pp_list)
        if strcmpi(di_pp_list{ii},dn2)
            output_activity=output_activity+[17.3399 1.2587 1.2158]*[pp1_values(ii); pp2_values(ii); pp3_values(ii)];
            break;
        end
    end
    
    %position (+4,+5)
    for ii=1:numel(di_pp_list)
        if strcmpi(di_pp_list{ii},dn3)
            output_activity=output_activity+[1.2935 0.0182 0.3655]*[pp1_values(ii); pp2_values(ii); pp3_values(ii)];
            break;
        end
    end
    
output_activity=0.620164*exp(-0.217512*output_activity)-0.0950818;
    
    pp_coeff=[-0.1046;
0.1299;
-0.1711;
0.0020; 0.0830;
-0.2473; -0.2154;
0.2038;
-0.2252; 14.4488;
-1.1574; -1.1446; -1.54;
0.0471;
-0.2765; 0.0094;
-0.042; 0.083];
output_activity=-7.07413;

%for jj=1:numel(sequence)-3
for jj=3:5
    if jj<4
        ss=sequence(jj:jj+1);
    else
        ss=sequence(jj+2:jj+3);
    end
    for ii=1:numel(di_pp_list)
        if strcmpi(di_pp_list{ii},ss)
            coeff_vector=pp_coeff(jj*3-2:jj*3);
            coeff_vector=coeff_vector';
            output_activity=output_activity+coeff_vector*[pp1_values(ii); pp2_values(ii); pp3_values(ii)];
            break;
        end
    end
end
output_activity=0.0184658*exp(0.759883*output_activity);
end