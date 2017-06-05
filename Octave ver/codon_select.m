function codon_table=codon_select(codon_file)

fileID=fopen(codon_file);
tline = fgetl(fileID);

assign_begin=0; % you can assign codons to amino acids

codon_table={};
while ischar(tline)
    tline=strtrim(tline); % remove white spaces
    
    if numel(tline)<1
        tline = fgetl(fileID);
        continue;
    end
    
    if strcmp(tline(1),'/')
        tline = fgetl(fileID);
        continue;        
    end
    
    %Assign codons
    if isstrprop(tline(1),'digit')
        for ii=2:numel(tline)
            AA=tline(ii);
            if isstrprop(AA,'alpha')
                assign_begin=1;
                codon_table.(AA)={};
                break;
            end
        end
    elseif isstrprop(tline(1),'alpha') && numel(tline)==3
        if assign_begin==1 % you know where to add these codons
            codon_table.(AA){end+1}=tline;
        end
    end
    tline = fgetl(fileID);
end

fclose(fileID);