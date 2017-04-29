% fitting script
clear maoSeq maoScore
load('maoInput.mat');
load('maoActivity.mat');
for i=1:size(maoInput,1)
    maoSeq{i,1}=[maoInput{i,1} maoInput{i,2} maoInput{i,3}];
    maoScore(i,1)=Mao_score(maoSeq{i,1});
end

fitX=maoScore;
fitY=maoActivity;

fitResult=gnuplotfit('a*exp(b*x)',{'a','b'},{'a=1.16','b=0.99'},fitX,fitY);

%% AA1 correction
clear maoScore
load('aa1seq.mat');
for i=1:size(aa1seq,1)
    aa1full{i,1}=[aa1seq{i} 'UAGAGU'];
    maoScore(i,1)=Mao_score(aa1full{i,1});
end