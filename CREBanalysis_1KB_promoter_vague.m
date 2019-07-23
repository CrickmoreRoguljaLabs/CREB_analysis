tic
mypool = parpool;
% 
data_TGACGTCA = zeros(17716,1);

data_AGACGTCA = zeros(17716,1);
data_GGACGTCA = zeros(17716,1);
data_CGACGTCA = zeros(17716,1);

data_TAACGTCA = zeros(17716,1);
data_TTACGTCA = zeros(17716,1);
data_TCACGTCA = zeros(17716,1);

data_TGTCGTCA = zeros(17716,1);
data_TGGCGTCA = zeros(17716,1);
data_TGCCGTCA = zeros(17716,1);
% 
disp('Initiate analysis')
parfor i = 1 : 17716
    if length(seq(i).Sequence) >= 2000
        data_TGACGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'TGACGTCA');

        data_AGACGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'AGACGTCA');
        data_GGACGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'GGACGTCA');
        data_CGACGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'CGACGTCA');

        data_TAACGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'TAACGTCA');
        data_TTACGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'TTACGTCA');
        data_TCACGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'TCACGTCA');

        data_TGTCGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'TGTCGTCA');
        data_TGGCGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'TGGCGTCA');
        data_TGCCGTCA(i) = seqwordcount(seq(i).Sequence(1001:2000),'TGCCGTCA');
    end
end
% 
delete(mypool)
% 

data = [data_TGACGTCA,...
    data_AGACGTCA,data_GGACGTCA,data_CGACGTCA,...
    data_TAACGTCA,data_TTACGTCA,data_TCACGTCA,...
    data_TGTCGTCA,data_TGGCGTCA,data_TGCCGTCA];

toc
%%

genenums = find(data(:,1)>0);

ngenes = length(genenums);

dataoutput = cell(ngenes,11);

for i = 1 : ngenes
    genenum = genenums(i);
    
    tempname = seq(genenum).Header;
    
    initiate = strfind(tempname,'name=') + 5;
    terminate = strfind(tempname,';');
    firstterminate = find(terminate>initiate);
    terminate_actual = terminate(firstterminate(1))-1;
      
    
    dataoutput{i,1} = tempname(initiate:terminate_actual);
    
    dataoutput{i,2} = data(genenum,1);
    dataoutput{i,3} = data(genenum,2);
    dataoutput{i,4} = data(genenum,3);
    dataoutput{i,5} = data(genenum,4);
    dataoutput{i,6} = data(genenum,5);
    dataoutput{i,7} = data(genenum,6);
    dataoutput{i,8} = data(genenum,7);
    dataoutput{i,9} = data(genenum,8);
    dataoutput{i,10} = data(genenum,9);
    dataoutput{i,11} = data(genenum,10);
end

disp('Saving data')
save TGACGTCAgenes_1kb_promoter_vague.mat

%% vague
data2 = [sum(data,2),data];
genenums2 = find(data2(:,1)>0);

ngenes2 = length(genenums2);

dataoutput_vague = cell(ngenes2,12);

for i = 1 : ngenes2
    genenum = genenums2(i);
    
    tempname2 = seq(genenum).Header;
    
    initiate = strfind(tempname2,'name=') + 5;
    terminate = strfind(tempname2,';');
    firstterminate = find(terminate>initiate);
    terminate_actual = terminate(firstterminate(1))-1;
      
    
    dataoutput_vague{i,1} = tempname2(initiate:terminate_actual);
    
    dataoutput_vague{i,2} = data2(genenum,1);
    dataoutput_vague{i,3} = data2(genenum,2);
    dataoutput_vague{i,4} = data2(genenum,3);
    dataoutput_vague{i,5} = data2(genenum,4);
    dataoutput_vague{i,6} = data2(genenum,5);
    dataoutput_vague{i,7} = data2(genenum,6);
    dataoutput_vague{i,8} = data2(genenum,7);
    dataoutput_vague{i,9} = data2(genenum,8);
    dataoutput_vague{i,10} = data2(genenum,9);
    dataoutput_vague{i,11} = data2(genenum,10);
    dataoutput_vague{i,12} = data2(genenum,11);
end

disp('Saving data')
save TGACGTCAgenes_1kb_promoter_vague_new.mat