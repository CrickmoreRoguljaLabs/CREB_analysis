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
    data_TGACGTCA(i) = seqwordcount([seq(i).Sequence],'TGACGTCA');
    
    data_AGACGTCA(i) = seqwordcount([seq(i).Sequence],'AGACGTCA');
    data_GGACGTCA(i) = seqwordcount([seq(i).Sequence],'GGACGTCA');
    data_CGACGTCA(i) = seqwordcount([seq(i).Sequence],'CGACGTCA');
    
    data_TAACGTCA(i) = seqwordcount([seq(i).Sequence],'TAACGTCA');
    data_TTACGTCA(i) = seqwordcount([seq(i).Sequence],'TTACGTCA');
    data_TCACGTCA(i) = seqwordcount([seq(i).Sequence],'TCACGTCA');

    data_TGTCGTCA(i) = seqwordcount([seq(i).Sequence],'TGTCGTCA');
    data_TGGCGTCA(i) = seqwordcount([seq(i).Sequence],'TGGCGTCA');
    data_TGCCGTCA(i) = seqwordcount([seq(i).Sequence],'TGCCGTCA');
end
% 
delete(mypool)
% 

data = [data_TGACGTCA,...
    data_AGACGTCA,data_GGACGTCA,data_CGACGTCA,...
    data_TAACGTCA,data_TTACGTCA,data_TCACGTCA,...
    data_TGTCGTCA,data_TGGCGTCA,data_TGCCGTCA];
disp('Saving data')
save TGACGTCAgenes_vague.mat
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
