% disp('Loading matrix')
% load('genome.mat')
% 
% mypool = parpool;
% % 
% data = zeros(17716,1);
% % 
% disp('Initiate analysis')
% parfor i = 1 : 17716
%     data(i) = seqwordcount([seq(i).Sequence],'TGACGTCA');
% end
% % 
% delete(mypool)
% % 
% disp('Saving data')
% save TGACGTCAgenes.mat

%%
genenums = find(data>0);

ngenes = length(genenums);

dataoutput = cell(ngenes,2);

for i = 1 : ngenes
    genenum = genenums(i);
    
    tempname = seq(genenum).Header;
    
    initiate = strfind(tempname,'name=') + 5;
    terminate = strfind(tempname,';');
    firstterminate = find(terminate>initiate);
    terminate_actual = terminate(firstterminate(1))-1;
      
    
    dataoutput{i,1} = tempname(initiate:terminate_actual);
    
    %dataoutput{i,2} = data(genenum);
end


