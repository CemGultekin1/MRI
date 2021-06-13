function server_bloch_multi_statistics(ID,cpu_num,extras)

[T,bss,SARMAX,readfile,extras1]=server_opt_ID(ID);
extras=strcat(extras,' ', extras1);
load(readfile,'history');
x=history.control;
if isfield(history,'basis')
    basis=history.basis;
else
    basis=[];
end

jj=strfind(extras,'tag=');
jj=jj(end)+4;
if ~isempty(extras)
    writeaddress=[readfile,'_analysis_folder_',extras(jj:end)];
else
    writeaddress=[readfile,'_analysis_folder'];
end
S=dir(writeaddress);
if isempty(S)
    mkdir(writeaddress);
end


samplenum=5;%000;
batch_size=1;%500;
if cpu_num>0
    parpool(cpu_num);
end
if ~isempty(basis)
    [QQ,RR]=qr(basis,0);
    basis={QQ,RR};
end

batchnum=ceil(samplenum/batch_size);
CRBvals_=zeros(3,samplenum);
ii=0;
for i=1:batchnum
    disp(['computing batch : ' , num2str(i),'/',num2str(batchnum), '...']);    
    [fingerprints,CRBvals]=Bloch.statistics((i-1)*batch_size,min([samplenum,i*batch_size]),x,basis,extras);
    disp('       done');
    saveloc=[writeaddress,'/','CRBvals_',num2str(i),'.mat'];
    
    save(saveloc,'CRBvals');
    saveloc=[writeaddress,'/','fingerprints_',num2str(i),'.mat'];
    save(saveloc,'fingerprints');
    CRBvals_(:, ii+1:ii+size(CRBvals,2))=CRBvals;
    ii=ii+size(CRBvals,2);
end

CRBvals_=CRBvals_(:,1:ii);
CRBvals=CRBvals_;
saveloc=[writeaddress,'/','CRBvals.mat'];
save(saveloc,'CRBvals');
for i=1:batchnum
    saveloc=[writeaddress,'/','CRBvals_',num2str(i),'.mat'];
    delete(saveloc)
end


end

