function server_bloch_multi_statistics(ID,cpu_num,extras)
[T,bss,SARMAX,readfile,~]=server_opt_ID(ID);
load(readfile,'history');
x=history.control;
if isfield(history,'basis')
    basis=history.basis;
else
    basis=[];
end

jj=strfind(extras,'tag=')+4;
if ~isempty(extras)
    writeaddress=[readfile,'_analysis_folder_',extras(jj:end)];
else
    writeaddress=[readfile,'_analysis_folder'];
end
S=dir(writeaddress);
if isempty(S)
    mkdir(writeaddress);
end


samplenum=5000;
batch_size=500;
if cpu_num>0
    parpool(cpu_num);
end

[QQ,RR]=qr(basis,0);
basis={QQ,RR};

batchnum=ceil(samplenum/batch_size);
for i=1:batchnum
    disp(['computing batch : ' , num2str(i),'/',num2str(batchnum), '...']);    
    CRBvals=Bloch.statistics((i-1)*batch_size,min(samplenum,i*batch_size),x,basis,extras);
    disp('       done');
    saveloc=[writeaddress,'/','CRBvals_',num2str(i),'.mat'];
    save(saveloc,'CRBvals');
end
end

