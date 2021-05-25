function server_bloch_fpg_gen(ID,sample_num,batch_size,cpu_num)
[T,bss,SARMAX,x,TR,TRFmax]=server_opt_ID(ID);
if cpu_num>0
    parpool(cpu_num);
end
writeaddress=['2020_ID_',num2str(ID),'_data'];
S=dir(writeaddress);
if isempty(S)
    mkdir(writeaddress);
end
batchnum=ceil(sample_num/batch_size);
for i=1:batchnum
    disp(['computing batch : ' , num2str(i),'/',num2str(batchnum), '...']);
    data=Bloch.fingerprint_gen((i-1)*batch_size,i*batch_size,x,TR,TRFmax);
    disp('       done');
    saveloc=[writeaddress,'/','fingerprints_',num2str(i),'.mat'];
    save(saveloc,'data');
end
end