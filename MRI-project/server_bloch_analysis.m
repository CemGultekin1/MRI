function server_bloch_analysis(x,basis,writeaddress,samplenum,batch_size,cpu_num,extras)

if cpu_num>0
    parpool(cpu_num);
end

[QQ,RR]=qr(basis,0);
basis={QQ,RR};

batchnum=ceil(samplenum/batch_size);
for i=1:batchnum
    disp(['computing batch : ' , num2str(i),'/',num2str(batchnum), '...']);
    if basis_opt
        CRBvals=Bloch.statistics((i-1)*batch_size,i*batch_size,x,basis,extras);
    else
        CRBvals=Bloch.statistics((i-1)*batch_size,i*batch_size,x,[],extras);
    end
    disp('       done');
    saveloc=[writeaddress,'/','CRBvals_',num2str(i),'.mat'];
    save(saveloc,'CRBvals');
end
end