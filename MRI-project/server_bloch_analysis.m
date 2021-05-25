function server_bloch_analysis(fileroot,samplenum,batch_size,cpu_num)
S=dir(fileroot);
if isempty(S)
    return;
end
if length(S)>1
    ii=0;
    for i=1:length(S)
        if contains(S(i).name,'finalized')
            ii=i;
            break;
        end
    end
    filename=S( max(ii,1)).name;
else
    filename=S(1).name;
end
if cpu_num>0
    parpool(cpu_num);
end
if ~contains(filename,'finalized')
    load(filename,'history');
    x=history.control;
    if isfield(history,'basis')
        basis=history.basis;
    end
else
    load(filename,'x');
end
name=strrep(fileroot,'*','');
ii=strfind(filename,'.')-1;
writeaddress=[filename(1:ii(end)),'_analysis_folder'];
S=dir(writeaddress);
if isempty(S)
    mkdir(writeaddress);
end
basis_opt=0;
if ~isempty(basis)
    basis_opt=1;
elseif basis_size<inf
    basis_opt=1;
    basis=Bloch.init_basis(numsamples,TR,TRFmax,control,basis_size);
end

if basis_opt
   [QQ,RR]=qr(basis,0);
   basis={QQ,RR};
end
batchnum=ceil(samplenum/batch_size);
for i=1:batchnum
    disp(['computing batch : ' , num2str(i),'/',num2str(batchnum), '...']);
    if basis_opt
        CRBvals=Bloch.statistics(name,(i-1)*batch_size,i*batch_size,x,basis);
    else
        CRBvals=Bloch.statistics(name,(i-1)*batch_size,i*batch_size,x);
    end
    disp('       done');
    ii=strfind(filename,'.')-1;
    saveloc=[writeaddress,'/',filename(1:ii(end)),'_CRBvals_',num2str(i),'.mat'];
    save(saveloc,'CRBvals');
end
totcrb=[];
for i=1:batchnum
    ii=strfind(filename,'.')-1;
    saveloc=[writeaddress,'/',filename(1:ii(end)),'_CRBvals_',num2str(i),'.mat'];
    load(saveloc,'CRBvals');
    
    if i==1
       totcrb=CRBvals;
    else
       totcrb.CRB=cat(2,totcrb.CRB,CRBvals.CRB);
    end
end
rmdir(writeaddress,'s');
CRBvals=totcrb;
save([filename(1:ii(end)),'_CRBvals.mat'],'CRBvals');
end