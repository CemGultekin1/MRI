function server_bloch_basis_optimization(ID,cpu_num,extras)
if cpu_num>0
    parpool(cpu_num);
end

[T,bss,SARMAX,readfile,extras1]=server_opt_ID(ID);
extras=strcat(extras,' ', extras1);
jj=strfind(extras,'tag=');
jj=jj(end)+4;
if ~isempty(extras)
    writeaddress=[readfile,'_analysis_folder_',extras(jj:end)];
else
    writeaddress=[readfile,'_analysis_folder'];
end
readfolder=writeaddress;
writefile=strcat(readfile,'_bopt',string(bss));


Bloch.fmincon_optimize_basis(...
            readfolder,writefile,...
            bss,T,extras);
end

