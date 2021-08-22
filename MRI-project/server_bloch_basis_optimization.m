function server_bloch_basis_optimization(ID,cpu_num,extras1)
if cpu_num>0
    parpool(cpu_num);
end

[T,bss,SARMAX,readfile,extras]=server_opt_ID(ID);
extras=strcat(extras,' ', extras1);
if ~isempty(extras)
    jj=strfind(extras,'tag=');
    jj=jj(end)+4;
    writeaddress=[readfile,'_analysis_folder_',extras(jj:end)];
else
    writeaddress=[readfile,'_analysis_folder'];
end
readfolder=writeaddress;
if ~isempty(extras1)
    jj=strfind(extras1,'tag=');
    jj=jj(end)+4;

    writefile=strcat(readfile,'_bopt',string(bss),'_',extras1(jj:end));
else
    writefile=strcat(readfile,'_bopt',string(bss));
end
Bloch.fmincon_optimize_basis(...
            readfolder,writefile,...
            bss,T,extras1);
end

