function server_bloch_basis_optimization(ID,bss,cpu_num)
if cpu_num>0
    parpool(cpu_num);
end
bsi=ceil(ID/100);
basis_size=bss(bsi);
ID=mod(ID,100);
readfolder=['2020_ID_',num2str(ID),'_data'];
sweepnum=2;
%switch sweepnum
    %case 3
        writefile=['2020_ID_',num2str(ID),'_data/','optimized_basis_',num2str(basis_size)];
    %case 2
        %writefile=['2020_ID_',num2str(ID),'_data/','optimized_basis_sw2_',num2str(basis_size)];
%end
Bloch.fmincon_optimize_basis(...
            readfolder,writefile,...
            basis_size);
end

