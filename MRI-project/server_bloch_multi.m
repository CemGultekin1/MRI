function server_bloch_multi(ID,cpu_num)
[T,bss,SARMAX,extras]=server_opt_ID(ID);
readfile=[];
Tstr=strrep(num2str(T),'.','p');
sarstr=num2str(SARMAX,'%.2e');
sarstr=strrep(strrep(sarstr,'.','p'),'+','');

writefile=['2021_05_27_T',Tstr,'_b',num2str(bss),'_SAR',sarstr];
if ~isempty(extras)
    addthis=strfind(extras,'tag=');
    writefile=strcat(writefile, ['_',extras(addthis+4:end)]);
end
sample_num=1;
server_bloch_optimization(T,bss,SARMAX,sample_num,cpu_num,readfile,writefile,extras)
end