function server_bloch_multi_statistics(ID,cpu_num)
[T,bss,SARMAX]=server_opt_ID(ID);
Tstr=strrep(num2str(T),'.','p');
sarstr=num2str(SARMAX,'%.2e');
sarstr=strrep(strrep(sarstr,'.','p'),'+','');
readfile=['2021_01_12_T',Tstr,'_b',num2str(bss),'_SAR',sarstr,'*'];
sample_num=5000;
server_bloch_analysis(readfile,sample_num,500,cpu_num);
end

