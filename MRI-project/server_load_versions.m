
past_versions=struct('TR',cell(4,1),'control',cell(4,1),'name',cell(4,1));
sprintf('v6')
load('/gpfs/data/asslaenderlab/20200928_Phantom_NewSweeping_MT_v6p2/control_MT_v6p2_TR3p5ms_discretized.mat');
past_versions(1).TR=3.5e-3;
past_versions(1).control=control;
past_versions(1).name='v6';
sprintf('v5')
load('/gpfs/data/asslaenderlab/20200917_InVivo_MT_1mm_MWI_1p7mm/control_MT_v5_TR3p5ms_discretized.mat');
past_versions(2).TR=3.5e-3;
past_versions(2).control=control;
past_versions(2).name='v5';
sprintf('v3')
load('/gpfs/data/asslaenderlab/20200806_MT_inVivo/control_MT_v3p2_TR3p5ms_discretized.mat');
past_versions(3).TR=3.5e-3;
past_versions(3).control=control;
past_versions(3).name='v3';

sprintf('v0')
load('+Bloch/defaultControl.mat');
past_versions(4).TR=4.5e-3;
past_versions(4).control=defaultControl;
past_versions(4).name='v0';

save('past_versions','past_versions');