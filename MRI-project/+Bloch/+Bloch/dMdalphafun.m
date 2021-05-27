function dMdalpha = dMdalphafun(alpha0,TRF,m0s,T1,T2f,R,T2s_,B0,B1,TR,RFflag)
%DMDALPHAFUN
%    DMDALPHA = DMDALPHAFUN(ALPHA0,TRF,M0S,T1,T2F,R,T2S_,B0,B1,TR,RFFLAG)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    02-Feb-2021 18:02:02

t2 = 1.0./TRF;
t3 = B1.*RFflag.*t2;
t4 = 1.0./T2s_;
t5 = TRF.*t4;
t6 = log(t5);
t7 = t6.^2;
t8 = t7.^2;
t9 = t8.^2;
t10 = t9.^2;
t11 = B1.^2;
t12 = RFflag.^2;
t13 = 1.0./TRF.^2;
t14 = t6.*4.734209035397098e-1;
t15 = t7.*2.142944865470107e-1;
t16 = t6.*t7.*6.168206376076885e-2;
t17 = t8.*3.405803678011424e-3;
t18 = t6.*t7.*t8.*2.510379993273342e-3;
t19 = t9.*7.252948004778164e-4;
t20 = t6.*t7.*t9.*3.462515593522486e-5;
t21 = t8.*t9.*5.75925759532023e-7;
t22 = t7.*t8.*t9.*1.597466862656782e-7;
t23 = t6.*t7.*t8.*t9.*4.950447003497e-8;
t24 = t7.*t10.*1.440717778396328e-10;
t25 = t8.*t10.*4.035255635129169e-13;
t28 = t6.*t8.*1.35200107780919e-2;
t29 = t7.*t8.*4.230807355630016e-3;
t30 = t6.*t9.*3.617855549952223e-4;
t31 = t7.*t9.*5.339997364891677e-5;
t32 = t6.*t8.*t9.*1.919567796635778e-6;
t33 = t10.*9.135634780140215e-9;
t34 = t6.*t10.*8.805182533780258e-11;
t35 = t6.*t7.*t10.*1.337098474521444e-11;
t26 = t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25-t28-t29-t30-t31-t32-t33-t34-t35+4.922812113644589e-1;
t27 = RFflag.*t2;
dMdalpha = reshape([0.0,0.0,-t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t27,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t27,0.0,0.0,0.0,0.0,0.0,0.0,0.0,T2s_.*alpha0.*t11.*t12.*t13.*t26.*-2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-T2s_.*(alpha0.*t11.*t12.*t13.*t26.*2.0-T2s_.*alpha0.*t11.*t12.*t13.*(t4.*4.734209035397098e-1+t4.*t6.*4.285889730940213e-1+t4.*t7.*1.850461912823065e-1-t4.*t8.*6.760005389045952e-2-t4.*t9.*3.256069994957001e-3-t4.*t10.*1.496881030742644e-9+t4.*t6.*t7.*1.36232147120457e-2-t4.*t6.*t8.*2.538484413378009e-2-t4.*t6.*t9.*5.339997364891677e-4+t4.*t7.*t8.*1.757265995291339e-2+t4.*t6.*t10.*2.593292001113391e-9+t4.*t7.*t9.*3.808767152874734e-4-t4.*t7.*t10.*2.540487101590744e-10-t4.*t8.*t9.*2.495438135626512e-5+t4.*t6.*t7.*t8.*5.802358403822531e-3+t4.*t6.*t7.*t9.*6.911109114384276e-6+t4.*t6.*t7.*t10.*8.070511270258337e-12+t4.*t6.*t8.*t9.*2.236453607719495e-6+t4.*t7.*t8.*t9.*7.4256705052455e-7-t4.*t6.*t7.*t8.*t9.*1.461701564822434e-7).*2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,B1.*T2s_.*alpha0.*t12.*t13.*t26.*-4.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[40,5]);
