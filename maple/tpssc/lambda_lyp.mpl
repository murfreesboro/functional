# GGA   lambda_deriv
#########################################################
# ref
# The lambda-dependent LYP functional is shown in the 
# MCY functional, see:
# Mori-Sanchez, P. and Cohen, A.J. and Yang, W.
# Self-interaction-free exchange-correlation functional 
# for thermochemistry and kinetics
# The Journal of chemical physics, 2006,124, 091102  
# end
# lambda = 0.69D0
#########################################################

rhoT  := RA + RB;
sigma := GAA + GBB + 2*GAB;
rhoGA := RA*GAA/rhoT;
rhoGB := RB*GBB/rhoT;
A     := 0.04918;
B     := 0.132:
c     := 0.2533;
Dd    := 0.349;
CF    := (3/10)*(3*PI^2)^(2/3);
OMEGA := exp(-c*lambda*rhoT^(-1/3))/(1 + Dd*lambda*rhoT^(-1/3))*rhoT^(-11/3);
DELTA := c*lambda*rhoT^(-1/3) + Dd*lambda*rhoT^(-1/3)/(1 + Dd*lambda*rhoT^(-1/3));
fac   := -A*B*OMEGA;
fac1  := -A*B*OMEGA*RA*RB;
LYP1  := -4*A*RA*RB/(1 + lambda*Dd*rhoT^(-1/3))*rhoT^(-1);
LYP2  :=  fac1*2^(11/3)*CF*(RA^(8/3)+RB^(8/3));
LYP3  :=  fac1*(47/18-(7/18)*DELTA)*sigma;
LYP4  := -fac1*(5/2-(1/18)*DELTA)*(GAA+GBB);
LYP5  := -fac1*((DELTA-11)/9)*(rhoGA + rhoGB);
LYP6  := -fac*(2/3)*(rhoT^2)*sigma;
LYP7  :=  fac*((2/3)*rhoT^2 - RA^2)*GBB;
LYP8  :=  fac*((2/3)*rhoT^2 - RB^2)*GAA;
# result
functional_RA_RB := LYP1 + LYP2 + LYP3 + LYP4 + LYP5 + LYP6 + LYP7 + LYP8;
functional_RA    := LYP6 + LYP8;
functional_RB    := LYP6 + LYP7;

