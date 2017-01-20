# GGA   
#########################################################
# ref
# LYP evaluates the Lee-Yang-Parr correlation functional     
# and its derivatives.
# C. Lee, W. Yang, and R. G. Parr, 
# “Development of the Colle-Salvetti correlation-energy 
#  formula into a functional of the electron density,” 
# Phys. Rev. B, 37 (1988) 785-89.
# However, this is the revised version of LYP:
# B. Miehlich, A. Savin, H. Stoll and H. Preuss, 
# “Results obtained with the correlation
#  energy density functionals of becke and lee, 
#  yang and parr”, 
# Chem. Phys. Lett.,157(3):200–206, 1989
# end
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
OMEGA := exp(-c*rhoT^(-1/3))/(1 + Dd*rhoT^(-1/3))*rhoT^(-11/3);
DELTA := c*rhoT^(-1/3) + Dd*rhoT^(-1/3)/(1 + Dd*rhoT^(-1/3));
fac   := -A*B*OMEGA;
fac1  := -A*B*OMEGA*RA*RB;
LYP1  := -4*A*RA*RB/(1 + Dd*rhoT^(-1/3))*rhoT^(-1);
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

