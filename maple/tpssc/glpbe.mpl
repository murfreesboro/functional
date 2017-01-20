#GGA  rho
###############################################################################
# REF
# this is the Gorling-Levy limit of PBE correlation energy
# see the paper of PSTS functional, the Appendix A:
# Density functional with full exact exchange, balanced nonlocality 
# of correlation, and constraint satisfaction
# Perdew, J.P. and Staroverov, V.N. and Tao, J. and Scuseria, G.E.
# note:
# for the constants used in this program, for example; the gammas, beta and 
# omega etc., I just used the value in this Appendix A. We note that they
# may not accurate enough.
# END
###############################################################################

pbec := proc (rhoT,xi,GAA,GBB,GAB)

# some constants
gammas:= (1-ln(2))/PI^2;
beta  := 0.066725;
omega := 0.046644;
c     := ((3*PI*PI)/16)^(1/3);
chi   := (beta/gammas)*c*c*exp((-omega/gammas));

# variables definition
z     := xi;
sigma := GAA + GBB + 2*GAB;
kf    := (3*(PI^2)*rhoT)^(1/3);
s     := sigma^(1/2)/(2*rhoT*kf);
phi   := (1/2)*((1+z)^(2/3)+(1-z)^(2/3));

# final expression
t1    := -gammas*phi^3;
t2    := chi*s*s/phi^2;
functional := t1*ln(1+1/(t2 + t2^2));
end proc;


# the real calculation
rhoT  := RA+RB;
xi    := (RA-RB)/rhoT;
pbe   := pbec(rhoT,xi,GAA,GBB,GAB);
pbe_a := pbec(RA,1,GAA,0,0);
pbe_b := pbec(RB,1,GBB,0,0);

functional_RA_RB  := pbe;
functional_RA     := pbe_a;
functional_RB     := pbe_b;


