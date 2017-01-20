#GGA
#########################################################
# REF
# J. P. Perdew, K. Burke, and M. Ernzerhof, 
# “Generalized gradient approximation made simple,” 
# Phys. Rev. Lett., 77 (1996) 3865-68.
# END
#########################################################

pbex := proc(RX,GXX)
kp     := 0.804;
mu     := 0.2195149727645171;
rhoA   := 2*RX;
gRhoa  := 4*GXX;
kfa    := (3*(PI^2)*rhoA)^(1/3);
exa    := (-3/4)*(kfa/PI);
sa     := (1/2)*(gRhoa^(1/2))/(kfa*rhoA);
FXa    := 1+kp-kp/(1 + mu*sa^2/kp);
f      := exa*FXa;
end proc;

# result
f1 := pbex (RA,GAA);         
f2 := pbex (RB,GBB);         
functional_RA_RB :=  RA*f1 + RB*f2;
functional_RA    :=  RA*f1;
functional_RB    :=  RB*f2;



