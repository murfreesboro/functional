#GGA
###############################################################################
# REF
# J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, M.R. Pederson, 
# D.J. Singh and C. Fiolhais, 
# “Atoms, molecules, solids and surfaces: Applications of the generalized
#  gradient approximation for exchange and correlation”, 
# Phys. Rev. B, 46(11):6671–6687, 1992
# END
###############################################################################

# the functional form is pretty simple so that we do not have any functional any more.....
rhoA   := 2*RA;
gRhoa  := 4*GAA;
kfa    := (3*(PI^2)*rhoA)^(1/3);
exa    := (-3/4)*(kfa/PI);
sa     := (1/2)*(gRhoa^(1/2))/(kfa*rhoA);
FXa    := (1 + (0.19645)*sa*arcsinh(7.7956*sa) + (0.2743-0.1508*exp(-100*sa^2))*sa^2)/(1+ 0.19645*sa*arcsinh(7.7956*sa) + 0.004*sa^4);
rhoB   := 2*RB;
gRhob  := 4*GBB;
kfb    := (3*(PI^2)*rhoB)^(1/3);
exb    := (-3/4)*(kfb/PI);
sb     := (1/2)*(gRhob^(1/2))/(kfb*rhoB);
FXb    := (1 + (0.19645)*sb*arcsinh(7.7956*sb) + (0.2743-0.1508*exp(-100*sb^2))*sb^2)/(1+ 0.19645*sb*arcsinh(7.7956*sb) + 0.004*sb^4);
f1     := RA*exa*FXa;
f2     := RB*exb*FXb;
# result
functional_RA_RB :=  f1 + f2;
functional_RA    :=  f1;
functional_RB    :=  f2;



