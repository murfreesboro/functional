#GGA  rho grho
######################################################################
# ref
# This is the C(zeta,xi) function used in TPSS functional
# Climbing the density functional ladder: Nonempirical 
# meta-generalized gradient approximation designed for 
# molecules and solids
# J. M. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria
# Phys. Rev. Lett., 91 (2003) 146401
# end
####################################################################### 


#######################################################################
#  The C(zeta,xi) function in TPSS functional
#######################################################################
c_zeta_xi := proc(zeta,xi2)

# some intermidiate expressions
Czeta0 :=  0.53 + 0.87*zeta^2 + 0.50*zeta^4 + 2.26*zeta^6;
t      :=  (1/2)*((1+zeta)^(-4/3)+(1-zeta)^(-4/3));
Czetaxi:=  Czeta0/(1+xi2*t)^4;
end proc;


#######################################################################
#  The C(zeta,0) function in TPSS functional
#  This is for the cases that only RA or RB exists
#  as zeta = 1 or -1, the (1+zeta)^{-4/3} will cause some error
#######################################################################
c_zeta_0 := proc(RA,RB)
Rho    :=  RA+RB;
zeta   :=  (RA-RB)/Rho;
Czeta0 :=  0.53 + 0.87*zeta^2 + 0.50*zeta^4 + 2.26*zeta^6;
end proc;



#######################################################################
#  the main body for the C(zeta) function
#######################################################################
Rho    :=  RA+RB;
zeta   :=  (RA-RB)/Rho;
onepz  :=  1+zeta;
onemz  :=  1-zeta;
gxi2   :=  (GAA*onemz*onemz + GBB*onepz*onepz -2*GAB*onemz*onepz)/Rho^2;
t      :=  (2*(3*PI^2*Rho)^(1/3))^2;
xi2    :=  gxi2/t;



# c(zeta,xi) function
functional_RA_RB   := c_zeta_xi(zeta,xi2);
functional_RA      := c_zeta_0(RA,RB);
functional_RB      := c_zeta_0(RA,RB);




