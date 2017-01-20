#META  
######################################################################
# ref
# This is the WT term used in the TPSS functional, which is:
# we note that it's slightly different with the WT form of functional
# used in PKZB
# see the original TPSS paper for more information
# Climbing the density functional ladder: Nonempirical 
# meta-generalized gradient approximation designed for 
# molecules and solids
# J. M. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria
# Phys. Rev. Lett., 91 (2003) 146401
# end
####################################################################### 


# here we have to multiply 1/2 to make the big tau into the small
# tau, which is its right expression
TTA    :=  (1/2)*TA;
TTB    :=  (1/2)*TB;
Rho    :=  RA+RB;
Grho   :=  GAA + GBB + 2*GAB;
WT     :=  (1/8)*Grho/Rho;     
Wa     :=  (1/8)*GAA/RA;
Wb     :=  (1/8)*GBB/RB;

# result
functional_RA_RB := WT/(TTA+TTB);
functional_RA    := WTa/TTA; 
functional_RB    := WTb/TTB; 



