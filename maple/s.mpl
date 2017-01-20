# LDA
#####################################################################################
# ref
# This is the first exchange functional
# see Slater's book:
# J. C. Slater, 
# "The Self-Consistent Field for Molecular and Solids, 
#  Quantum Theory of Molecular and Solids", Vol. 4 
# McGraw-Hill, New York, 1974
# end
#####################################################################################
x := (-1)*(2/3)*(9/4)*(3/(4*PI))^(1/3);
functional_RA_RB := RA^(4/3)*x + RB^(4/3)*x;
functional_RA    := RA^(4/3)*x;
functional_RB    := RB^(4/3)*x;


