# GGA
#########################################################
# ref
# A. D. Becke, 
# “Density-functional exchange-energy approximation 
#  with correct asymptotic-behavior,” 
# Phys. Rev. A, 38 (1988) 3098.
# end
#########################################################

# becke 88 exchange functional for each of spin components
becke88 := proc(RX,GXX)
fac := -(3/2)*(3/(4*PI))^(1/3); 
lda := RX^(4/3)*fac;
g   := -(b*x^2)/(1+6*b*x*arcsinh(x));
SGXX:= GXX^(1/2);
xa  := SGXX/RX^(4/3);   
g_xa:= subs(b=0.0042, x=xa, g);   
f   := lda + RX^(4/3)*g_xa; 
end proc;

#result
alpha := becke88(RA,GAA);
beta  := becke88(RB,GBB);
functional_RA_RB  := alpha + beta;
functional_RA     := alpha;
functional_RB     := beta;


