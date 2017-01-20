# LDA  rho
#########################################################
# ref
# S. H. Vosko, L. Wilk, and M. Nusair, 
# “Accurate spin-dependent electron liquid correlation 
#  energies for local spin density calculations: 
#  A critical analysis,”  the 1th functional
# Can. J. Phys., 58 (1980) 1200-11.
#                                                                
# This uses the Random Phase Approximation (RPA) parameters     
# rather than those for the Ceperley-Alder solution,        
# and is commonly used in the B3LYP hybrid method.   
# end
#########################################################

vwn1rpa := proc (rhoT, xi)
x    := (3/(4*Pi*rhoT))^(1/6);
g_xi := (1/(2*2^(1/3)-2))*( (1+xi)^(4/3) + (1-xi)^(4/3) - 2);
Q    := (4*c-b^2)^(1/2);
X    := t^2 + b*t + c;
eps  := A*( ln(x^2/subs(t=x, X)) + (2*b/Q)*arctan(Q/(2*x+b)) - \
	   b*x0/subs(t=x0, X)*(ln((x-x0)^2/subs(t=x,X)) + (2*(2*x0+b)/Q)*arctan(Q/(2*x+b))));
ECP  := subs(A=0.03109070000,  x0=-0.409286,   b=13.0720, c=42.7198, eps);
ECF  := subs(A=0.01554535000,  x0=-0.743294,   b=20.1231, c=101.578, eps);
DEC  := (ECF-ECP)*g_xi; 
EC   := ECP + DEC;
end proc;


# result
rhoT := RA + RB;
xi   := (RA - RB)/rhoT;
f1   := vwn1rpa(rhoT,xi);
f2   := vwn1rpa(RA,1);
f3   := vwn1rpa(RB,1);

functional_RA_RB := rhoT*f1;
functional_RA    := RA*f2;
functional_RB    := RB*f3;



