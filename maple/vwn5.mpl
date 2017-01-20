# LDA  rho
#########################################################
# ref
# S. H. Vosko, L. Wilk, and M. Nusair, 
# “Accurate spin-dependent electron liquid correlation 
#  energies for local spin density calculations: 
#  A critical analysis,”  the 5th functional
# Can. J. Phys., 58 (1980) 1200-11.
# end
#########################################################
vwn5 := proc(rhoT,xi)
x    := (3/(4*Pi*rhoT))^(1/6);
g_xi := (9/8)*( (1+xi)^(4/3) + (1-xi)^(4/3) - 2);
Q    := (4*c-b^2)^(1/2);
X    := t^2 + b*t + c;
eps  := A*( ln(x^2/subs(t=x, X)) + (2*b/Q)*arctan(Q/(2*x+b)) - \
	   b*x0/subs(t=x0, X)*(ln((x-x0)^2/subs(t=x,X)) + (2*(2*x0+b)/Q)*arctan(Q/(2*x+b))));
ECP  := subs(A=0.03109070000,  x0=-0.10498,    b=3.72744, c=12.9352, eps);
ECF  := subs(A=0.01554535000,  x0=-0.32500,    b=7.06042, c=18.0578, eps);
EAC  := subs(A=-1/(6*Pi^2),    x0=-0.00475840, b=1.13107, c=13.0045, eps);
fac  := 1 + (4/(9*(2^(1/3)-1)) * (ECF-ECP)/EAC -1)*xi^4;
fac1 := ECP + EAC*g_xi*fac;
end proc;

# result
rhoT  := RA+RB;
xi    := (RA-RB)/rhoT;
f1    := vwn5(rhoT,xi); 
f2    := vwn5(RA,  1);
# here if only RB exists, RB should replace RA's position
f3    := vwn5(RB,  1); 
functional_RA_RB := rhoT*f1;
functional_RA    := RA*f2;
functional_RB    := RB*f3;


