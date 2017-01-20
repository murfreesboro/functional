#META
######################################################################
# This is the exchange part for PKZB functional
# ref
# This is the exchange part for PKZB functional
# "Accurate density functional with correct formal properties: 
# A step beyond the generalized gradient approximation"
# J. P. Perdew, S. Kurth, A. Zupan, and P. Blaha
# Phys. Rev. Lett., 82 (1999) 2544-47
# end
####################################################################### 



#######################################################################
#
# this is the exchange functional evaluated at each spin component
#
#######################################################################
pkzb := proc (RX, GXX, tauX)

# define the variables
rho   := 2*RX;
gRho  := 4*GXX; 
# we use the big tau, hence here we do not multiply
# the factor of 1/2. However, in principle we should
tau   := tauX;

# parameters
d     := 0.113; 
kappa := 0.804;

# lSD part
kf     := (3*(PI^2)*rho)^(1/3);
lsd    := (-3/4)*(kf/PI);

# the correlation to the LSD
rho83 := rho^(8/3);
rho53 := rho^(5/3);
t1    := 4*(3*PI^2)^(2/3);
t2    := 2*(3*PI^2)^(2/3);
p     := gRho/(t1*rho83);
q     := 3*tau/(t2*rho53) - 9/20 - p/12;
x     := (10/81)*p + (146/2025)*q^2 - (73/405)*q*p + (d + (1/kappa)*(10/81)^2)*p^2;
Fxpq  := 1 + kappa - kappa/(1 + x/kappa); 

# result
functional := rho*lsd*Fxpq;
end proc;

f1 := pkzb(RA,GAA,TA);
f2 := pkzb(RB,GBB,TB);
functional_RA_RB := 0.5*(f1 + f2);
functional_RA    := 0.5*f1;
functional_RB    := 0.5*f2;


