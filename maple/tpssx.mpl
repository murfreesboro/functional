#META
######################################################################
# ref
# This is the exchange part for TPSS functional
# "Climbing the density functional ladder: Nonempirical 
#  meta-generalized gradient approximation designed for 
#  molecules and solids"
#  J. M. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria
#  Phys. Rev. Lett., 91 (2003) 146401
# end
####################################################################### 



#######################################################################
#
# this is the exchange functional evaluated at each spin component
#
#######################################################################
tpssx := proc (RX, GXX, tauX)

# parameters
c     := 1.59096; 
kappa := 0.804;
e     := 1.537;
b     := 0.40;
mu    := 0.21951;


# define the variables
rho   := 2*RX;
gRho  := 4*GXX; 
# we use the big tau, hence here we do not multiply
# the factor of 1/2. However, in principle we should
tau   := tauX; 


# some intermidiate expressions
rho83 := rho^(8/3);
w     := (1/8)*gRho/rho;
z     := w/tau;
p     := gRho/(4*(3*PI^2)^(2/3)*rho83);
alpha := (5/3)*p*(1/z-1);
qb    := ((9/20)*(alpha-1))/(1+b*alpha*(alpha-1))^(1/2) + 2*p/3;


# lSD part
kf     := (3*(PI^2)*rho)^(1/3);
lsd    := (-3/4)*(kf/PI);

# compoenents in x expression
t1    := (10/81 + c*z^2/(1+z^2)^2)*p;
t2    := (146/2025)*qb^2;
t3    := -(73/405)*qb*(1/2*((3/5)*z)^2 + (1/2)*p^2)^(1/2);
t4    := (1/kappa)*((10/81)^2)*p^2;
t5    := 2*(e^(1/2))*(10/81)*((3/5)*z)^2;
t6    := e*mu*p^3;
t7    := (1 + p*e^(1/2))^2;
x     := (t1+t2+t3+t4+t5+t6)/t7;

# the correlation to the LSD
Fxpq  := 1 + kappa - kappa/(1 + x/kappa); 

# result
functional := rho*lsd*Fxpq;
end proc;

f1 := tpssx(RA,GAA,TA);
f2 := tpssx(RB,GBB,TB);
functional_RA_RB  := 0.5*(f1 + f2);
functional_RA     := 0.5*f1;
functional_RB     := 0.5*f2;


