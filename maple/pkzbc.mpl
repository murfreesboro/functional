#META
######################################################################
# ref
# This is the correlation part for PKZB functional
# "Accurate density functional with correct formal properties: 
# A step beyond the generalized gradient approximation"
# J. P. Perdew, S. Kurth, A. Zupan, and P. Blaha
# Phys. Rev. Lett., 82 (1999) 2544-47
# end
####################################################################### 



#######################################################################
#
# this is the correlation functional for PBE
#
#######################################################################
pbec := proc(RA,RB,GAA,GBB,GAB)
pa  := 1;
Aa  := 0.0168869;
a1a := 0.11125;
b1a := 10.357;
b2a := 3.6231;
b3a := 0.88026;
b4a := 0.49671;
pe  := 1;
c0p := 0.0310907;
a1p := 0.21370;
b1p := 7.5957;
b2p := 3.5876;
b3p := 1.6382;
b4p := 0.49294;
c0f := 0.01554535;
a1f := 0.20548;
b1f := 14.1189;
b2f := 6.1977;
b3f := 3.3662;
b4f := 0.62517; 
rhoT  := RA + RB;
z     := (RA - RB)/rhoT;
rs    := (3/4)^(1/3)*(PI*rhoT)^(-1/3);
fz    := ( (1+z)^(4/3) + (1-z)^(4/3) - 2)/(2^(4/3) - 2);
d2fz0 := 4/(9*(2^(1/3) - 1));
GC    := (-2)*A1*(1 + a1*r)*ln(1 + 1/(2*A1*(b1*r^(1/2) + b2*r + b3*r^(3/2) + b4*r^(p+1))));
AC    := subs(r=rs,A1=Aa, a1=a1a,b1=b1a,b2=b2a,b3=b3a,b4=b4a,p=pa, GC);
ECP   := subs(r=rs,A1=c0p,a1=a1p,b1=b1p,b2=b2p,b3=b3p,b4=b4p,p=pe, GC);
ECF   := subs(r=rs,A1=c0f,a1=a1f,b1=b1f,b2=b2f,b3=b3f,b4=b4f,p=pe, GC);
EC    := ECP - (AC*fz*(1-z^4))/d2fz0 + (ECF-ECP)*fz*z^4;
kf    := (3*(PI^2)*rhoT)^(1/3);
ks    := 2*(kf/PI)^(1/2);
phi   := (1/2)*((1+z)^(2/3)+(1-z)^(2/3));
sigma := GAA + GBB + 2*GAB;
T     := (1/2)*(sigma^(1/2))/(phi*ks*rhoT);
beta  := 0.06672455060314922;
gammas:= (1-ln(2))/PI^2;
A     := (beta/gammas)/(exp((-EC)/(gammas*(phi^3))) -1);
H     := (phi^3*gammas)*ln(1 + (beta/gammas)*T^2*(1 + A*T^2)/(1 + A*T^2 + A^2*T^4));
# result
functional := EC + H;
end proc;



#######################################################################
#  the main body for the correlation functional of PKZB
#######################################################################
# constant
c      :=  0.53;

# here we have to multiply 1/2 to make the big tau into the small
# tau, which is its right expression
TTA    :=  (1/2)*TA;
TTB    :=  (1/2)*TB;

Rho    :=  RA+RB;
Wa     :=  (1/8)*GAA/RA;
Wb     :=  (1/8)*GBB/RB;
WTa    :=  (Wa/TTA)^2;
WTb    :=  (Wb/TTB)^2;
pbe    :=  pbec(RA,RB,GAA,GBB,GAB);

# here for spin polarized part, the alpha and beta are different. However,
# the expression for pbeb is not correct if we write it as:
# pbeb   :=  pbec(0,RB, 0, GBB,0) 
pbea   :=  pbec(RA,0,GAA,0,0); 
pbeb   :=  pbec(RB,0,GBB,0,0); 


WT     :=  ((Wa+Wb)/(TTA+TTB))^2;
term1  :=  Rho*pbe*(1+c*WT);
term1RA:=  RA*pbea*(1+c*WTa);
term1RB:=  RB*pbeb*(1+c*WTb);
term2  :=  (1+c)*(WTa*RA*pbea+WTb*RB*pbeb);
term2RA:=  (1+c)*WTa*RA*pbea;
term2RB:=  (1+c)*WTb*RB*pbeb;

# result
functional_RA_RB := term1 - term2;
functional_RA    := term1RA - term2RA;
functional_RB    := term1RB - term2RB;


