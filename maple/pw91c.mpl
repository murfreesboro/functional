#GGA  rho
###########################################################################
#  REF
#  J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, 
#  M.R. Pederson, D.J. Singh and C. Fiolhais, 
#  “Atoms, molecules, solids and surfaces: Applications of the generalized
#   gradient approximation for exchange and correlation”, 
#  Phys. Rev. B, 46(11):6671–6687, 1992
#  END
###########################################################################

# here below is the functional per electron
pw91c := proc (rhoT,xi,GAA,GBB,GAB)
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
z     := xi;
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
gs    := (1/2)*((1+z)^(2/3)+(1-z)^(2/3));
sigma := GAA + GBB + 2*GAB;
T     := (1/2)*(sigma^(1/2))/(gs*ks*rhoT);
alpha := 0.09;
Cc0   := 0.004235;
Cx    := -0.001667;
nu    := (16/PI)*(3*PI^2)^(1/3);
beta  := nu*Cc0;
A     := (2*alpha/beta)*1/(exp((-2)*alpha*EC/(gs^3*beta^2)) - 1);
Cc    := (1/1000)*(2.568 + 23.266*rs + 0.007389*rs^2)/(1 + 8.723*rs + 0.472*rs^2 + 0.073890*rs^3) - Cx;
H0    := (1/2)*(gs^3*beta^2/alpha)*ln(1 + (2*alpha/beta)*(T^2 + A*T^4)/(1 + A*T^2 + A^2*T^4));
H1    := nu*(Cc - Cc0 - (3/7)*Cx)*gs^3*T^2*exp(-100*(gs^4*ks^2*T^2)/(kf^2));
functional1:= EC;
functional2:= H0 + H1;
functional := functional1 + functional2;
end proc;


# final results
rhoT    :=  RA+RB;
xi      :=  (RA-RB)/rhoT;
f1      :=  pw91c(rhoT,xi,GAA,GBB,GAB);
f2      :=  pw91c(RA  , 1,GAA,  0,  0);
f3      :=  pw91c(RB  , 1,GBB,  0,  0);

functional_RA_RB  := rhoT*f1;
functional_RA     := RA*f2;
functional_RB     := RB*f3;



