#LDA  rho
######################################################################
#  REF
#  J.P. Perdew and Y. Wang, 
#  “Accurate and simple analytic representation of the
#   electron-gas correlation energy”, 
#  Phys. Rev. B, 45(23):13244–13249, 1992
#  END
######################################################################

pw92c  := proc (rhoT, xi) 
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
# result
functional:= EC;
end proc;

# final functionals
rhoT  :=  RA + RB;
xi    :=  (RA - RB)/rhoT;
f1    :=  pw92c(rhoT,xi);
f2    :=  pw92c(RA, 1);
f3    :=  pw92c(RB, 1);
 
functional_RA_RB  := rhoT*f1;
functional_RA     := RA*f2;
functional_RB     := RB*f3;





