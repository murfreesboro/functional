#GGA  rho 
####################################################################################
# ref
# J.-D. Chai and M. Head-Gordon, 
# Systematic optimization of long-range corrected hybrid density functionals, 
# J. Chem. Phys., 128 (2008) 084106.
# This the the correlation part for WB97 functional
# end
####################################################################################

###############################################################################
# rho is the total density; and z is the polarized density z= ra-rb/(ra+rb)
# this is actually the pw92 correlation
###############################################################################
lsda := proc(rho, z)
pa  := 1;
Ba  := 0.0168869;
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
rs    := (3/4)^(1/3)*(PI*rho)^(-1/3);
fz    := ( (1+z)^(4/3) + (1-z)^(4/3) - 2)/(2^(4/3) - 2);
d2fz0 := 4/(9*(2^(1/3) - 1));
GC    := (-2)*B1*(1 + a1*r)*ln(1 + 1/(2*B1*(b1*r^(1/2) + b2*r + b3*r^(3/2) + b4*r^(p+1))));
BC    := subs(r=rs,B1=Ba, a1=a1a,b1=b1a,b2=b2a,b3=b3a,b4=b4a,p=pa, GC);
ECP   := subs(r=rs,B1=c0p,a1=a1p,b1=b1p,b2=b2p,b3=b3p,b4=b4p,p=pe, GC);
ECF   := subs(r=rs,B1=c0f,a1=a1f,b1=b1f,b2=b2f,b3=b3f,b4=b4f,p=pe, GC);
EC    := ECP - (BC*fz*(1-z^4))/d2fz0 + (ECF-ECP)*fz*z^4;
functional:= rho*EC;
end proc;


###############################################################################
# this is actually the pw92 correlation if only alpha or beta exist
# we note that in this case, in the original pw92 the z is 1 (only alpha) 
# or -1 (only beta); however; the functional is symmetrical between 1 and -1
# so that we set z = 1. Furthermore, fz = 1, and final EC is ECF
###############################################################################
ss_lsda := proc(RX)

# parameters
A    := 0.015545;
a1f  := 0.20548;
b1f  := 14.1189;
b2f  := 6.1977;
b3f  := 3.3662;
b4f  := 0.62517;
pe   := 1;

rho   := RX;
rs    := (3/4)^(1/3)*(PI*rho)^(-1/3);
GC    := (-2)*B1*(1 + a1*r)*ln(1 + 1/(2*B1*(b1*r^(1/2) + b2*r + b3*r^(3/2) + b4*r^(p+1))));
ECF   := subs(r=rs,B1=A,a1=a1f,b1=b1f,b2=b2f,b3=b3f,b4=b4f,p=pe, GC);
EC    := ECF;
functional:= rho*EC;
end proc;



###############################################################################
# GGA part of opsisite spin for wb97_c
###############################################################################
os_wb97_c := proc (xa2,xb2)

F12   := 1/2; 
beta  := 0.006;
S2    := F12*(xa2+xb2);
PlbS2 := 1+beta*S2;
u     := beta*S2/PlbS2;
h     := 1.00000+2.37031*u-11.399500*u^2+6.584050*u^3-3.78132*u^4;

# final result
functional := h;
end proc;


###############################################################################
# same spin for wb97_c
###############################################################################
ss_wb97_c := proc (RX,GXX)


# first step, calculate the ss_LSDA part
EC  := ss_lsda(RX);

# second step, calculate the GGA part
# constant
F13   := 1/3; 
F12   := 1/2; 
beta  := 0.2;

# real calculation
RX13  := RX^F13;
RX43  := RX*RX13;
srGXX := sqrt(GXX);
xx    := srGXX/RX43;
S2    := xx^2;
PlbS2 := 1+beta*S2;
u     := beta*S2/PlbS2;
h     := 1.00000-4.338790*u+18.230800*u^2-31.74300*u^3+17.2901*u^4;

# final result
functional := EC*h;
end proc;




###############################################################################
# the main part of correlation functional
# 
###############################################################################
# LDA part for opsite spin
# here we note that for the spin-polarized state, zeta =1 or -1 are qeuivalent
# in the final expression
rhoT       := RA+RB;
zeta       := (RA-RB)/rhoT;
ECT        := lsda(rhoT, zeta);
ECA        := ss_lsda(RA);
ECB        := ss_lsda(RB);
EC         := ECT-ECA-ECB;
EC_RA      := 0;
EC_RB      := 0;

# GGA part for opsite spin state
xa2        := GAA/RA^(8/3);
xb2        := GBB/RB^(8/3);
f0         := os_wb97_c(xa2,xb2);
f1         := EC*f0;

# same spin part
f2    := ss_wb97_c(RA,GAA);
f3    := ss_wb97_c(RB,GBB);


functional_RA_RB  := f1 + f2 + f3; 
functional_RA     := f2;
functional_RB     := f3;









