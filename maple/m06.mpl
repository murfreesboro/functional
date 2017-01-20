# META   SP
###############################################################################
# REF
# This is the M06-L functional, here below is the reference for M06 family
#
# Y. Zhao and D. G. Truhlar, 
# “A new local density functional for main-group thermochemistry, 
#  transition metal bonding, thermochemical kinetics, and noncovalent 
#  interactions,” 
# J. Chem. Phys., 125 (2006), 194101: 1-18. 
#
# Y. Zhao and D. G. Truhlar, 
# “Comparative DFT study of van der Waals complexes: Rare-gas dimers, 
#  alkaline-earth dimers, zinc dimer, and zinc-rare-gas dimers,” 
# J. Phys. Chem., 110 (2006) 5121-29.
#
# Y. Zhao and D. G. Truhlar, 
# “Density Functional for Spectroscopy: No Long-Range Self-Interaction 
#  Error, Good Performance for Rydberg and Charge-Transfer States, 
#  and Better Performance on Average than B3LYP for Ground States,” 
#  J. Phys. Chem. A, 110 (2006) 13126-30.
#
#  Y. Zhao and D. G. Truhlar, 
#  “The M06 suite of density functionals for main group thermochemistry, 
#   thermochemical kinetics, noncovalent interactions, excited states, 
#   and transition elements: two new functionals and systematic testing 
#   of four M06-class functionals and 12 other functionals,” 
#  Theor. Chem. Acc., 120 (2008) 215-41.
#
# END
###############################################################################


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
# zgvt4
#taken from the vs98 exchange fortran file
###############################################################################
gvt := proc(xg,zg,gama,a,b,c,d,e,f)   
   gvt4 := (a + b*xg + c*zg + d*xg*xg + e*zg*xg + f*zg*zg)/gama;
end proc;


###############################################################################
#vs98ss
#taken from the vs98 correlation fortran file
###############################################################################
vs98ss := proc(PX, GX, TX)
# parameters in vs98
r13 :=   0.4905945;
r14 :=  -0.1437348;
r15 :=   0.2357824;
r16 :=   0.001871015;
r17 :=  -0.003788963;
r18 :=   0.0000000;
gcc:=  0.00515088;
cf :=  9.115599720;
#cf :=  (3/5)*(3*PI*PI)^(2/3);



rhoo := PX;  
rho53 := rhoo^(5/3);
rho83 := rhoo^(8/3);

# lsda calculation;
EUEG  := lsda(rhoo, 1);

Chii := GX/rho83;
Z    := (TX/rho53) - cf;
kc   := 1 + gcc*(Chii + Z);
xk   := Chii/kc;
zk   := Z/kc;
DD   := 1 - Chii/(4*(Z + cf));
gac  := gvt(xk,zk,kc,r13,r14,r15,r16,r17,r18);
E    := DD*EUEG*gac;
end proc;


###############################################################################
# the exchange functional, we just put all things together
###############################################################################
vs98x := proc(RX,TX,GX)
# parameters
cf     := 9.115599720; 
Axlsda := -0.9305257363491;
gg     := 0.00186726;
f13    := 1/3;
f43    := 4/3;
f53    := 5/3;
f83    := 8/3; 
F113   := 11/3;

# adjust parameters
r1  := 0.1422057*Axlsda;
r2  := 0.0007370319*Axlsda;
r3  :=-0.01601373*Axlsda;
r4  := 0.00000000;
r5  := 0.00000000;
r6  := 0.00000000;

# variables for alpha
rhoo  := RX;
rho43 := rhoo^f43;  
rrho  := 1/rhoo;   
rho13 := rho43*rrho;
rho53 := rhoo^f53;
rho83 := rho53*rhoo;
tauu  := TX;
Gamma := GX;
x     := Gamma/rho83;
z     := tauu/rho53 - cf;
kx    := 1 + gg*x + gg*z;
xk    := x/kx;
zk    := z/kx;

#zgvt4
xg   := xk;
zg   := zk;
gama := kx;
ct   := gg;
ct2  := gg;
a    := r1;
b    := r2;
c    := r3;
d    := r4;
e    := r5;
f    := r6;
g    := gama;

functional  := rho43*a/g + rho43*b*xg/g + rho43*c*zg/g + rho43*d*xg*xg/g \
             + rho43*e*zg*xg/g + rho43*f*zg*zg/g ;
end proc;



###############################################################################
# vs98 correlation for both alpha and beta
###############################################################################
vs98c_ab := proc (RA,RB,GAA,GBB,TA,TB)

# constant
r7  :=-2.741539;
r8  :=-0.6720113;
r9  :=-0.07932688;
r10 := 0.001918681;
r11 :=-0.002032902;
r12 := 0.000000;
gab := 0.00304966;
cf  := 9.115599720;

rhoT := RA+RB;
zeta := (RA-RB)/rhoT;
ECT  := lsda(rhoT, zeta);
ECA  := lsda(RA, 1);
ECB  := lsda(RB, 1);
EC   := ECT-ECA-ECB;

rho53A := RA^(5/3);
rho83A := RA^(8/3);
rho53B := RB^(5/3);
rho83B := RB^(8/3);

XAB    := GAA/rho83A + GBB/rho83B;
ZAB    := (TA/rho53A) + (TB/rho53B)- 2*cf;
kab    := 1 + gab*(ZAB + XAB);
xk     := XAB/kab;
zk     := ZAB/kab;
gcab   := gvt(xk,zk,kab,r7,r8,r9,r10,r11,r12);
f5     := gcab*EC;

functional := f5; 
end proc;



###############################################################################
#M06ss
#taken from the M06 correlation fortran file
###############################################################################
M06ss := proc(PX, GX, TX)
# parameters in M06L
sss0 :=  0.5094055;
sss1 :=  -1.491085;
sss2 :=  17.23922;
sss3 :=  -38.59018;
sss4 :=  28.45044;
Css  :=  0.06;

# the real codes
EUEG  := lsda(PX, 1);
DD    := TX-(1/4)*GX/PX;
Chii  := GX/PX^(8/3);
U     := Css*Chii/(1 + Css*Chii);
W     := sss0+U*(sss1+U*(sss2+U*(sss3+U*sss4)));
Fscc  := DD/TX;
E     := Fscc*W*EUEG;
end proc;


###############################################################################
# m06x
# take from the M06 exchange fortran file
###############################################################################
m06x := proc(RX, GX, TX)   
# parameters
at0 :=  0.5877943;
at1 := -0.1371776;
at2 :=  0.2682367;
at3 := -2.5158980;
at4 := -2.9788920;
at5 :=  8.7106790;
at6 := 16.8819500;
at7 := -4.4897240;
at8 :=-32.9998300;
at9 :=-14.4905000;
at10:= 20.4374700;
at11:= 12.5650400;
fac := 1.0;
at  := 1.0;
C1  := 0.00336116;
C2  := 0.00449267;

# codes
fL  := fac;
fNL := fac;
Ax  := -(3/2)*((4/3)*PI)^(-(1/3));
rhoo:= RX;
rho43:= rhoo^(4/3);
rho13:= rhoo^(1/3);
rho53:= rhoo^(5/3);
tauN:= TX;
tauu:= tauN;
TauUEG:= (3/5)*((6*PI*PI)^(2/3))*rho53;
Tsig := TauUEG/tauN;
Wsig := (Tsig-1)/(Tsig+1);
W1:= Wsig;
W2:= Wsig*W1;
W3:= Wsig*W2;
W4:= Wsig*W3;
W5:= Wsig*W4;
W6:= Wsig*W5;
W7:= Wsig*W6;
W8:= Wsig*W7;
W9:= Wsig*W8;
W10:= Wsig*W9;
W11:= Wsig*W10;
Fsig:= at*(at0+ at1*W1+ at2*W2 + at3*W3 + at4*W4 + at5*W5 + at6*W6 + at7*W7 \
	   + at8*W8 + at9*W9 + at10*W10 + at11*W11);

Gamma2:= GX;
Gamma := (Gamma2)^(1/2);
x:= Gamma/rho43;
x2:= x*x;
En:= C1*x2;
Ed:= 1 + C2*x2;
E:= -En/Ed;

functional := rho43*(fL*Ax+fNL*E)*Fsig;
end proc;

###############################################################################
# here is the main routine for calculating the functional of M06 
###############################################################################
# alpha and beta part contribution in seperate way
f1  := M06ss(RA, GAA, TA);
f2  := M06ss(RB, GBB, TB);

f3  := m06x(RA,GAA,TA);
f4  := m06x(RB,GBB,TB);

# contribution from the VS98 functional
# correlation part
vs981  := vs98ss(RA, GAA, TA);
vs982  := vs98ss(RB, GBB, TB);
vs983  := vs98c_ab(RA,RB,GAA,GBB,TA,TB);

# exchange part
vs984  := vs98x(RA,TA,GAA);
vs985  := vs98x(RB,TB,GBB);

VS98_RA_RB := vs981 + vs982 + vs983 + vs984 + vs985;
VS98_RA    := vs981 + vs984;
VS98_RB    := vs982 + vs985;


# now consider the spin density, lsda part
rhoT := RA+RB;
zeta := (RA-RB)/rhoT;
ECT  := lsda(rhoT, zeta);
ECA  := lsda(RA, 1);
ECB  := lsda(RB, 1);
EC   := ECT-ECA-ECB;

# parameters
sopp0 := 3.741539;
sopp1 := 218.7098;
sopp2 :=-453.1252;
sopp3 := 293.6479;
sopp4 :=-62.87470;
COpp  := 0.0031;

#codes
ChiA  := GAA/RA^(8/3);
ChiB  := GBB/RB^(8/3);
U     := COpp*(ChiA+ChiB)/(1 + COpp*(ChiA+ChiB));
W     := sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)));
f5    := EC*W;

#result
functional_RA_RB := f1 + f2 + f3 + f4 + f5 + VS98_RA_RB;
functional_RA    := f1 + f3 + VS98_RA;
functional_RB    := f2 + f4 + VS98_RB;





