# META  SP
###############################################################################
#  REF
#  T. Van Voorhis and G. E. Scuseria, 
#  “A never form for the exchange-correlation energy functional,” 
#  J. Chem. Phys., 109 (1998) 400-10.
#  END
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
r13 :=   0.3270912;
r14 :=  -0.03228915;
r15 :=  -0.02942406;
r16 :=   0.002134222;
r17 :=  -0.005451559;
r18 :=   0.01577575;
gcc :=  0.00515088;
cf  :=  9.115599720;
#cf :=  (3/5)*(3*PI*PI)^(2/3);

rhoo  := PX;  
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
# this is for the pure spin density, X=alpha or beta
###############################################################################
vs98x  := proc(RX,TX,GX)
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
r1  :=-0.9800683;
r2  :=-0.003556788;
r3  := 0.006250326;
r4  :=-0.00002354518;
r5  :=-0.0001282732;
r6  := 0.0003574822;

# variables
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
             + rho43*e*zg*xg/g + rho43*f*zg*zg/g; 
end proc;

###############################################################################
# here is the main routine for calculating the functional of VS98 correlation
###############################################################################
# alpha and beta part contribution in seperate way
f1  := vs98ss(RA, GAA, TA);
f2  := vs98ss(RB, GBB, TB);

# exchange part
f3  := vs98x(RA,TA,GAA);
f4  := vs98x(RB,TB,GBB);

# now consider the spin density, the parameters are different from
# the ones in the pure density evaluation!
r7  :=0.7035010;
r8  :=0.007694574;
r9  :=0.05152765;
r10 :=0.00003394308;
r11 :=-0.001269420;
r12 :=0.001296118;
gab :=0.00304966;
cf  :=9.115599720;

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

# result
functional_RA_RB := f1 + f2 + f3 + f4 + f5;
functional_RA    := f1 + f3;
functional_RB    := f2 + f4;

