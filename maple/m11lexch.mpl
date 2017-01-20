# META
###############################################################################
# REF
# This is the M11L exchange functional, here below is the reference 
# for M08/M11 functionals family:
#
#  Ref: (a) Zhao, Y.  and Truhlar, D. G. JCTC, 2008, 4 , 1849          
#       (b) Peverati, R. and Truhlar, D. G. J.P.C.Lett. 2011, 2, 2810  
#       (c) Peverati, R. and Truhlar, D. G. J.P.C.Lett. 2012, 3, 117   
#
# END
###############################################################################

#######################################################################
#
# this is the short-range exchange functional for LSDA 
#
#######################################################################
lrclsda := proc(RX)

# constants
pi    := 3.1415926535897932384626433832795;
F1    := 1.0;
F2    := 2.0;
F3    := 3.0;
F4    := 4.0;
F5    := 5.0;
F6    := 6.0;
F7    := 7.0;
F8    := 8.0;
F9    := 9.0;
F10   := 10.0;
F11   := 11.0;
F48   := 48.0;
F81   := 81.0;
F1o2  := F1/F2;
F1o3  := F1/F3; 
F1o4  := F1/F4;
F2o3  := F2/F3;
F3o5  := F3/F5;
F4o3  := F4/F3; 
F4o9  := F4/F9;
F5o3  := F5/F3;
F8o3  := F8/F3;
PI12  := sqrt(pi);

# this is the omega value used in the range-separated functional
Emu   := 0.25;   

# real codes
AX    := -(F3/F2)*(F4o3*pi)**(-F1o3);
Cmu   := (F6*pi**F2)**F1o3;   
Rho13 := RX**F1o3;
Rho43 := RX**F4o3;
tmu   := Emu/(F2*Cmu*Rho13);
tmu2  := tmu*tmu;
tmu3  := tmu*tmu2;
tmu4  := tmu*tmu3;
W     := exp(-F1o4/tmu2);
ERFV  := erf(F1o2/tmu);
Fsr   := F1-F4o3*tmu*(-F6*tmu+F8*tmu3+W*(F4*tmu-F8*tmu3)+F2*PI12*ERFV);
functional  := AX*Rho43*Fsr;
end proc;

#######################################################################
# this is the exchange functional evaluated at each spin component
#######################################################################
m08m11exch := proc(RX, GXX, TX)

# constants
pi    := 3.1415926535897932384626433832795;
F1    := 1.0;
F2    := 2.0;
F3    := 3.0;
F4    := 4.0;
F5    := 5.0;
F6    := 6.0;
F7    := 7.0;
F8    := 8.0;
F9    := 9.0;
F10   := 10.0;
F11   := 11.0;
F48   := 48.0;
F81   := 81.0;
F1o3  := F1/F3; 
F2o3  := F2/F3;
F3o5  := F3/F5;
F4o3  := F4/F3; 
F5o3  := F5/F3;
Ax    := -(F3/F2)*(F4o3*pi)**(-F1o3); 

# PBE parameters
Mus   := F10/F81;
kapas := 0.552;
Mu    := 0.21951;
kapa  := 0.804;

# Parameters for M11-L
at00  :=   8.121131E-01;
at01  :=   1.738124E+01;
at02  :=   1.154007E+00;
at03  :=   6.869556E+01;
at04  :=   1.016864E+02;
at05  :=  -5.887467E+00;
at06  :=   4.517409E+01;
at07  :=  -2.773149E+00;
at08  :=  -2.617211E+01;
at09  :=   0.000000E+00;
at10  :=   0.000000E+00;
at11  :=   0.000000E+00;
bt00  :=   1.878869E-01;
bt01  :=  -1.653877E+01;
bt02  :=   6.755753E-01;
bt03  :=  -7.567572E+01;
bt04  :=  -1.040272E+02;
bt05  :=   1.831853E+01;
bt06  :=  -5.573352E+01;
bt07  :=  -3.520210E+00;
bt08  :=   3.724276E+01;
bt09  :=   0.000000E+00;
bt10  :=   0.000000E+00;
bt11  :=   0.000000E+00;
ct00  :=  -4.386615E-01;
ct01  :=  -1.214016E+02;
ct02  :=  -1.393573E+02;
ct03  :=  -2.046649E+00;
ct04  :=   2.804098E+01;
ct05  :=  -1.312258E+01;
ct06  :=  -6.361819E+00;
ct07  :=  -8.055758E-01;
ct08  :=   3.736551E+00;
ct09  :=   0.000000E+00;
ct10  :=   0.000000E+00;
ct11  :=   0.000000E+00;
dt00  :=   1.438662E+00;
dt01  :=   1.209465E+02;
dt02  :=   1.328252E+02;
dt03  :=   1.296355E+01;
dt04  :=   5.854866E+00;
dt05  :=  -3.378162E+00;
dt06  :=  -4.423393E+01;
dt07  :=   6.844475E+00;
dt08  :=   1.949541E+01;
dt09  :=   0.000000E+00;
dt10  :=   0.000000E+00;
dt11  :=   0.000000E+00;

# begin the codes
rhoo  := RX;
rho43 := rhoo**F4o3;
rho13 := rho43/rhoo;
rho53 := rhoo**F5o3;
TauN  := TX;
TauUEG:= F3o5*((F6*pi*pi)**F2o3)*rho53;
TSIG  := TauUEG/TauN;
Wsig  := (TSIG - F1)/(TSIG + F1);
W1    := Wsig;
W2    := Wsig*W1;
W3    := Wsig*W2;
W4    := Wsig*W3;
W5    := Wsig*W4;
W6    := Wsig*W5;
W7    := Wsig*W6;
W8    := Wsig*W7;
W9    := Wsig*W8;
W10   := Wsig*W9;
W11   := Wsig*W10;
Fsig1 := (at00 + at01*W1 + at02*W2 + at03*W3 + at04*W4 + at05*W5 + at06*W6 + at07*W7 + at08*W8 + at09*W9 + at10*W10+ at11*W11);
Fsig2 := (bt00 + bt01*W1 + bt02*W2 + bt03*W3 + bt04*W4 + bt05*W5 + bt06*W6 + bt07*W7 + bt08*W8 + bt09*W9 + bt10*W10+ bt11*W11);
Fsig3 := (ct00 + ct01*W1 + ct02*W2 + ct03*W3 + ct04*W4 + ct05*W5 + ct06*W6 + ct07*W7 + ct08*W8 + ct09*W9 + ct10*W10+ ct11*W11);
Fsig4 := (dt00 + dt01*W1 + dt02*W2 + dt03*W3 + dt04*W4 + dt05*W5 + dt06*W6 + dt07*W7 + dt08*W8 + dt09*W9 + dt10*W10+ dt11*W11); 
Gam   := GXX;
Gam12 := sqrt(Gam);
x     := Gam12/rho43;
s     := x/(F48*pi*pi)**F1o3;
y     := s*s;
Deno  := (F1 + Mu*y/kapa);
fx1   := F1+kapa*(F1-F1/Deno);
fx2   := F1+kapas*(F1-exp(-Mus*y/kapas));

# here different functionals may be different
ElSR  := lrclsda(rhoo);
ElLR  := Ax*rho43-ElSR;

# result
GGA1  := ElSR*fx1;
GGA2  := ElSR*fx2;
GGA3  := ElLR*fx1;
GGA4  := ElLR*fx2;
functional := GGA1*Fsig1+GGA2*Fsig2+GGA3*Fsig3+GGA4*Fsig4;
end proc;


#######################################################################
#  result
#######################################################################
f1 := m08m11exch(RA,GAA,TA);         
f2 := m08m11exch(RB,GBB,TB);         
functional_RA_RB :=  f1 + f2;
functional_RA    :=  f1;
functional_RB    :=  f2;


