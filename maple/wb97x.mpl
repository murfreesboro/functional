# GGA 
####################################################################################
# ref
# J.-D. Chai and M. Head-Gordon, 
# “Systematic optimization of long-range corrected hybrid density functionals,” 
# J. Chem. Phys., 128 (2008) 084106.
# This the the exchange part for WB97 functional
# end
####################################################################################

wb97x  := proc (RX, GXX)

# constants
F12   := 1/2;
F13   := 1/3; 
F14   := 1/4;
F23   := F13*2;
F43   := F13*4;
F49   := 4/9;
F83   := F43*2;
Cs    := (-1.5)*(F43*PI)^(-F13);
beta  := 0.004;
betaw := beta;
sqpi  := sqrt(PI);
Ckf   := (6*PI^2)^F13;
omega := 0.4;

# the LDA part
Rhoa   := RX;
wmua   := omega/(2*Ckf*Rhoa^F13);
wmu2a  := wmua*wmua;
wmu3a  := wmua*wmu2a;
wmu4a  := wmua*wmu3a;
Wa     := exp(-F14/wmu2a);
ERFVa  := erf( F12/wmua);
Fna    := 1 - F43*wmua*(-6*wmua+8*wmu3a+Wa*(4*wmua-8*wmu3a)+2*sqpi*ERFVa);

# the GGA part
srGa   := sqrt(GXX);
Gia    := 1/GXX;
srGia  := 1/srGa;
xa     := srGa/RX^(4/3);
S2a    := xa*xa;
PlbS2a := 1+betaw*S2a;
ua     := betaw*S2a/PlbS2a;
ha     := 1 + 1.13116*ua-2.749150*ua^2+12.09000*ua^3-5.71642*ua^4; 

# result
Exa    := Cs*RX^(4/3);
Gna    := Exa*ha;
Fa     := Gna*Fna;
f      := Fa;
end proc; 


# the final functional
f1     := wb97x(RA,GAA);
f2     := wb97x(RB,GBB);

functional_RA_RB := f1 + f2;
functional_RA    := f1;
functional_RB    := f2;


