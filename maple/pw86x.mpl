# GGA
#########################################################
# ref
# J. P. Perdew, Y. Wang, Phys. Rev. B, 33, 8800 (1986) 
# end
#########################################################

# PW86 exchange functional for each of spin components
pw86x := proc(RXX,GXX)
GX  := 4*GXX;
RX  := 2*RXX;
Ax  := -(3/4)*((3/PI)^(1/3));
SGX := GX^(1/2);
xa  := SGX/RX^(4/3);   
s   := (0.5/(3*PI*PI)^(1/3))*xa;
m   := 1/15;
b   := 14;
c   := 0.2;
f   := 1+(0.0864/m)*s*s+b*s*s*s*s+c*s^6;
F   := f^m;
fun := 0.5*Ax*RX^(4/3)*F;
end proc;

#result
alpha := pw86x(RA,GAA);
beta  := pw86x(RB,GBB);
functional_RA_RB  := alpha + beta;
functional_RA     := alpha;
functional_RB     := beta;


