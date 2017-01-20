#LDA  

Cs      :=  -(3/4)*(3/PI)^(1/3);
VA      :=  (2*RA)^(4/3);
VB      :=  (2*RB)^(4/3);
Rho     :=  2*(RA+RB);
Rhoa    :=  2*RA;
Rhob    :=  2*RB;

# result
functional_RA_RB := Cs*(VA+VB)/Rho;
functional_RA    := Cs*VA/Rhoa; 
functional_RB    := Cs*VB/Rhob; 



