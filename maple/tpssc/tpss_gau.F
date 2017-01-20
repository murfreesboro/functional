*Deck TPSSc
      Subroutine gaussian_TPSSc(MaxDer,NGrid,sRhoA,sRhoB,
     $  D1RhoA,D1RhoB,sTauA,sTauB,F,D1F,D2F)
      Implicit Real*8(A-H,O-Z)
C
C     Correlation part of the TPSS meta-GGA functional
C     [1] J. Tao, J.P. Perdew, V.N. Staroverov, and G.E. Scuseria,
C         Phys. Rev. Lett., to be published
      Logical AccPar, DoA, DoB, DoAB, SpUA, SpUB
      Dimension sRhoA(*), sRhoB(*), D1RhoA(NGrid,3), D1RhoB(NGrid,3),
     $  sTauA(*), sTauB(*), F(*), D1F(NGrid,*), D2F(NGrid,*), ScaDFX(4)
      Save CfD, Zero, One, Two, Three, Four, Six, F18, TauTol
      Data CfD/2.8d0/, Zero/0.d0/, One/1.d0/, Two/2.d0/, Three/3.d0/,
     $  Four/4.d0/, Six/6.d0/, F18/0.125d0/, TauTol/1.d-13/
C
      Tiny = 1.0D-14
      MaxDer = 1
      If(MaxDer.lt.0.or.MaxDer.gt.2)
     $  Call QCrash('Invalid value of MaxDer in TPSSc')
      ScaDFX(1) = 1.0D0
      ScaDFX(2) = 1.0D0
      ScaDFX(3) = 1.0D0
      ScaDFX(4) = 1.0D0
      If(Abs(ScaDFX(3)-ScaDFX(4)).gt.Tiny)
     $  Call QCrash('TPSS correlation requires equal scale factors')
      AccPar = .True.
      Pi = Four*ATan(One)
      Pi34 = Three/(Four*Pi)
      F13 = One/Three
      F43 = Four/Three

C     determine the shell type case
      IOPCL = 0
      CALL RemGet(NALPHA, REM_NALPHA)
      CALL RemGet(NBETA,  REM_NBETA)
      NDEN=1
      CALL RemGet(JUSTAL,  REM_JUSTAL)
      IF (NALPHA .NE. NBETA .OR. JUSTAL .EQ. 1) NDEN=2
      IF (NDEN .EQ. 2) IOPCL = 1


C
      Do 100 i = 1, NGrid
        RhoA = sRhoA(i)
        RhoB = sRhoB(i)
        Rho = RhoA + RhoB
        If(Rho.gt.Tiny) then
          Rs = (Pi34/Rho)**F13
          Zeta = (RhoA - RhoB)/Rho
          TauA = sTauA(i)/Two
          TauB = sTauB(i)/Two
          Tau = TauA + TauB
          YA = D1RhoA(i,1)**2 + D1RhoA(i,2)**2 + D1RhoA(i,3)**2
          YB = D1RhoB(i,1)**2 + D1RhoB(i,2)**2 + D1RhoB(i,3)**2
          YC = D1RhoA(i,1)*D1RhoB(i,1) + D1RhoA(i,2)*D1RhoB(i,2)
     $       + D1RhoA(i,3)*D1RhoB(i,3)
          Y = (D1RhoA(i,1) + D1RhoB(i,1))**2
     $      + (D1RhoA(i,2) + D1RhoB(i,2))**2
     $      + (D1RhoA(i,3) + D1RhoB(i,3))**2
          GRhoA = Sqrt(YA)
          GRhoB = Sqrt(YB)
          GRho = Sqrt(Y)
          TauW = F18*Y/Rho
C
          DoA = RhoA.gt.Tiny .and. TauA.gt.TauTol
          DoB = RhoB.gt.Tiny .and. TauB.gt.TauTol
          DoAB = DoA .and. DoB
C
C         Evaluate coefficient C(Zeta,xi) and its derivatives
C
          Call EvCfC(MaxDer,IOpCl,Rho,Zeta,YA,YB,YC,C,dCdR,dCdZ,
     $      dCdYA,dCdYB,dCdYC,C2RR,C2RZ,C2RYA,C2RYB,C2RYC,C2ZZ,C2ZYA,
     $      C2ZYB,C2ZYC,C2YAYA,C2YAYB,C2YAYC,C2YBYB,C2YBYC,C2YCYC)
C
C         Spin-unpolarized contribution to the rev-PKZB eps_c
C
          Call PW91Lc(AccPar,Rs,Zeta,PLoc,dLdS,dLdZ,d2LdSS,
     $         d2LdSZ,d2LdZZ,xx,xx,xx,xx)
          Call PBEH0(MaxDer,Rho,GRho,Rs,Zeta,PLoc,dLdS,dLdZ,
     $         d2LdSS,d2LdSZ,d2LdZZ,H,dHdR,dHdG,dHdZ,
     $         d2HdRR,d2HdRG,d2HdGG,d2HdRZ,d2HdGZ,d2HdZZ)
          P = PLoc + H
C
C         Spin-polarized contributions to the rev-PKZB eps_c
C
          S = Zero
          If(DoA) then
            RsA = (Pi34/RhoA)**F13
            ZetaA = One
            Call PW91Lc(AccPar,RsA,ZetaA,PLocA,dLAdS,dLAdZ,d2LASS,
     $           d2LASZ,d2LAZZ,xx,xx,xx,xx)
            Call PBEH0(MaxDer,RhoA,GRhoA,RsA,ZetaA,PLocA,dLAdS,dLAdZ,
     $           d2LASS,d2LASZ,d2LAZZ,HA,dHAdR,dHAdG,dHAdZ,
     $           HA2RR,HA2RG,HA2GG,HA2RZ,HA2GZ,HA2ZZ)
            PA = PLocA + HA
            SpUA = PA.lt.P
            If(SpUA) then
              S = S + RhoA*P
            else
              S = S + RhoA*PA
            endIf
          endIf
          If(DoB) then
            RsB = (Pi34/RhoB)**F13
            ZetaB = -One
            Call PW91Lc(AccPar,RsB,ZetaB,PLocB,dLBdS,dLBdZ,d2LBSS,
     $           d2LBSZ,d2LBZZ,xx,xx,xx,xx)
            Call PBEH0(MaxDer,RhoB,GRhoB,RsB,ZetaB,PLocB,dLBdS,dLBdZ,
     $           d2LBSS,d2LBSZ,d2LBZZ,HB,dHBdR,dHBdG,dHBdZ,
     $           HB2RR,HB2RG,HB2GG,HB2RZ,HB2GZ,HB2ZZ)
            PB = PLocB + HB
            SpUB = PB.lt.P
            If(SpUB) then
              S = S + RhoB*P
            else
              S = S + RhoB*PB
            endIf
          endIf
C
C         Evaluate rev-PKZB eps_c and its derivatives
C
          RT = TauW/Tau
          U = RT*RT
          Q = One + C*U
          D = -(One + C)/Rho
          X = D*U
          E = P*Q + S*X
C
          If(MaxDer.ge.1) then
            dRsdR = -F13*Rs/Rho
            dWdR = -TauW/Rho
            dWdY = F18/Rho
            dUdW = Two*RT/Tau
            dUdR = dUdW*dWdR
            dUdY = dUdW*dWdY
            dUdT = -Two*U/Tau
            dGdY = One/(Two*GRho)
            dPdY = dHdG*dGdY
            dQdT = C*dUdT
            dXdT = D*dUdT
            dEdT = P*dQdT + S*dXdT
            If(DoA) then
              dZdRA = Two*RhoB/(Rho*Rho)
              dLdRA = dLdS*dRsdR + dLdZ*dZdRA
              dHdRA = dHdR + dHdZ*dZdRA
              dPdRA = dLdRA + dHdRA
              dCdRA = dCdR + dCdZ*dZdRA
              dQdRA = dCdRA*U + C*dUdR
              dQdYA = dCdYA*U + C*dUdY
              dDdRA = -(dCdRA + D)/Rho
              dDdYA = -dCdYA/Rho
              If(SpUA) then
                dSdRA = P + RhoA*dPdRA
                dSdYA = RhoA*dPdY
              else
                dRsAdR = -F13*RsA/RhoA
                dPAdR = dLAdS*dRsAdR + dHAdR
                dSdRA = PA + RhoA*dPAdR
                dGAdY = One/(Two*GRhoA)
                dPAdY = dHAdG*dGAdY
                dSdYA = RhoA*dPAdY
              endIf
              If(SpUB) then
                dSdRA = dSdRA + RhoB*dPdRA
                dSdYA = dSdYA + RhoB*dPdY
              endIf
              dXdRA = dDdRA*U + D*dUdR
              dXdYA = dDdYA*U + D*dUdY
              dEdRA = dPdRA*Q + P*dQdRA + dSdRA*X + S*dXdRA
              dEdYA = dPdY*Q + P*dQdYA + dSdYA*X + S*dXdYA
            endIf
            If(DoB) then
              dZdRB = -Two*RhoA/(Rho*Rho)
              dLdRB = dLdS*dRsdR + dLdZ*dZdRB
              dHdRB = dHdR + dHdZ*dZdRB
              dPdRB = dLdRB + dHdRB
              dCdRB = dCdR + dCdZ*dZdRB
              dQdRB = dCdRB*U + C*dUdR
              dQdYB = dCdYB*U + C*dUdY
              dDdRB = -(dCdRB + D)/Rho
              dDdYB = -dCdYB/Rho
              If(SpUB) then
                dSdRB = P + RhoB*dPdRB
                dSdYB = RhoB*dPdY
              else
                dRsBdR = -F13*RsB/RhoB
                dPBdR = dLBdS*dRsBdR + dHBdR
                dSdRB = PB + RhoB*dPBdR
                dGBdY = One/(Two*GRhoB)
                dPBdY = dHBdG*dGBdY
                dSdYB = RhoB*dPBdY
              endIf
              If(SpUA) then
                dSdRB = dSdRB + RhoA*dPdRB
                dSdYB = dSdYB + RhoA*dPdY
              endIf
              dXdRB = dDdRB*U + D*dUdR
              dXdYB = dDdYB*U + D*dUdY
              dEdRB = dPdRB*Q + P*dQdRB + dSdRB*X + S*dXdRB
              dEdYB = dPdY*Q + P*dQdYB + dSdYB*X + S*dXdYB
            endIf
            If(DoAB) then
              dQdYC = dCdYC*U + Two*C*dUdY
              dDdYC = -dCdYC/Rho
              dXdYC = dDdYC*U + Two*D*dUdY
              dSdYC = Zero
              If(SpUA) dSdYC = dSdYC + Two*RhoA*dPdY
              If(SpUB) dSdYC = dSdYC + Two*RhoB*dPdY
              dEdYC = Two*dPdY*Q + P*dQdYC + dSdYC*X + dXdYC*S
            endIf
C
            If(MaxDer.ge.2) then
              d2RsRR = -F43*dRsdR/Rho
              Rho3 = Rho*Rho*Rho
              d2UdWW = Two/(Tau*Tau)
              d2UdWT = -Four*RT/(Tau*Tau)
              d2UdTT = -Three*dUdT/Tau
              d2GdYY = -dGdY/(Two*GRho*GRho)
              d2PdYY = d2HdGG*dGdY**2 + dHdG*d2GdYY
              d2QdTT = C*d2UdTT
              d2WdRR = Y/(Four*Rho3)
              d2WdRY = -dWdY/Rho
              d2UdRR = d2UdWW*dWdR**2 + dUdW*d2WdRR
              d2UdRY = d2UdWW*dWdR*dWdY + dUdW*d2WdRY
              d2UdRT = d2UdWT*dWdR
              d2UdYY = d2UdWW*dWdY**2
              d2UdYT = d2UdWT*dWdY
              d2XdTT = D*d2UdTT
              d2EdTT = P*d2QdTT + S*d2XdTT
              If(DoA) then
                d2ZdAA = -Four*RhoB/Rho3
                d2LdAA = d2LdSS*dRsdR**2 + Two*d2LdSZ*dRsdR*dZdRA
     $                 + dLdS*d2RsRR + d2LdZZ*dZdRA**2 + dLdZ*d2ZdAA
                d2HdAA = d2HdRR + Two*d2HdRZ*dZdRA
     $                 + d2HdZZ*dZdRA**2 + dHdZ*d2ZdAA
                P2RARA = d2LdAA + d2HdAA
                C2RARA = C2RR + Two*C2RZ*dZdRA + C2ZZ*dZdRA**2
     $                 + dCdZ*d2ZdAA
                Q2RARA = C2RARA*U + Two*dCdRA*dUdR + C*d2UdRR
                D2RARA = -(C2RARA + Two*dDdRA)/Rho
                X2RARA = D2RARA*U + Two*dDdRA*dUdR + D*d2UdRR
                If(SpUA) then
                  S2RARA = Two*dPdRA + RhoA*P2RARA
                else
                  RsA2RR = -F43*dRsAdR/RhoA
                  d2LARR = d2LASS*dRsAdR**2 + dLAdS*RsA2RR
                  PA2RR = d2LARR + HA2RR
                  S2RARA = Two*dPAdR + RhoA*PA2RR
                endIf
                If(SpUB) S2RARA = S2RARA + RhoB*P2RARA
                E2RARA = P2RARA*Q + Two*dPdRA*dQdRA + P*Q2RARA
     $                 + S2RARA*X + Two*dSdRA*dXdRA + S*X2RARA
                P2RAY = (d2HdRG + d2HdGZ*dZdRA)*dGdY
                C2RAYA = C2RYA + C2ZYA*dZdRA
                Q2RAYA = C2RAYA*U + dCdRA*dUdY + dCdYA*dUdR + C*d2UdRY
                D2RAYA = -(C2RAYA + dDdYA)/Rho
                X2RAYA = D2RAYA*U + dDdRA*dUdY + dDdYA*dUdR + D*d2UdRY
                If(SpUA) then
                  S2RAYA = dPdY + RhoA*P2RAY
                else
                  PA2RY = HA2RG*dGAdY
                  S2RAYA = dPAdY + RhoA*PA2RY
                endIf
                If(SpUB) S2RAYA = S2RAYA + RhoB*P2RAY
                E2RAYA = P2RAY*Q + dPdRA*dQdYA + dPdY*dQdRA + P*Q2RAYA
     $                 + S2RAYA*X + dSdRA*dXdYA + dSdYA*dXdRA + S*X2RAYA
                Q2RAT = dCdRA*dUdT + C*d2UdRT
                X2RAT = dDdRA*dUdT + D*d2UdRT
                E2RAT = dPdRA*dQdT + P*Q2RAT + dSdRA*dXdT + S*X2RAT
                Q2YAYA = C2YAYA*U + Two*dCdYA*dUdY + C*d2UdYY
                D2YAYA = -C2YAYA/Rho
                X2YAYA = D2YAYA*U + Two*dDdYA*dUdY + D*d2UdYY
                If(SpUA) then
                  S2YAYA = RhoA*d2PdYY
                else
                  GA2YY = -dGAdY/(Two*YA)
                  PA2YY = HA2GG*dGAdY**2 + dHAdG*GA2YY
                  S2YAYA = RhoA*PA2YY
                endIf
                If(SpUB) S2YAYA = S2YAYA + RhoB*d2PdYY
                E2YAYA = d2PdYY*Q + Two*dPdY*dQdYA + P*Q2YAYA
     $                 + S2YAYA*X + Two*dSdYA*dXdYA + S*X2YAYA
                Q2YAT = dCdYA*dUdT + C*d2UdYT
                X2YAT = dDdYA*dUdT + D*d2UdYT
                E2YAT = dPdY*dQdT + P*Q2YAT + dSdYA*dXdT + S*X2YAT
              endIf
              If(DoB) then
                d2ZdBB = Four*RhoA/Rho3
                d2LdBB = d2LdSS*dRsdR**2 + Two*d2LdSZ*dRsdR*dZdRB
     $                 + dLdS*d2RsRR + d2LdZZ*dZdRB**2 + dLdZ*d2ZdBB
                d2HdBB = d2HdRR + Two*d2HdRZ*dZdRB
     $                 + d2HdZZ*dZdRB**2 + dHdZ*d2ZdBB
                P2RBRB = d2LdBB + d2HdBB
                C2RBRB = C2RR + Two*C2RZ*dZdRB + C2ZZ*dZdRB**2
     $                 + dCdZ*d2ZdBB
                Q2RBRB = C2RBRB*U + Two*dCdRB*dUdR + C*d2UdRR
                D2RBRB = -(C2RBRB + Two*dDdRB)/Rho
                X2RBRB = D2RBRB*U + Two*dDdRB*dUdR + D*d2UdRR
                If(SpUB) then
                  S2RBRB = Two*dPdRB + RhoB*P2RBRB
                else
                  RsB2RR = -F43*dRsBdR/RhoB
                  d2LBRR = d2LBSS*dRsBdR**2 + dLBdS*RsB2RR
                  PB2RR = d2LBRR + HB2RR
                  S2RBRB = Two*dPBdR + RhoB*PB2RR
                endIf
                If(SpUA) S2RBRB = S2RBRB + RhoA*P2RBRB
                E2RBRB = P2RBRB*Q + Two*dPdRB*dQdRB + P*Q2RBRB
     $                 + S2RBRB*X + Two*dSdRB*dXdRB + S*X2RBRB
                P2RBY = (d2HdRG + d2HdGZ*dZdRB)*dGdY
                C2RBYB = C2RYB + C2ZYB*dZdRB
                Q2RBYB = C2RBYB*U + dCdRB*dUdY + dCdYB*dUdR + C*d2UdRY
                D2RBYB = -(C2RBYB + dDdYB)/Rho
                X2RBYB = D2RBYB*U + dDdRB*dUdY + dDdYB*dUdR + D*d2UdRY
                If(SpUB) then
                  S2RBYB = dPdY + RhoB*P2RBY
                else
                  PB2RY = HB2RG*dGBdY
                  S2RBYB = dPBdY + RhoB*PB2RY
                endIf
                If(SpUA) S2RBYB = S2RBYB + RhoA*P2RBY
                E2RBYB = P2RBY*Q + dPdRB*dQdYB + dPdY*dQdRB + P*Q2RBYB
     $                 + S2RBYB*X + dSdRB*dXdYB + dSdYB*dXdRB + S*X2RBYB
                Q2RBT = dCdRB*dUdT + C*d2UdRT
                X2RBT = dDdRB*dUdT + D*d2UdRT
                E2RBT = dPdRB*dQdT + P*Q2RBT + dSdRB*dXdT + S*X2RBT
                Q2YBYB = C2YBYB*U + Two*dCdYB*dUdY + C*d2UdYY
                D2YBYB = -C2YBYB/Rho
                X2YBYB = D2YBYB*U + Two*dDdYB*dUdY + D*d2UdYY
                If(SpUB) then
                  S2YBYB = RhoB*d2PdYY
                else
                  GB2YY = -dGBdY/(Two*YB)
                  PB2YY = HB2GG*dGBdY**2 + dHBdG*GB2YY
                  S2YBYB = RhoB*PB2YY
                endIf
                If(SpUA) S2YBYB = S2YBYB + RhoA*d2PdYY
                E2YBYB = d2PdYY*Q + Two*dPdY*dQdYB + P*Q2YBYB
     $                 + S2YBYB*X + Two*dSdYB*dXdYB + S*X2YBYB
                Q2YBT = dCdYB*dUdT + C*d2UdYT
                X2YBT = dDdYB*dUdT + D*d2UdYT
                E2YBT = dPdY*dQdT + P*Q2YBT + dSdYB*dXdT + S*X2YBT
              endIf
              If(DoAB) then
                d2ZdAB = Two*Zeta/(Rho*Rho)
                d2LdAB = d2LdSS*dRsdR**2 + d2LdSZ*dRsdR*(dZdRA + dZdRB)
     $                 + dLdS*d2RsRR + d2LdZZ*dZdRA*dZdRB + dLdZ*d2ZdAB
                d2HdAB = d2HdRR + d2HdRZ*(dZdRA + dZdRB)
     $                 + d2HdZZ*dZdRA*dZdRB + dHdZ*d2ZdAB
                P2RARB = d2LdAB + d2HdAB
                C2RARB = C2RR + C2RZ*(dZdRA + dZdRB) + C2ZZ*dZdRA*dZdRB
     $                 + dCdZ*d2ZdAB
                Q2RARB = C2RARB*U + (dCdRA + dCdRB)*dUdR + C*d2UdRR
                D2RARB = -(C2RARB + dDdRA + dDdRB)/Rho
                X2RARB = D2RARB*U + (dDdRA + dDdRB)*dUdR + D*d2UdRR
                S2RARB = Zero
                If(SpUA) S2RARB = S2RARB + dPdRB + RhoA*P2RARB
                If(SpUB) S2RARB = S2RARB + dPdRA + RhoB*P2RARB
                E2RARB = P2RARB*Q + dPdRA*dQdRB + dPdRB*dQdRA + P*Q2RARB
     $                 + S2RARB*X + dSdRA*dXdRB + dSdRB*dXdRA + S*X2RARB
                C2RAYB = C2RYB + C2ZYB*dZdRA
                Q2RAYB = C2RAYB*U + dCdRA*dUdY + dCdYB*dUdR + C*d2UdRY
                D2RAYB = -(C2RAYB + dDdYB)/Rho
                X2RAYB = D2RAYB*U + dDdRA*dUdY + dDdYB*dUdR + D*d2UdRY
                S2RAYB = Zero
                If(SpUA) S2RAYB = S2RAYB + dPdY + RhoA*P2RAY
                If(SpUB) S2RAYB = S2RAYB + RhoB*P2RAY
                E2RAYB = P2RAY*Q + dPdRA*dQdYB + dPdY*dQdRA + P*Q2RAYB
     $                 + S2RAYB*X + dSdRA*dXdYB + dSdYB*dXdRA + S*X2RAYB
                C2RBYA = C2RYA + C2ZYA*dZdRB
                Q2RBYA = C2RBYA*U + dCdRB*dUdY + dCdYA*dUdR + C*d2UdRY
                D2RBYA = -(C2RBYA + dDdYA)/Rho
                X2RBYA = D2RBYA*U + dDdRB*dUdY + dDdYA*dUdR + D*d2UdRY
                S2RBYA = Zero
                If(SpUA) S2RBYA = S2RBYA + RhoA*P2RBY
                If(SpUB) S2RBYA = S2RBYA + dPdY + RhoB*P2RBY
                E2RBYA = P2RBY*Q + dPdRB*dQdYA + dPdY*dQdRB + P*Q2RBYA
     $                 + S2RBYA*X + dSdRB*dXdYA + dSdYA*dXdRB + S*X2RBYA
                C2RAYC = C2RYC + C2ZYC*dZdRA
                Q2RAYC = C2RAYC*U + dCdYC*dUdR
     $                 + Two*(dCdRA*dUdY + C*d2UdRY)
                D2RAYC = -(C2RAYC + dDdYC)/Rho
                X2RAYC = D2RAYC*U + dDdYC*dUdR
     $                 + Two*(dDdRA*dUdY + D*d2UdRY)
                S2RAYC = Two*S2RAYB
                E2RAYC = Two*(P2RAY*Q + dPdY*dQdRA) + dPdRA*dQdYC
     $                 + P*Q2RAYC + S2RAYC*X + dSdRA*dXdYC
     $                 + dSdYC*dXdRA + S*X2RAYC
                C2RBYC = C2RYC + C2ZYC*dZdRB
                Q2RBYC = C2RBYC*U + dCdYC*dUdR
     $                 + Two*(dCdRB*dUdY + C*d2UdRY)
                D2RBYC = -(C2RBYC + dDdYC)/Rho
                X2RBYC = D2RBYC*U + dDdYC*dUdR
     $                 + Two*(dDdRB*dUdY + D*d2UdRY)
                S2RBYC = Two*S2RBYA
                E2RBYC = Two*(P2RBY*Q + dPdY*dQdRB) + dPdRB*dQdYC
     $                 + P*Q2RBYC + S2RBYC*X + dSdRB*dXdYC
     $                 + dSdYC*dXdRB + S*X2RBYC
                Q2YAYB = C2YAYB*U + (dCdYA + dCdYB)*dUdY + C*d2UdYY
                D2YAYB = -C2YAYB/Rho
                X2YAYB = D2YAYB*U + (dDdYA + dDdYB)*dUdY + D*d2UdYY
                S2YAYB = Zero
                If(SpUA) S2YAYB = S2YAYB + RhoA*d2PdYY
                If(SpUB) S2YAYB = S2YAYB + RhoB*d2PdYY
                E2YAYB = d2PdYY*Q + dPdY*(dQdYA + dQdYB) + P*Q2YAYB
     $                 + S2YAYB*X + dSdYA*dXdYB + dSdYB*dXdYA + S*X2YAYB
                Q2YAYC = C2YAYC*U + dCdYC*dUdY
     $                 + Two*(dCdYA*dUdY + C*d2UdYY)
                D2YAYC = -C2YAYC/Rho
                X2YAYC = D2YAYC*U + dDdYC*dUdY
     $                 + Two*(dDdYA*dUdY + D*d2UdYY)
                S2YAYC = Two*S2YAYB
                E2YAYC = Two*(d2PdYY*Q + dPdY*dQdYA) + dPdY*dQdYC
     $                 + P*Q2YAYC + S2YAYC*X + dSdYA*dXdYC
     $                 + dSdYC*dXdYA + S*X2YAYC
                Q2YBYC = C2YBYC*U + dCdYC*dUdY
     $                 + Two*(dCdYB*dUdY + C*d2UdYY)
                D2YBYC = -C2YBYC/Rho
                X2YBYC = D2YBYC*U + dDdYC*dUdY
     $                 + Two*(dDdYB*dUdY + D*d2UdYY)
                S2YBYC = S2YAYC
                E2YBYC = Two*(d2PdYY*Q + dPdY*dQdYB) + dPdY*dQdYC
     $                 + P*Q2YBYC + S2YBYC*X + dSdYB*dXdYC
     $                 + dSdYC*dXdYB + S*X2YBYC
                Q2YCYC = C2YCYC*U + Four*(dCdYC*dUdY + C*d2UdYY)
                D2YCYC = -C2YCYC/Rho
                X2YCYC = D2YCYC*U + Four*(dDdYC*dUdY + D*d2UdYY)
                S2YCYC = Four*S2YAYB
                E2YCYC = Four*(d2PdYY*Q + dPdY*dQdYC) + P*Q2YCYC
     $                 + S2YCYC*X + Two*dSdYC*dXdYC + S*X2YCYC
                Q2YCT = dCdYC*dUdT + Two*C*d2UdYT
                X2YCT = dDdYC*dUdT + Two*D*d2UdYT
                E2YCT = Two*dPdY*dQdT + P*Q2YCT + dSdYC*dXdT + S*X2YCT
              endIf
            endIf
          endIf
C
C         TPSS correlation energy density and derivatives
C
          U = CfD*RT**3
          Q = One + E*U
          V = E*Q
          Fc = Rho*V
          F(i) = F(i) + ScaDFX(3)*Fc
          If(MaxDer.ge.1) then
            dUdW = Three*CfD*RT*RT/Tau
            dUdR = dUdW*dWdR
            dUdY = dUdW*dWdY
            dUdT = -Three*U/Tau
            dQdT = dEdT*U + E*dUdT
            dVdT = dEdT*Q + E*dQdT
            dFdT = Rho*dVdT
            If(DoA) then
              dQdRA = dEdRA*U + E*dUdR
              dQdYA = dEdYA*U + E*dUdY
              dVdRA = dEdRA*Q + E*dQdRA
              dFdRA = V + Rho*dVdRA
              dVdYA = dEdYA*Q + E*dQdYA
              dFdYA = Rho*dVdYA
              D1F(i,POS_RA)  = D1F(i,POS_RA)  + ScaDFX(3)*dFdRA
              D1F(i,POS_GAA) = D1F(i,POS_GAA) + ScaDFX(3)*dFdYA
              D1F(i,POS_TA)  = D1F(i,POS_TA)  + ScaDFX(3)*dFdT/Two
            endIf
            If(DoB) then
              dQdRB = dEdRB*U + E*dUdR
              dQdYB = dEdYB*U + E*dUdY
              dVdRB = dEdRB*Q + E*dQdRB
              dFdRB = V + Rho*dVdRB
              dVdYB = dEdYB*Q + E*dQdYB
              dFdYB = Rho*dVdYB
              D1F(i,POS_RB)  = D1F(i,POS_RB)  + ScaDFX(3)*dFdRB
              D1F(i,POS_GBB) = D1F(i,POS_GBB) + ScaDFX(3)*dFdYB
              D1F(i,POS_TB)  = D1F(i,POS_TB)  + ScaDFX(3)*dFdT/Two
            endIf
            If(DoAB) then
              dQdYC = dEdYC*U + Two*E*dUdY
              dVdYC = dEdYC*Q + E*dQdYC
              dFdYC = Rho*dVdYC
              D1F(i,POS_GAB) = D1F(i,POS_GAB) + ScaDFX(3)*dFdYC
            endIf
C
C  the second derivatives will be left for future
C            
C            
C            If(MaxDer.ge.2) then
C              d2UdWW = Six*CfD*RT/(Tau*Tau)
C              d2UdWT = -Three*dUdW/Tau
C              d2UdRR = d2UdWW*dWdR**2 + dUdW*d2WdRR
C              d2UdRY = d2UdWW*dWdR*dWdY + dUdW*d2WdRY
C              d2UdRT = d2UdWT*dWdR
C              d2UdYY = d2UdWW*dWdY**2
C              d2UdYT = d2UdWT*dWdY
C              d2UdTT = -Four*dUdT/Tau
C              d2QdTT = d2EdTT*U + Two*dEdT*dUdT + E*d2UdTT
C              d2VdTT = d2EdTT*Q + Two*dEdT*dQdT + E*d2QdTT
C              F2TT = Rho*d2VdTT
C              If(DoA) then
C                Q2RARA = E2RARA*U + Two*dEdRA*dUdR + E*d2UdRR
C                V2RARA = E2RARA*Q + Two*dEdRA*dQdRA + E*Q2RARA
C                F2RARA = Two*dVdRA + Rho*V2RARA
C                Q2RAYA = E2RAYA*U + dEdRA*dUdY + dEdYA*dUdR + E*d2UdRY
C                V2RAYA = E2RAYA*Q + dEdRA*dQdYA + dEdYA*dQdRA + E*Q2RAYA
C                F2RAYA = dVdYA + Rho*V2RAYA
C                Q2RAT = E2RAT*U + dEdRA*dUdT + dEdT*dUdR + E*d2UdRT
C                V2RAT = E2RAT*Q + dEdRA*dQdT + dEdT*dQdRA + E*Q2RAT
C                F2RAT = dVdT + Rho*V2RAT
C                Q2YAYA = E2YAYA*U + Two*dEdYA*dUdY + E*d2UdYY
C                V2YAYA = E2YAYA*Q + Two*dEdYA*dQdYA + E*Q2YAYA
C                F2YAYA = Rho*V2YAYA
C                Q2YAT = E2YAT*U + dEdYA*dUdT + dEdT*dUdY + E*d2UdYT
C                V2YAT = E2YAT*Q + dEdYA*dQdT + dEdT*dQdYA + E*Q2YAT
C                F2YAT = Rho*V2YAT
C                D2F(i,RARA) = D2F(i,RARA) + ScaDFX(3)*F2RARA
C                D2F(i,RAGA) = D2F(i,RAGA) + ScaDFX(3)*F2RAYA
C                D2F(i,RAKA) = D2F(i,RAKA) + ScaDFX(3)*F2RAT/Two
C                D2F(i,GAGA) = D2F(i,GAGA) + ScaDFX(3)*F2YAYA
C                D2F(i,GAKA) = D2F(i,GAKA) + ScaDFX(3)*F2YAT/Two
C                D2F(i,KAKA) = D2F(i,KAKA) + ScaDFX(3)*F2TT/Four
C              endIf
C              If(DoB) then
C                Q2RBRB = E2RBRB*U + Two*dEdRB*dUdR + E*d2UdRR
C                V2RBRB = E2RBRB*Q + Two*dEdRB*dQdRB + E*Q2RBRB
C                F2RBRB = Two*dVdRB + Rho*V2RBRB
C                Q2RBYB = E2RBYB*U + dEdRB*dUdY + dEdYB*dUdR + E*d2UdRY
C                V2RBYB = E2RBYB*Q + dEdRB*dQdYB + dEdYB*dQdRB + E*Q2RBYB
C                F2RBYB = dVdYB + Rho*V2RBYB
C                Q2RBT = E2RBT*U + dEdRB*dUdT + dEdT*dUdR + E*d2UdRT
C                V2RBT = E2RBT*Q + dEdRB*dQdT + dEdT*dQdRB + E*Q2RBT
C                F2RBT = dVdT + Rho*V2RBT
C                Q2YBYB = E2YBYB*U + Two*dEdYB*dUdY + E*d2UdYY
C                V2YBYB = E2YBYB*Q + Two*dEdYB*dQdYB + E*Q2YBYB
C                F2YBYB = Rho*V2YBYB
C                Q2YBT = E2YBT*U + dEdYB*dUdT + dEdT*dUdY + E*d2UdYT
C                V2YBT = E2YBT*Q + dEdYB*dQdT + dEdT*dQdYB + E*Q2YBT
C                F2YBT = Rho*V2YBT
C                D2F(i,RBRB) = D2F(i,RBRB) + ScaDFX(3)*F2RBRB
C                D2F(i,RBGB) = D2F(i,RBGB) + ScaDFX(3)*F2RBYB
C                D2F(i,RBKB) = D2F(i,RBKB) + ScaDFX(3)*F2RBT/Two
C                D2F(i,GBGB) = D2F(i,GBGB) + ScaDFX(3)*F2YBYB
C                D2F(i,GBKB) = D2F(i,GBKB) + ScaDFX(3)*F2YBT/Two
C                D2F(i,KBKB) = D2F(i,KBKB) + ScaDFX(3)*F2TT/Four
C              endIf
C              If(DoAB) then
C                Q2RARB = E2RARB*U + (dEdRA + dEdRB)*dUdR + E*d2UdRR
C                V2RARB = E2RARB*Q + dEdRA*dQdRB + dEdRB*dQdRA + E*Q2RARB
C                F2RARB = dVdRA + dVdRB + Rho*V2RARB
C                Q2RAYB = E2RAYB*U + dEdRA*dUdY + dEdYB*dUdR + E*d2UdRY
C                V2RAYB = E2RAYB*Q + dEdRA*dQdYB + dEdYB*dQdRA + E*Q2RAYB
C                F2RAYB = dVdYB + Rho*V2RAYB
C                Q2RAYC = E2RAYC*U + dEdYC*dUdR
C     $                 + Two*(dEdRA*dUdY + E*d2UdRY)
C                V2RAYC = E2RAYC*Q + dEdRA*dQdYC + dEdYC*dQdRA + E*Q2RAYC
C                F2RAYC = dVdYC + Rho*V2RAYC
C                Q2RBYA = E2RBYA*U + dEdRB*dUdY + dEdYA*dUdR + E*d2UdRY
C                V2RBYA = E2RBYA*Q + dEdRB*dQdYA + dEdYA*dQdRB + E*Q2RBYA
C                F2RBYA = dVdYA + Rho*V2RBYA
C                Q2RBYC = E2RBYC*U + dEdYC*dUdR
C     $                 + Two*(dEdRB*dUdY + E*d2UdRY)
C                V2RBYC = E2RBYC*Q + dEdRB*dQdYC + dEdYC*dQdRB + E*Q2RBYC
C                F2RBYC = dVdYC + Rho*V2RBYC
C                Q2YAYB = E2YAYB*U + (dEdYA + dEdYB)*dUdY + E*d2UdYY
C                V2YAYB = E2YAYB*Q + dEdYA*dQdYB + dEdYB*dQdYA + E*Q2YAYB
C                F2YAYB = Rho*V2YAYB
C                Q2YAYC = E2YAYC*U + (dEdYC + Two*dEdYA)*dUdY
C     $                 + Two*E*d2UdYY
C                V2YAYC = E2YAYC*Q + dEdYA*dQdYC + dEdYC*dQdYA + E*Q2YAYC
C                F2YAYC = Rho*V2YAYC
C                Q2YBYC = E2YBYC*U + (dEdYC + Two*dEdYB)*dUdY
C     $                 + Two*E*d2UdYY
C                V2YBYC = E2YBYC*Q + dEdYB*dQdYC + dEdYC*dQdYB + E*Q2YBYC
C                F2YBYC = Rho*V2YBYC
C                Q2YCYC = E2YCYC*U + Four*(dEdYC*dUdY + E*d2UdYY)
C                V2YCYC = E2YCYC*Q + Two*dEdYC*dQdYC + E*Q2YCYC
C                F2YCYC = Rho*V2YCYC
C                Q2YCT = E2YCT*U + dEdYC*dUdT + Two*dEdT*dUdY
C     $                 + Two*E*d2UdYT
C                V2YCT = E2YCT*Q + dEdYC*dQdT + dEdT*dQdYC + E*Q2YCT
C                F2YCT = Rho*V2YCT
C                D2F(i,RARB) = D2F(i,RARB) + ScaDFX(3)*F2RARB
C                D2F(i,RAGB) = D2F(i,RAGB) + ScaDFX(3)*F2RAYB
C                D2F(i,RAGC) = D2F(i,RAGC) + ScaDFX(3)*F2RAYC
C                D2F(i,RAKB) = D2F(i,RAKB) + ScaDFX(3)*F2RAT/Two
C                D2F(i,RBGA) = D2F(i,RBGA) + ScaDFX(3)*F2RBYA
C                D2F(i,RBGC) = D2F(i,RBGC) + ScaDFX(3)*F2RBYC
C                D2F(i,RBKA) = D2F(i,RBKA) + ScaDFX(3)*F2RBT/Two
C                D2F(i,GAGB) = D2F(i,GAGB) + ScaDFX(3)*F2YAYB
C                D2F(i,GAGC) = D2F(i,GAGC) + ScaDFX(3)*F2YAYC
C                D2F(i,GAKB) = D2F(i,GAKB) + ScaDFX(3)*F2YAT/Two
C                D2F(i,GBGC) = D2F(i,GBGC) + ScaDFX(3)*F2YBYC
C                D2F(i,GBKA) = D2F(i,GBKA) + ScaDFX(3)*F2YBT/Two
C                D2F(i,GCGC) = D2F(i,GCGC) + ScaDFX(3)*F2YCYC
C                D2F(i,GCKA) = D2F(i,GCKA) + ScaDFX(3)*F2YCT/Two
C                D2F(i,GCKB) = D2F(i,GCKB) + ScaDFX(3)*F2YCT/Two
C                D2F(i,KAKB) = D2F(i,KAKB) + ScaDFX(3)*F2TT/Four
C              endIf
CC
C
C
C            endIf
          endIf
        endIf
 100  Continue
      Return
      End




*Deck EvCfC
      Subroutine EvCfC(MaxDer,IOpCl,Rho,Z,YA,YB,YC,C,dCdR,dCdZ,dCdYA,
     $  dCdYB,dCdYC,C2RR,C2RZ,C2RYA,C2RYB,C2RYC,C2ZZ,C2ZYA,C2ZYB,C2ZYC,
     $  C2YAYA,C2YAYB,C2YAYC,C2YBYB,C2YBYC,C2YCYC)
      Implicit Real*8(A-H,O-Z)
C
C     Evaluates C(zeta,xi) and its derivatives for the TPSSc functional
C     Input: Rho(=RhoA+RhoB), Z(=zeta), YA = |GradA|^2,
C            YB = |GradB|^2, YC = GradA * GradB
C     Output: C, dCdZ, C2RR, C2RZ, ...
C
      Real*8 MDCutO
      Parameter (NPS=4)
      Save Zero, One, Two, Three, Four, Six, Seven, Eight, D11, D12,
     $  D30, Pt5, cf0, cf1, cf2, cf3, cfnu
      Data Zero/0.d0/, One/1.d0/, Two/2.d0/, Three/3.d0/, Four/4.d0/,
     $  Six/6.d0/, Seven/7.d0/, Eight/8.d0/, D11/11.d0/, D12/12.d0/,
     $  D30/30.d0/, Pt5/0.5d0/, cf0/0.53d0/, cf1/0.87d0/, cf2/0.5d0/,
     $  cf3/2.26d0/, cfnu/1.d0/
C
      C = cf0
      If(MaxDer.ge.1) then
        dCdR = Zero
        dCdZ = Zero
        dCdYA = Zero
        dCdYB = Zero
        dCdYC = Zero
        If(MaxDer.ge.2) then
          C2RR = Zero
          C2RZ = Zero
          C2RYA = Zero
          C2RYB = Zero
          C2RYC = Zero
          C2ZZ = Zero
          C2ZYA = Zero
          C2ZYB = Zero
          C2ZYC = Zero
          C2YAYA = Zero
          C2YAYB = Zero
          C2YAYC = Zero
          C2YBYB = Zero
          C2YBYC = Zero
          C2YCYC = Zero
        endIf
      endIf
      If(IOpCl.eq.0) Return
C
      Z2 = Z*Z
      A = cf0 + (cf1 + (cf2 + cf3*Z2)*Z2)*Z2
C
      OMZ = One - Z
      OPZ = One + Z
C      Small = MDCutO(0)
      Small = 1.0D-14
      If(OMZ.lt.Small .or. OPZ.lt.Small) then
        C = A
        Return
      endIf
C
      F13 = One/Three
      F43 = Four/Three
      F73 = Seven/Three
      Pi = Four*ATan(One)
      OMZ43 = OMZ**F43
      OPZ43 = OPZ**F43
      UM = Pt5/OMZ43
      UP = Pt5/OPZ43
      U = UM + UP
      dUMdZ = F43*UM/OMZ
      dUPdZ = -F43*UP/OPZ
      dUdZ = dUMdZ + dUPdZ
      d2UdZ = F73*(dUMdZ/OMZ - dUPdZ/OPZ)
      OMZ2 = OMZ**2
      OPZ2 = OPZ**2
      CFX = One/(Two*(Three*Pi**2*Rho)**F13*Rho)**2
      X = CFX*(OMZ2*YA - Two*OMZ*OPZ*YC + OPZ2*YB)
      S = One + cfnu*X*U
      B = S**NPS
      C = A/B
C
      If(MaxDer.ge.1) then
        F83 = Eight/Three
        PS = GFloat(NPS)
        dAdZ = (Two*cf1 + (Four*cf2 + Six*cf3*Z2)*Z2)*Z
        dXdR = -F83*X/Rho
        dXdZ = Two*CFX*(-OMZ*YA + Two*Z*YC + OPZ*YB)
        dCdS = -PS*C/S
        dSdX = cfnu*U
        dCdX = dCdS*dSdX
        dCdR = dCdX*dXdR
        dSdU = cfnu*X
        dSdZ = dSdX*dXdZ + dSdU*dUdZ
        dBdS = PS*S**(NPS-1)
        dBdZ = dBdS*dSdZ
        dCdZ = (dAdZ - C*dBdZ)/B
        dXdYA = CFX*OMZ2
        dXdYB = CFX*OPZ2
        dXdYC = -Two*CFX*OMZ*OPZ
        dCdYA = dCdX*dXdYA
        dCdYB = dCdX*dXdYB
        dCdYC = dCdX*dXdYC
C
        If(MaxDer.ge.2) then
          A2ZZ = Two*cf1 + (D12*cf2 + D30*cf3*Z2)*Z2
          X2RR = -F13*D11*dXdR/Rho
          C2SS = -GFloat(NPS+1)*dCdS/S
          C2XX = C2SS*dSdX**2
          C2RR = C2XX*dXdR**2 + dCdX*X2RR
          X2ZZ = Two*CFX*(YA + Two*YC + YB)
          S2ZZ = dSdX*X2ZZ + Two*cfnu*dXdZ*dUdZ + dSdU*d2UdZ
          B2SS = GFloat(NPS*(NPS-1))*S**(NPS-2)
          B2ZZ = B2SS*dSdZ*dSdZ + dBdS*S2ZZ
          C2ZZ = (A2ZZ - C*B2ZZ - Two*dBdZ*dCdZ)/B
          X2RZ = -F83*dXdZ/Rho
          C2SZ = -(PS*dCdZ + dCdS*dSdZ)/S
          S2XZ = cfnu*dUdZ
          C2XZ = C2SZ*dSdX + dCdS*S2XZ
          C2RZ = C2XZ*dXdR + dCdX*X2RZ
          C2YAYA = C2XX*dXdYA**2
          C2YBYB = C2XX*dXdYB**2
          C2YCYC = C2XX*dXdYC**2
          dSdYA = cfnu*U*dXdYA
          dSdYB = cfnu*U*dXdYB
          dSdYC = cfnu*U*dXdYC
          C2SYA = -(PS*dCdYA + dCdS*dSdYA)/S
          C2SYB = -(PS*dCdYB + dCdS*dSdYB)/S
          C2SYC = -(PS*dCdYC + dCdS*dSdYC)/S
          X2RYA = -F83*dXdYA/Rho
          X2RYB = -F83*dXdYB/Rho
          X2RYC = -F83*dXdYC/Rho
          C2RYA = (C2SYA*dXdR + dCdS*X2RYA)*dSdX
          C2RYB = (C2SYB*dXdR + dCdS*X2RYB)*dSdX
          C2RYC = (C2SYC*dXdR + dCdS*X2RYC)*dSdX
          X2ZYA = -Two*CFX*OMZ
          X2ZYB = Two*CFX*OPZ
          X2ZYC = Four*CFX*Z
          C2ZYA = C2XZ*dXdYA + dCdX*X2ZYA
          C2ZYB = C2XZ*dXdYB + dCdX*X2ZYB
          C2ZYC = C2XZ*dXdYC + dCdX*X2ZYC
          C2SYB = -(PS*dCdYB + dCdS*dSdYB)/S
          C2SYC = -(PS*dCdYC + dCdS*dSdYC)/S
          C2XYB = C2SYB*dSdX
          C2XYC = C2SYC*dSdX
          C2YAYB = C2XYB*dXdYA
          C2YAYC = C2XYC*dXdYA
          C2YBYC = C2XYC*dXdYB
        endIf
      endIf
C
      Return
      End


*Deck PW91Lc
      Subroutine PW91Lc(GoodPar,RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $  d2LdZZ,d3LSSS,d3LSSZ,d3LSZZ,d3LZZZ)
      Implicit Real*8(A-H,O-Z)
C
C     Evaluate the Perdew 91 local correlation functional and its
C     derivatives at one point.  GoodPar indicates whether to use
C     the mathematically exact values for parameters or the
C     values used in the definition of PW91 and PBE.
C
      Logical GoodPar
      Dimension ECLPar(6,3)
      Save F1, F2, F3, F4, F6, F8, F9, F12, F24, F36, GammaB, ECLPar,
     $  BadFZZ, BadEC13
      Data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/, F6/6.0d0/,
     $  F8/8.0d0/, F9/9.0d0/, F12/12.0d0/, F24/24.0d0/, F36/36.0d0/,
     $  GammaB/0.5198421D0/, BadFZZ/1.709921D0/,BadEC13/0.01688690D0/,
     $  ECLPar/0.03109070D0,0.21370D0, 7.5957D0,3.5876D0,1.6382D0,
     $         0.49294D0,
     $         0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     $         0.62517D0,
     $         0.01688690D0,0.11125D0,10.3570D0,3.6231D0,0.88026D0,
     $         0.49671D0/
C
      Pi = F4*ATan(F1)
      Third = F1 / F3
C     The exact values of FZZI and EClPar(1,3) yield results which
C     differ from the sloppy value used here and in the standard
C     PW91 example from Perdew -- MJF.
      If(GoodPar) then
        EClPar(1,3) = F1 / (F6*Pi*Pi)
        FZZI = F9*(F2**Third-F1) / F4
        GammaI = F1 / (F2*F2**Third-F2)
      else
        EClPar(1,3) = BadEC13
        FZZI = F1 / BadFZZ
        GammaI = F1 / GammaB
        endIf
      Call EvFZet(GammaI,Zeta,FZeta,dfZdz,d2fZdz,d3fZdz)
      Call EvPWLC(ECLPar(1,1),ECLPar(2,1),ECLPar(3,1),ECLPar(4,1),
     $  ECLPar(5,1),ECLPar(6,1),RS,EU,dEUdRS,d2UdRS,d3UdRS)
      Call EvPWLC(ECLPar(1,2),ECLPar(2,2),ECLPar(3,2),ECLPar(4,2),
     $  ECLPar(5,2),ECLPar(6,2),RS,EP,dEPdRS,d2PdRS,d3PdRS)
      Call EvPWLC(ECLPar(1,3),ECLPar(2,3),ECLPar(3,3),ECLPar(4,3),
     $  ECLPar(5,3),ECLPar(6,3),RS,AlphaM,dAMdRS,d2AdRS,d3AdRS)
      Z2 = Zeta*Zeta
      Z3 = Zeta*Z2
      Z4 = Zeta*Z3
      GZ = FZeta*Z4
      HZ = FZZI*(FZeta-GZ)
      PotLC = EU*(F1-GZ) + EP*GZ - AlphaM*HZ
      dLdS = dEUdRS*(F1-GZ) + dEPdRS*GZ - dAMdRS*HZ
      dGZ = dfZdz*Z4 + F4*FZeta*Z3
      dHZ = FZZI*(dFZdz-dGZ)
      dLdz = (EP-EU)*dGZ - AlphaM*dHZ
      d2GZ = d2fZdz*Z4 + F8*Z3*dfZdz + F12*FZeta*Z2
      d2HZ = FZZI*(d2FZdz-d2GZ)
      d2LdSS = d2UdRS*(F1-GZ) + d2PdRS*GZ - d2AdRS*HZ
      d2LdSZ = (dEPdRS-dEUdRS)*dGZ - dAMdRS*dHZ
      d2LdZZ = (EP-EU)*d2GZ - AlphaM*d2HZ
      d3GZ = d3fZdz*Z4+F12*d2fZdz*Z3+F36*Z2*dfZdz+F24*FZeta*Zeta
      d3HZ =  FZZI*(d3FZdz-d3GZ)
      d3LSSS = d3UdRS*(F1-GZ) + d3PdRS*GZ - d3AdRS*HZ
      d3LSSZ = (d2PdRS-d2UdRS)*dGZ - d2AdRS*dHZ
      d3LSZZ = (dEPdRS-dEUdRS)*d2GZ - dAMdRS*d2HZ
      d3LZZZ = (EP-EU)*d3GZ - AlphaM*d3HZ
      Return
      End


*Deck EvFZet
      Subroutine EvFZet(S,Zeta,FZeta,dfZdz,d2fZdz,d3fZdz)
      Implicit Real*8(A-H,O-Z)
C
C     Evaluate f(Zeta) and its derivatives for DFT.
C
      Real*8 MDCutO, Nine
      Save Zero, One, Two, Three, Four, Nine, F8, F27
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/,
     $  Four/4.0d0/, Nine/9.0d0/, F8/8.0D0/, F27/27.0D0/
C
C      Small = MDCutO(0)
      Small = 1.0D-14
      FZeta = -Two
      dfZdz = Zero
      d2fZdz = Zero
      d3fZdz = Zero
      OMZ = One - Zeta
      OPZ = One + Zeta
      OMZ2 = OMZ**2
      OPZ2 = OPZ**2
      F13 = One / Three
      F43 = Four / Three
      F49 = Four / Nine
      F827 = F8 / F27
      If(OMZ.gt.Small) then
        OMZ3 = OMZ ** F13
        fZeta = fZeta + OMZ*OMZ3
        dfZdz = dfZdz - OMZ3
        d2fZdz = d2fZdz + OMZ3/OMZ
        d3fZdz = d3fZdz + OMZ3/OMZ2
        endIf
      If(OPZ.gt.Small) then
        OPZ3 = OPZ ** F13
        fZeta = fZeta + OPZ*OPZ3
        dfZdz = dfZdz + OPZ3
        d2fZdz = d2fZdz + OPZ3/OPZ
        d3fZdz = d3fZdz - OPZ3/OPZ2
        endIf
      fZeta = fZeta * S
      dfZdz = dfZdz * F43 * S
      d2fZdz = d2fZdz * F49 * S
      d3fZdz = d3fZdz * F827 * S
      Return
      End


*Deck EvPWLC
      Subroutine EvPWLC(A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS,d3VdRS)
      Implicit Real*8(A-H,O-Z)
C
C     Evaluate the interpolation function for PW91 local correlation.
C
      Save F1, F2, F3, F4
      Data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/
C
      Q0 = -F2*A*(F1+A1*RS)
      RS12 = Sqrt(RS)
      RS32 = RS*RS12
      Q1 = F2*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RS)
      Q2 = Log(F1+F1/Q1)
      V = Q0*Q2
C
      dQ0dRS = -F2*A*A1
      dQ1dRS = A*(B1/RS12+F2*B2+F3*B3*RS12+F4*B4*RS)
      dQ2dRS = -dQ1dRS/(Q1+Q1**2)
      dVdRS = dQ0dRS*Q2 + Q0*dQ2dRS
C
      d2Q1dS = A*(F3*B3/(RS12*F2)-B1/(RS32*F2)+F4*B4)
      d2Q2dS = (F2*Q1+F1)*(dQ1dRS/(Q1+Q1**2))**2 - d2Q1dS/(Q1+Q1**2)
      d2VdRS = F2*dQ0dRS*dQ2dRS + Q0*d2Q2dS
C
      d3Q1dS = A*F3*(B1/RS-B3)/(RS32*F4)
      d3Q2dS = -F2*(F3*(Q1+F1)*Q1+F1)*(dQ1dRS/(Q1+Q1**2))**3
     $ + F3*(F2*Q1+F1)*dQ1dRS*d2Q1dS/(Q1+Q1**2)**2
     $ - d3Q1dS/(Q1+Q1**2)
      d3VdRS = F3*dQ0dRS*d2Q2dS + Q0*d3Q2dS
      Return
      End


*Deck pbeH0
      Subroutine pbeH0(MaxDrF,Rho,GRho,RS,Zeta,PotLC,dLdS,dLdz,dLdSS,
     $  dLdSZ,dLdZZ,H0,d0dR,d0dG,d0dZ,d0dRR,d0dRG,d0dGG,d0dRZ,d0dGZ,
     $  d0dZZ)
      Implicit Real*8(A-H,O-Z)
C
C     Non-local part of the PBE correlation functional.
C     This routine is a simple modification of the original pw91 routine
C
C     Matthias Ernzerhof
C
      Real*8 MDCutO
      Save XNu,CC0,Pt5,F0,F1,F2,F3,F4,F6,F7,F8,F12
      Data XNu/15.75592D0/,CC0/0.004235D0/,Pt5/0.5d0/,F1/1.0d0/,
     $  F2/2.0d0/,F3/3.0d0/,F4/4.0d0/,F6/6.0d0/,F7/7.0d0/,F8/8.0d0/,
     $  F12/12.0d0/,F0/0.0d0/
C
C      Small = MDCutO(0)
      Small = 1.0D-14
      F13 = F1 / F3
      F16 = Pt5 * F13
      F23 = F13 + F13
      F43 = F23 + F23
      F76 = F7 / F6
      Beta = XNu*CC0
      Pi = F4*ATan(F1)
      GPBE = (F1-Log(F2))/Pi**2
      Alpha =  Beta**2/(F2*GPBE)
      Delta = F2*Alpha/Beta
      BOvD = Beta/Delta
      FKFac = (F3*Pi*Pi)**F13
      SKFac = Sqrt(FKFac/Pi) * F2
C
      Rho16 = Rho**F16
      Rho76 = Rho*Rho16
      OPZ = F1 + Zeta
      OMZ = F1 - Zeta
      If(OPZ.gt.Small) then
        OPZ23 = OPZ**F23
        G = Pt5*OPZ23
        dGdZ = G*F23/OPZ
        dGdZZ = -dGdZ*F13/OPZ
      else
        G = F0
        dGdZ = F0
        dGdZZ = F0
        endIf
      If(OMZ.gt.Small) then
        OMZ23 = OMZ**F23
        GM = Pt5*OMZ23
        G = G + GM
        dGdZM =  -F23*GM/OMZ
        dGdZ = dGdZ + dGdZM
        dGdZZ = dGdZZ + F13*dGdZM/OMZ
        endIf
      T = GRho/(F2*SKFac*Rho76*G)
      G2 = G*G
      G3 = G*G2
      PON = -Delta*PotLC/(G3*Beta)
      ExpPON = Exp(PON)
      B = Delta/(ExpPON-F1)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      Q4 = F1 + B*T2
      Q5 = Q4 + B2*T4
      Ar = F1 + Delta*Q4*T2/Q5
      H0X = Log(Ar)
      H0 = G3*BOvD*H0X
C
      If(MaxDrF.eq.0) Return
      dRsdR = (-(F13*Rs/Rho))
      dPodR = (dLdS*dRsdR*PON/PotLC)
      dPodZ = (-F3*dGdZ*PON/G + dLdZ*PON/PotLC)
      dBdR = (B*dPodR*ExpPON/(F1 - ExpPON))
      dBdZ = (B*dPodZ*ExpPON/(F1 - ExpPON))
      dT2dR = (-F2*F76*T2/Rho)
      dT2dG = (F2*T2/GRho)
      dT2dZ = (-F2*dGdZ*T2/G)
      dQ4dR = (B*dT2dR + dBdR*T2)
      dQ4dG = (B*dT2dG)
      dQ4dZ = (B*dT2dZ + dBdZ*T2)
      dQ5dR = (dQ4dR + F2*B**2*dT2dR*T2 + F2*B*dBdR*T2**2)
      dQ5dG = (dQ4dG + F2*B**2*dT2dG*T2)
      dQ5dZ = (dQ4dZ + F2*B**2*dT2dZ*T2 + F2*B*dBdZ*T2**2)
      dArdR = (Delta*(dT2dR*Q4*Q5 - dQ5dR*Q4*T2 + dQ4dR*Q5*T2)/Q5**2)
      dArdG = (Delta*(dT2dG*Q4*Q5 - dQ5dG*Q4*T2 + dQ4dG*Q5*T2)/Q5**2)
      dArdZ = (Delta*(dT2dZ*Q4*Q5 - dQ5dZ*Q4*T2 + dQ4dZ*Q5*T2)/Q5**2)
      d0dR = (BOvD*dArdR*G**3/Ar)
      d0dG = (BOvD*dArdG*G**3/Ar)
      d0dZ = (BOvD*G**2*(dArdZ*G + F3*Ar*dGdZ*Log(Ar))/Ar)
C
      If(MaxDrF.le.1) Return
      dRsdRR = (-F43*dRsdR/Rho)
      dPodRR = ((dLdSS*dRsdR**2 + dLdS*dRsdRR)*PON/PotLC)
      dPodRZ = ((-F3*dGdZ*dLdS*dRsdR + dLdSZ*dRsdR*G)*PON/(G*PotLC))
      dPodZZ = (PON*(-F6*dGdZ*dLdZ*G + dLdZZ*G**2 + F12*dGdZ**2*PotLC -
     $  F3*dGdZZ*G*PotLC)/(G**2*PotLC))
      dBdRR = (B*ExpPON*(dPodR**2 + dPodRR + dPodR**2*ExpPON -
     $  dPodRR*ExpPON)/(-F1 + ExpPON)**2)
      dBdRZ = (B*ExpPON*(dPodRZ + dPodR*dPodZ - dPodRZ*ExpPON +
     $  dPodR*dPodZ*ExpPON)/(-F1 + ExpPON)**2)
      dBdZZ = (B*ExpPON*(dPodZ**2 + dPodZZ + dPodZ**2*ExpPON -
     $  dPodZZ*ExpPON)/(-F1 + ExpPON)**2)
      dT2dRR = (F2*F76*(F1 + F2*F76)*T2/Rho**2)
      dT2dRG = (-F4*F76*T2/(GRho*Rho))
      dT2dRZ = (F4*dGdZ*F76*T2/(G*Rho))
      dT2dGG = (F2*T2/GRho**2)
      dT2dGZ = (-F4*dGdZ*T2/(G*GRho))
      dT2dZZ = (F2*(F3*dGdZ**2 - dGdZZ*G)*T2/G**2)
      dQ4dRR = (F2*dBdR*dT2dR + B*dT2dRR + dBdRR*T2)
      dQ4dRG = (dBdR*dT2dG + B*dT2dRG)
      dQ4dRZ = (dBdZ*dT2dR + B*dT2dRZ + dBdR*dT2dZ + dBdRZ*T2)
      dQ4dGG = (B*dT2dGG)
      dQ4dGZ = (dBdZ*dT2dG + B*dT2dGZ)
      dQ4dZZ = (F2*dBdZ*dT2dZ + B*dT2dZZ + dBdZZ*T2)
      dQ5dRR = (dQ4dRR + F2*B**2*dT2dR**2 + F8*B*dBdR*dT2dR*T2 +
     $  F2*B**2*dT2dRR*T2 + F2*dBdR**2*T2**2 + F2*B*dBdRR*T2**2)
      dQ5dRG = (dQ4dRG + F2*B**2*dT2dG*dT2dR + F4*B*dBdR*dT2dG*T2 +
     $  F2*B**2*dT2dRG*T2)
      dQ5dRZ = (dQ4dRZ + F2*B**2*dT2dR*dT2dZ + F4*B*dBdZ*dT2dR*T2 +
     $  F2*B**2*dT2dRZ*T2 + F4*B*dBdR*dT2dZ*T2 + F2*B*dBdRZ*T2**2 +
     $  F2*dBdR*dBdZ*T2**2)
      dQ5dGG = (dQ4dGG + F2*B**2*dT2dG**2 + F2*B**2*dT2dGG*T2)
      dQ5dGZ = (dQ4dGZ + F2*B**2*dT2dG*dT2dZ + F4*B*dBdZ*dT2dG*T2 +
     $  F2*B**2*dT2dGZ*T2)
      dQ5dZZ = (dQ4dZZ + F2*B**2*dT2dZ**2 + F8*B*dBdZ*dT2dZ*T2 +
     $  F2*B**2*dT2dZZ*T2 + F2*dBdZ**2*T2**2 + F2*B*dBdZZ*T2**2)
      dQ5dR = dQ5dR/Q5
      dQ5dG = dQ5dG/Q5
      dQ5dZ = dQ5dZ/Q5
      dArdRR = Delta*(-F2*dQ5dR*dT2dR*Q4 + F2*dQ4dR*dT2dR +
     $  dT2dRR*Q4 + F2*dQ5dR**2*Q4*T2 - F2*dQ4dR*dQ5dR*T2 -
     $  (dQ5dRR/Q5)*Q4*T2 + dQ4dRR*T2)/Q5
      dArdRG = Delta*(-dQ5dR*dT2dG*Q4 - dQ5dG*dT2dR*Q4 +
     $  dQ4dR*dT2dG + dQ4dG*dT2dR + dT2dRG*Q4 +
     $  F2*dQ5dG*dQ5dR*Q4*T2 - dQ4dR*dQ5dG*T2 - dQ4dG*dQ5dR*T2 -
     $  (dQ5dRG/Q5)*Q4*T2 + dQ4dRG*T2)/Q5
      dArdRZ = Delta*(-dQ5dZ*dT2dR*Q4 - dQ5dR*dT2dZ*Q4 +
     $  dQ4dZ*dT2dR + dQ4dR*dT2dZ + dT2dRZ*Q4 +
     $  F2*dQ5dR*dQ5dZ*Q4*T2 - dQ4dZ*dQ5dR*T2 - dQ4dR*dQ5dZ*T2 -
     $  (dQ5dRZ/Q5)*Q4*T2 + dQ4dRZ*T2)/Q5
      dArdGG = Delta*(-F2*dQ5dG*dT2dG*Q4 + F2*dQ4dG*dT2dG +
     $  dT2dGG*Q4 + F2*dQ5dG**2*Q4*T2 - F2*dQ4dG*dQ5dG*T2 -
     $  (dQ5dGG/Q5)*Q4*T2 + dQ4dGG*T2)/Q5
      dArdGZ = Delta*(-dQ5dZ*dT2dG*Q4 - dQ5dG*dT2dZ*Q4 +
     $  dQ4dZ*dT2dG + dQ4dG*dT2dZ + dT2dGZ*Q4 +
     $  F2*dQ5dG*dQ5dZ*Q4*T2 - dQ4dZ*dQ5dG*T2 - dQ4dG*dQ5dZ*T2 -
     $  (dQ5dGZ/Q5)*Q4*T2 + dQ4dGZ*T2)/Q5
      dArdZZ = Delta*(-F2*dQ5dZ*dT2dZ*Q4 + F2*dQ4dZ*dT2dZ +
     $  dT2dZZ*Q4 + F2*dQ5dZ**2*Q4*T2 - F2*dQ4dZ*dQ5dZ*T2 -
     $  (dQ5dZZ/Q5)*Q4*T2 + dQ4dZZ*T2)/Q5
      d0dRR = (BOvD*(-dArdR**2 + Ar*dArdRR)*G**3/Ar**2)
      d0dRG = (BOvD*(-(dArdG*dArdR) + Ar*dArdRG)*G**3/Ar**2)
      d0dRZ = (BOvD*G**2*(F3*Ar*dArdR*dGdZ + Ar*dArdRZ*G -
     $  dArdR*dArdZ*G)/Ar**2)
      d0dGG = (BOvD*(-dArdG**2 + Ar*dArdGG)*G**3/Ar**2)
      d0dGZ = (BOvD*G**2*(F3*Ar*dArdG*dGdZ + Ar*dArdGZ*G -
     $  dArdG*dArdZ*G)/Ar**2)
      d0dZZ = (BOvD*G*(F6*Ar*dArdZ*dGdZ*G - dArdZ**2*G**2 +
     $  Ar*dArdZZ*G**2 + F6*Ar**2*dGdZ**2*Log(Ar) +
     $  F3*Ar**2*dGdZZ*G*Log(Ar))/Ar**2)
      Return
      End







