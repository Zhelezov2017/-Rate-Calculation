function men = RateCalculation(k0, q, rho, EE, GG, HH )

p = sqrt(1 - q.^2);
q = q.* (2*(imag(q) <= 0)-1);
Q = k0.* rho * q;

    mainq = EE.^2 - GG.^2 + EE.* HH - (HH + EE).* p.^2;
    radq = sqrt((HH - EE).^2 * p.^4 + 2 * ((GG.^2).* (HH + EE) - EE.* (HH - EE).^2) * p.^2 +...
        (EE.^2 - GG.^2 - EE.* HH).^2);
    q1 = sqrt(0.5 * (mainq - radq)./ EE);
    q2 = sqrt(0.5 * (mainq + radq)./ EE);

    n1 = -(EE.* (p.* GG).^(-1)).* (p.^2 + q1.^2 + (GG.^2)./ EE - EE);
    n2 = -(EE.* (p * GG).^(-1)).* (p.^2 + q2.^2 + (GG.^2)./ EE - EE);
    
    Q1 = (q1.* rho).* k0;
    Q2 = (q2.* rho).* k0;
    
    D1 = q1.* k0;
    D2 = q2.* k0;
    
    Jm0_1 = besselj(0, Q1);
    Jm1_1 = besselj(1, Q1);
    Jm2_1 = besselj(2, Q1);
    
    Jm0_2 = besselj(0, Q2);
    Jm1_2 = besselj(1, Q2);
    Jm2_2 = besselj(2, Q2);
    
    H1m  = besselh(0, 1, Q);
    dH1m = (H1m.* 0)./ Q  - besselh(0 + 1, 1, Q);
    
    Int_1 = rho* (D1 * Jm2_1 * Jm1_2 - D2* Jm1_1*Jm2_2) /(D1^2 - D2^2);
    Int_2 = rho^2*(Jm1_1^2- Jm0_1* Jm2_1)/2;
    Int_3 = rho^2*(Jm1_2^2- Jm0_2* Jm2_2)/2;
    
    Int_4=rho^2*dH1m^2/2 -1/2 *(rho^2+1/(k0.* q)^2) * H1m^2;
    
    %константы
    Bs1 = 0.709459854502189 - 3.13930306809721i;
    B_s1 = -0.709459854502189 - 3.13930306809721i;
    Bs2 = 1;
    B_s2 = 6.98387166058172e-18 + 0.00169038056806685i;
    
    
    men = 4 * pi * (EE ^ (-1) * Int_2 * ( p * n1^2 * Bs1 * Bs1 + GG * Bs1 * Bs1 * n1) - Int_2 * (p * Bs1 * Bs1) + EE^(-1) * Int_3 * ( p * n2^2 * Bs2 * Bs2 + GG * Bs2 * Bs2 * n2) -...
    Int_3 * ( p * Bs2 * Bs2) + EE^(-1)* Int_1 * ( p * Bs2 * Bs1 * n1 * n2 + p * Bs2 * Bs1 * n2 * n1 + GG * n2 * Bs2 * Bs1 + GG * n1 * Bs2 * Bs1) - Int_1 * (p * (Bs2 * Bs1 + Bs2 * Bs1))- Int_4 *  (p*B_s2^2 + B_s1*B_s2-  p*B_s1^2));





end

