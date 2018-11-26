function ron = RateCalculation2(k0, q, rho, EE, GG, HH )

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
    
  
    Km  = -2*(besselk(1, Q *1i))/pi;
    dKm = Km./ Q   + 2*besselk( 2, Q * 1i)/pi;
    
    Jm_outH1  = besselj(1, Q);
    Jm_outH0  = besselj(0, Q);
    Jm_outH2  = besselj(2, Q);
    dJm_outH = Jm_outH1.* 1 / Q  - besselj(2, Q);
    
    Ym_outH0 = bessely(0, Q);
    Ym_outH1 = bessely(1, Q);
    Ym_outH2 = bessely(2, Q);
    
    
    
    
    Jm_out1  = besselj(1, Q1);
    Jm_out12  = besselj(2, Q1);
    dJm_out1 = Jm_out1.* 1 / Q1  - besselj(2, Q1);
    
    Jm_out2  = besselj(1, Q2);
    Jm_out22  = besselj(2, Q2);
    dJm_out2 = Jm_out2.* 1 / Q2  - besselj(2, Q2);
    
    Int_1 = rho* (D1 * Jm_out12 * Jm_out2 - D2* Jm_out1*Jm_out22) /(D1^2 - D2^2);
    Int_2 = (rho^2*dJm_out1^2)/2 + ((rho^2-1/ (k0.* q1)^2)* Jm_out1^2)/2;
    Int_3 = (rho^2*dJm_out2^2)/2 + ((rho^2-1/ (k0.* q2)^2)* Jm_out2^2)/2;
    
    
    %Int_4 = (rho^2*dKm^2)/2 - ((rho^2+1/ (k0.* q)^2)* Km^2)/2;;
    Int_6 = rho^2*(Ym_outH1^2-Ym_outH0*Ym_outH2)/2;
    Int_5 = rho^2*(2*Jm_outH1*Ym_outH1-Jm_outH2 *Ym_outH0-Ym_outH2*Jm_outH0)/4;
    Int_4=rho^2*(Jm_outH1^2-Jm_outH0*Jm_outH2)/2;
    
    
    Int_7=Int_4-2i*Int_5-Int_6;
    %константы
    Bs1 = 0.709459854502189 - 3.13930306809721i;
    B_s1 = -0.709459854502189 - 3.13930306809721i;
    Bs2 = 1;
    B_s2 = 6.98387166058172e-18 + 0.00169038056806685i;
    
    
    %men = 4 * pi * (EE ^ (-1) * Int_2 * ( p * n1^2 * Bs1 * Bs1 + GG * Bs1 * Bs1 * n1) - Int_2 * (p * Bs1 * Bs1) + EE^(-1) * Int_3 * ( p * n2^2 * Bs2 * Bs2 + GG * Bs2 * Bs2 * n2) -...
    %Int_3 * ( p * Bs2 * Bs2) + EE^(-1)* Int_1 * ( p * Bs2 * Bs1 * n1 * n2 + p * Bs2 * Bs1 * n2 * n1 + GG * n2 * Bs2 * Bs1 + GG * n1 * Bs2 * Bs1) - Int_1 * (p * (Bs2 * Bs1 + Bs2 * Bs1))- Int_4 *  (p*B_s2^2 + B_s1*B_s2-  p*B_s1^2));
    
    men1 = 4*pi *EE ^ (-1)  * ( p * n1^2 * Bs1 * Bs1 + GG * Bs1 * Bs1 * n1) * Int_2;
     
    ron = -4*pi * Int_7 *  (p*B_s2^2 + B_s1*B_s2-  p*B_s1^2);





end