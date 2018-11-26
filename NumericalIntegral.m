function Integral = NumericalIntegral(k0, q, rho, EE, GG, HH )
p = sqrt(1 - q.^2);
q = q.* (2*(imag(q) <= 0)-1);
Q = k0.* rho * q;
    Qlow =  k0.*q;

    mainq = EE.^2 - GG.^2 + EE.* HH - (HH + EE).* p.^2;
    radq = sqrt((HH - EE).^2 * p.^4 + 2 * ((GG.^2).* (HH + EE) - EE.* (HH - EE).^2) * p.^2 +...
        (EE.^2 - GG.^2 - EE.* HH).^2);
    q1 = sqrt(0.5 * (mainq - radq)./ EE);
    q2 = sqrt(0.5 * (mainq + radq)./ EE);

    n1 = -(EE.* (p.* GG).^(-1)).* (p.^2 + q1.^2 + (GG.^2)./ EE - EE);
    n2 = -(EE.* (p * GG).^(-1)).* (p.^2 + q2.^2 + (GG.^2)./ EE - EE);
    
   
    
    D1 = q1.* k0;
    D2 = q2.* k0;
    
    ValueOf1PointBessel=rho/2-rho/(2*sqrt(3));
    ValueOf2PointBessel=rho/2+rho/(2*sqrt(3));
    
    ValueOf1PointHankel=rho/2+rho/(2*sqrt(3));
    ValueOf2PointHankel=rho/2-rho/(2*sqrt(3));
    
    Jm1_11 = besselj(1, D1*ValueOf1PointBessel);
    Jm1_21 = besselj(1, D2*ValueOf1PointBessel);
    Jm1_12 = besselj(1, D1*ValueOf2PointBessel);
    Jm1_22 = besselj(1, D2*ValueOf2PointBessel);

    
    H1m1  = besselh(0, 1, Qlow * ValueOf1PointHankel);
   
    H1m2  = besselh(0, 1, Qlow * ValueOf2PointHankel);
   
    
    %константы
    Bs1 = 0.709459854502189 - 3.13930306809721i;
    B_s1 = -0.709459854502189 - 3.13930306809721i;
    Bs2 = 1;
    B_s2 = 6.98387166058172e-18 + 0.00169038056806685i;
    
    
    f1Bessel =4 * pi * ValueOf1PointBessel *(EE ^ (-1)*   ( p * n1^2 * Bs1 * Bs1 + GG * Bs1 * Bs1 * n1)* Jm1_11^2 - Jm1_11^2 * (p * Bs1 * Bs1) + EE^(-1) *Jm1_21^2 * ( p * n2^2 * Bs2 * Bs2 + GG * Bs2 * Bs2 * n2) -...
    Jm1_21^2  * ( p * Bs2 * Bs2) + EE^(-1)*  Jm1_11*Jm1_21 * ( p * Bs2 * Bs1 * n1 * n2 + p * Bs2 * Bs1 * n2 * n1 + GG * n2 * Bs2 * Bs1 + GG * n1 * Bs2 * Bs1) - Jm1_11*Jm1_21 * (p * (Bs2 * Bs1 + Bs2 * Bs1)));
    
    f2Bessel =4 * pi * ValueOf2PointBessel *(EE ^ (-1)*   ( p * n1^2 * Bs1 * Bs1 + GG * Bs1 * Bs1 * n1)* Jm1_12^2 - Jm1_12^2 * (p * Bs1 * Bs1) + EE^(-1) * Jm1_22^2 * ( p * n2^2 * Bs2 * Bs2 + GG * Bs2 * Bs2 * n2) -...
    Jm1_22^2 * ( p * Bs2 * Bs2) + EE^(-1)* Jm1_12*Jm1_22 * ( p * Bs2 * Bs1 * n1 * n2 + p * Bs2 * Bs1 * n2 * n1 + GG * n2 * Bs2 * Bs1 + GG * n1 * Bs2 * Bs1) - Jm1_12*Jm1_22 * (p * (Bs2 * Bs1 + Bs2 * Bs1)));
    
    f3Hankel= 4 * pi *H1m1^2 *  (p*B_s2^2 + B_s1*B_s2-  p*B_s1^2)*ValueOf1PointHankel;
    f4Hankel= 4 * pi *H1m2^2 *  (p*B_s2^2 + B_s1*B_s2-  p*B_s1^2)*ValueOf2PointHankel;

    %Integral = rho/2*(f1Bessel+f2Bessel)-rho/2*(f3Hankel+f4Hankel);
    Integral = -rho*(f3Hankel+f4Hankel)/2;


end

