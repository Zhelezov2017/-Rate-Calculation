function Integral = NumericalIntegralTrapz(k_0, q, rho, EE, GG, HH )
p = sqrt(1 - q.^2);

Q = k_0.* rho * q;
    Qlow =  k_0.*q;

    mainq = EE.^2 - GG.^2 + EE.* HH - (HH + EE).* p.^2;
    radq = sqrt((HH - EE).^2 * p.^4 + 2 * ((GG.^2).* (HH + EE) - EE.* (HH - EE).^2) * p.^2 +...
        (EE.^2 - GG.^2 - EE.* HH).^2);
    q1 = sqrt(0.5 * (mainq - radq)./ EE);
    q2 = sqrt(0.5 * (mainq + radq)./ EE);

    n1 = -(EE.* (p.* GG).^(-1)).* (p.^2 + q1.^2 + (GG.^2)./ EE - EE);
    n2 = -(EE.* (p * GG).^(-1)).* (p.^2 + q2.^2 + (GG.^2)./ EE - EE);
    
    
    %%%«адаем параметры интегрировани€
    n = 100;
    a = rho;  
    b = 1111111111111111111111122222222222211;  
    h = (b-a) / n;  
    x=a:h:b; 
   
    
    D1 = q1.* k_0;
    D2 = q2.* k_0; 
   
    Jm1_1 = besselj(1, D1*x)    
    Jm1_2 = besselj(1, D2*x);
    
    H1m = besselh(1,2,Qlow*x)
    
    %константы
    Bs1 = 0.709459854502189 - 3.13930306809721i;
    B_s1 = -0.709459854502189 - 3.13930306809721i;
    Bs2 = 1;
    B_s2 = 6.98387166058172e-18 + 0.00169038056806685i;
    
    
    %y= 4 * pi * (EE ^ (-1) * Jm1_1.^2 * ( p * n1^2 * Bs1 * Bs1 + GG * Bs1 * Bs1 * n1) - Jm1_1.^2 * (p * Bs1 * Bs1) + EE^(-1) * Jm1_2.^2 * ( p * n2^2 * Bs2 * Bs2 + GG * Bs2 * Bs2 * n2) -...
    %Jm1_2.^2 * ( p * Bs2 * Bs2) + EE^(-1)* Jm1_1.* Jm1_2.* ( p * Bs2 * Bs1 * n1 * n2 + p * Bs2 * Bs1 * n2 * n1 + GG * n2 * Bs2 * Bs1 + GG * n1 * Bs2 * Bs1) - Jm1_1.* Jm1_2 * (p * (Bs2 * Bs1 + Bs2 * Bs1))- H1m.^2 *  (p*B_s2^2 + B_s1*B_s2-  p*B_s1^2));

    y = 4*pi * x *  (p*B_s2^2 + B_s1*B_s2-  p*B_s1^2).* H1m.^2  ;

    Integral = trapz( x, y);




end
