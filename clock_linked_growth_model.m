% function Model_ABA_art_fin returns Func 
% in our function definition, we are defining all parameters 

%  

function Func = clock_linked_growth_model(t,y,light_hours, params_clock,Temp,Day, mut, params_growth,dusk)
    
    function output = is_day(t, Daylength)
        t1=mod(t, 24-Daylength);
        if t1==0
            output = 1;
        else
            output = 1-heaviside(t1);
        end
    end
% what is tw? 
tw=0.05;

params_clockCell = num2cell(params_clock);
[q1,q2,q3,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22,g23,g24,g25,g26,g27,g28,g29,a,b,c,d,e,f,g,h,i,j,period,Loffset,Lamplitude,A0] = params_clockCell{:};
params_growthCell = num2cell(params_growth);
[pB28, kr22, kr28, pE122, pE128, pE222, pE228, dE, pPE22, pPE28, dP, kPC, dPB, pCL28, pCD, dC, pG, kG, pGP, pGE, pGB, pGH, pHC, mutBox,mutEox, mutPox, mutPko1, mutPko2, mutCox, mutCko1,mutCko2] = params_growthCell{:};
 
L=1;
dusk=light_hours;
L=Loffset+0.5*Lamplitude*((1+tanh((t-period*floor(t/period))/tw))-(1+tanh((t-period*floor(t/period)-dusk)/tw))+(1+tanh((t-period*floor(t/period)-period)/tw)));   

% here we are instantiating the variable Func as an array of 0's 
Func = zeros(40, 1);
% the y we are passing to the model is a vector as described by below 

%Variables
B=y(36);%PHYB
E=y(37);%ELF3
P=y(38);%PIF
C=y(39);%COP1
G=y(40);%Hypocotyl

%Parameters
% t1 is the hour in the day 
t1=mod(t, 24);
%L_growth=is_day(t1,Day);
pB=10.0;
kr=kr22; %0.232 datos de Casal;
pE1=pE122; %adimensional
pE2=pE222;
pP=1.0;%adimensional
pPE=pPE22;
pCL=1.0;%adimensional
mB=1.0;%maximum PHYB value
if Temp==28
    pB=pB28;
    kr=kr28;%0.411 datos Casal
    pE1=pE128;
    pE2=pE228;
    pPE=pPE28;
    pCL=pCL28;
end


%if 'PHYBox' in mut
if any(strcmp(mut, 'PHYBox'))
    mB = mB*mutBox;
end

%if 'ELF3ox' in mut
if any(strcmp(mut, 'ELF3ox'))
    pE1 = pE1*mutEox;
end

%if 'PIF4ox' in mut
if any(strcmp(mut, 'PIF4ox'))
    pP = pP*mutPox;
end

%if 'pif4' in mut
if any(strcmp(mut, 'pif4'))
    pP = pP*mutPko1;
end

%if 'pifq' in mut
if any(strcmp(mut, 'pifq'))
    pP = pP*mutPko2;
end

%if 'COP1' in mut
if any(strcmp(mut, 'COP1'))
    pCL= pCL*mutCox;
    pCD= pCD*mutCox;
end

%if 'cop1-4' in mut
 if any(strcmp(mut, 'cop1-4'))
    pCL= pCL*mutCko1;
    pCD= pCD*mutCko1;
end

%if 'cop1-6' in mut
if any(strcmp(mut, 'cop1-6'))
    pCL= pCL*mutCko2;
    pCD= pCD*mutCko2;
end

%if 'hy5' in mut
if any(strcmp(mut, 'hy5'))
    pGH=0;
end

%ODEs
% y(1)  LHY mRNA                   
% y(2)  P
% y(3)  GI-ZTL 
% y(4)  GI-ELF3 cytoplasm
% y(5)  LHY prot
% y(6)  TOC1 mRNA
% y(7)  PRR9 prot
% y(8)  PRR5 (NI) mRNA
% y(9)  PRR5 (NI) prot
% y(10) GI prot cytoplasm
% y(11) TOC1 prot
% y(12) ZTL
% y(13) EC
% y(14) GI mRNA
% y(15) PRR9 mRNA
% y(16) PRR7 mRNA
% y(17) PRR7 prot
% y(18) ELF4 mRNA
% y(19) ELF4 prot
% y(20) LHY prot modif.
% y(21) ABAR mRNA
% y(22) COP1 cytoplasm
% y(23) ELF3 mRNA
% y(24) ELF3 cytoplasm
% y(25) ELF3 nuclear
% y(26) COP1 nuclear night
% y(27) COP1 nuclear day
% y(28) LUX mRNA
% y(29) LUX prot
% y(30) PP2C prot
% y(31) SnRK2 prot
% y(32) stomata
% y(33) 
% y(34) 
% y(35) 
% y(36) phyB
% y(37) ELF3
% y(38) PIF proteins 
% y(39) COP1 light and dark proteins
% y(40) growth 

Gn=p28*y(10)/(p29+m19+p17*y(25));
EGn=(p18*y(4)+p17*y(25)*Gn)/(m10*y(26)+m9*y(27)+p31);
e34=p25*y(19)*y(25)/(p26*y(29)+p21+m10*y(26)+m9*y(27));
ar=A0*y(21)/(y(21)+g29);
ar=0.5*(A0+y(21)+g29-sqrt((A0+y(21)+g29)^2-4*A0*y(21)));

% here we define each element of the array previously instantiated. in this
% case, each element of the array corresponds to one of the ODE's
% describing the system 

Func(1) = 1*(q1*L*y(2)+n1*g1^a/(g1^a+(y(7)+y(17)+y(9)+y(11))^a))-y(1)*(m1*L+m2*(1-L));
Func(2) = p7*(1-L)*(1-y(2))-m11*y(2)*L;
Func(3) = p12*L*y(12)*y(10)-p13*y(3)*(1-L)-m21*y(3);
Func(4) = p17*y(24)*y(10)-m10*y(4)*y(22)-p18*y(4)+p31*EGn;
Func(5) = (p2+p1*L)*(y(1))-m3*y(5)-p3*y(5)^c/(y(5)^c+g3^c);
Func(6) = 1*n2/(1+(y(5)/(g5*(1+(y(31)/g25)^j)))^e)*g4/(g4+y(13))-y(6)*m5;
Func(7) = p8*y(15)-(m13+m22*(1-L))*y(7);
Func(8) = 1*g23^g/(g23^g+y(11)^g)*(n10*y(20)^e/(g12^e+y(20)^e)+n11*y(17)^b/(g13^b+y(17)^b))-m16*y(8);
Func(9) = p10*y(8)-(m17+m24*(1-L))*y(9);
Func(10)= p11*y(14)-m19*y(10)-p12*L*y(12)*y(10)+p13*y(3)*(1-L)-p17*y(24)*y(10)-p28*y(10)+p29*Gn;
Func(11)= p4*(y(6)+n16)-m8*y(11)-(m6+m7*(1-L))*y(11)*(p5*y(12)+y(3));
Func(11)= p4*y(6)-m8*y(11)-(m6+m7*(1-L))*y(11)*(p5*y(12)+y(3));
Func(12)= 1*p14-m20*y(12)-p12*L*y(12)*y(10)+p13*y(3)*(1-L);
Func(13)= p26*y(29)*e34-m10*y(13)*y(26)-m9*y(13)*y(27)-m32*y(13)*(1+p24*L*(EGn+Gn)^d/(g7^d+(EGn+Gn)^d));
Func(14)= 1*g17^g/(g17^g+y(11)^g)*(q2*L*y(2)+g15^e/(g15^e+y(5)^e)*g14/(g14+y(13))*n12)-y(14)*m18;
Func(15)= 1*g18^g/(g18^g+y(11)^g)*(q3*L*y(2)+g8/(g8+y(13))*(n4+n7*y(5)^e/(y(5)^e+g9^e)))-m12*y(15);
Func(16)= 1*g22^g/(g22^g+y(11)^g)*(n8*(y(5)+y(20))^e/(g10^e+(y(5)+y(20))^e)+n9*y(7)^f/(g11^f+y(7)^f))-m14*y(16);
Func(17)= p9*y(16)-y(17)*(m15+m23*(1-L));
Func(18)= n15*g21^g/(g21^g+y(11)^g)*g6^e/(g6^e+y(5)^e)*g20/(g20+y(13))-y(18)*m34;
Func(19)= p23*y(18)-m35*y(19)-p25*y(25)*y(19)+p21*e34;
Func(20)= p3*y(5)^c/(y(5)^c+g3^c)-m4*y(20);
Func(21)= n17*y(5)^e/(y(5)^e+g28^e)*g24^g/(g24^g+y(11)^g)-m37*y(21);
Func(22)= 1*n5-p6*y(22)-m27*y(22)*(1+p15*L);
Func(23)= 1*n3*g16^e/(g16^e+y(5)^e)-m26*y(23);
Func(24)= p16*y(23)-m9*y(24)*y(22)-p17*y(24)*y(10)-p19*y(24)+p20*y(25);
Func(25)= p19*y(24)-p20*y(25)-m10*y(25)*y(26)-m9*y(25)*y(27)-p25*y(25)*y(19)+p21*e34-p17*y(25)*Gn;
Func(26)= p6*y(22)-n6*L*y(2)*y(26)-n14*y(26)-m27*y(26)*(1+p15*L);
Func(27)= 1*(n14*y(26)+n6*L*y(2)*y(26))-m31*(1+m33*(1-L))*y(27);
Func(28)= n13*g19^g/(g19^g+y(11)^g)*g6^e/(g6^e+y(5)^e)*g2/(g2+y(13))-y(28)*m34;
Func(29)= p27*y(28)-m36*y(29)-p26*y(29)*e34;
Func(30)= p33*g27^h/(ar^h+g27^h)-m39*y(30);
Func(31)= p32-m30*y(31)*y(30);
Func(32)= (n19+n18*L)*(1-y(32))*g26^i/(g26^i+y(31)^i)-m29*y(32);
Func(33)= 0;
Func(34)= p22-m38*y(34)*y(27)-m25*y(34)*y(26);
Func(35)= p30-m28*y(35)*y(26);

%dBdt 
Func(36) = pB*L*(mB-B)-kr*B;
% dEdt 
Func(37) = Func(25);
% dPdt
Func(38) = pP/(1+pPE*E)-dP*P/(1+kPC*C)-dPB*P*B;

%dCdt
Func(39) = Func(26) + Func(27);

%dGdt
Func(40) = pG+kG*pGP*P/(1+pGP*P+pGE*E+pGB*B+pGH/(1+pHC*C));


%if 'elf3-8' in mut
    if any(strcmp(mut, 'elf3-8'))
        Func(37)=0;
    end
    
    
    %if 'phyB' in mut
     if any(strcmp(mut, 'phyB'))
        Func(36)=0;
     end

end
