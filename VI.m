 function F = VI(x,Pri,Nc,Vrev,T,r1,r2,s1,s2,s3,t1,t2,t3,A)
 Vr = x(1);
 Ir = x(2);
 fa = Vr - Vrev - (r1+r2*T)*Ir./A - (s1+s2*T+s3*T^2)*log10((t1+t2/T+t3/T^2)*Ir./A+1);
 fb = (Nc*Vr)*Ir-Pri;
 F=[fa;fb];
 end