% PEM Electrolyzer
% Parameters
%Vrev = 1.229; % (V) Reversible voltaje
T=25; % (°C) Electrolyzer operating temperature
Pe=101000; % (pa)
Pe_out = 50000000; % (Pa) state of maximun charge in the tank
eta_c = 0.8; % compressor's isentropic efficiency
%Vr = 2.0; % (V) rated voltage
%Ir = 80; %(A) rated current
A = 0.25; % (m^2) area of electrode
Nc = 12; % Number of cells conected in series
F = 96485.34; % (C/mol) Faraday's constant
ne = 2; % 
DG_25 = 237000;
DG_80 = 228480;
DG = DG_25-(T-25)/55*(DG_25-DG_80);
Vrev = DG/ne/F;
DH = 286000;
Vtn = DH/ne/F; % thermoneutral voltage
min_charge = 0.2*Pe_out; % minumun state of charge in the hydrogen tank
max_charge = 0.95*Pe_out; % maximun state of charge in the hydrogen tank
soc_i=0.19*Pe_out; % initial state of charge
n_i = soc_i*0.3/8.3145/(25+273.15);
R = 8.31445; % (J/mol-K) universal constant of gases
Ne = 50; % number of electrolyzers
Vtank = 0.3; % (m^3) volume of the tank

r1 = 7.331e-5; %(ohm m^2) ri parameter for ohmic resistance of electrolyte
r2 = -1.107e-7; % (ohm m2 °C^-1)
r3 = 0; 
s1 = 1.586e-1; % (V) si and ti parameters for overvoltage on electrodes
s2 = 1.378e-3; % (V°C^-1)
s3 = -1.606e-5; % V °C^-2)
t1 = 1.599e-2;
t2 = -1.302;
t3 = 4.213e2;

Iinitial = 0;
Ifinal = 870;
Istep = 1;
nn = (Ifinal-Iinitial)/Istep;

I = [0:1:870];

Vgraph=zeros(nn+1,1);
Id=zeros(nn+1,1);
 for i = 1:nn+1
%V(i) = Vrev+(r1+r2*T)*I(i)/A+(s1+s2*T+s3*T^2)*log((t1+t2/T+t3/T^2)*I(i)/A+1);
Vgraph(i) = Vrev + (r1+r2*T)*I(i)./A+(s1+s2*T+s3*T^2)*log10((t1+t2/T+t3/T^2)*I(i)./A+1);
%Id(i) = I(i)*1000/A*0.0001;
Id(i) = I(i)/2500;
 end

%plot(Id,Vgraph)
 
tvector = importdata('time.txt');
Pvector = importdata ('power.txt');
X=tvector;
Y=Pvector;
p=polyfit(X,Y,21);
timespan = X(length(X));

Pl=zeros(timespan-1,1);
Pr=zeros(timespan-1,1);
Ptot=zeros(timespan-1,1);
Ir=zeros(timespan-1,1);
P=zeros(timespan-1,1);
V=zeros(timespan-1,1);
ne=zeros(timespan-1,1);
Qh2_V=zeros(timespan-1,1);
Qh2_m=zeros(timespan-1,1);
m_dotH2=zeros(timespan-1,1);
Ptank=zeros(timespan-1,1);
Tout=zeros(timespan-1,1);
W_c=zeros(timespan-1,1);
moles=zeros(timespan,1);
soc = soc_i;
moles(1,1)=n_i;

for i = 1:timespan-1
t(i)=i;
% Power profile (kW)
%Pl1(i)=p(1)*t(i)^21+p(2)*t(i)^20+p(3)*t(i)^19+p(4)*t(i)^18+p(5)*t(i)^17+p(6)*t(i)^16+p(7)*t(i)^15+p(8)*t(i)^14+p(9)*t(i)^13+p(10)*t(i)^12+p(11)*t(i)^11+p(12)*t(i)^10+p(13)*t(i)^9+p(14)*t(i)^8+p(15)*t(i)^7+p(16)*t(i)^6+p(17)*t(i)^5+p(18)*t(i)^4+p(19)*t(i)^3+p(20)*t(i)^2+p(21)*t(i)+p(22);
Pl(i) = 1.1*(-1.66886501583441e-91*t(i)^21+1.47641487070072e-85*t(i)^20-6.04737306005779e-80*t(i)^19+1.5215691572903e-74*t(i)^18-2.63102213044733e-69*t(i)^17+3.31408128540438e-64*t(i)^16-3.14500220264774e-59*t(i)^15+2.29376348189984e-54*t(i)^14-1.29999214364872e-49*t(i)^13+5.75117958605101e-45*t(i)^12-1.98427144802997e-40*t(i)^11+5.30725470343592e-36*t(i)^10-1.08802897081482e-31*t(i)^9+1.67983622356027e-27*t(i)^8-1.90337588931417e-23*t(i)^7+1.52346864091771e-19*t(i)^6-8.11379421971959e-16*t(i)^5+2.58376941897956e-12*t(i)^4-3.77390895140414e-09*t(i)^3-1.85238743353635e-07*t(i)^2-0.00331412861521689*t(i)+645.405712039294)/10;
Pr(i) = 1000*Pl(i)/Ne; % (W) to power one stack of electrolyzer

Pri=Pr(i);

Vinitial=round(Vrev,2);
Vfinal=2.2;
stepV=0.001;
mm=round((Vfinal-Vinitial)/stepV);

co=0;
while co ~= 1
for j = 1:mm
Va=round(Vrev,2)+ stepV*(j-1);
B = Vrev + (r1+r2*T)*Pri/A/Nc/Va+(s1+s2*T+s3*T^2)*log10((t1+t2/T+t3/T^2)*Pri/A/Nc/Va+1);
diff = abs(B-Va);
%V(i) = Vrev + (r1+r2*T)*Ir(i)./A+(s1+s2*T+s3*T^2)*log10((t1+t2/T+t3/T^2)*Ir(i)./A+1);
if diff < 0.005
co=1;
V(i)=B;
end
end
end

%while (soc <= min_charge && soc <= max_charge)

Ir(i) = Pr(i)/Nc/V(i);
P(i) = Ne*Nc*V(i)*Ir(i);

 % Faraday Efficiency
 a1 = 0.995; % 99.5%
 a2 = -9.5788; % (m^2*A^-1)
 a3 = -0.0555; % (m^2*A^-1*°C)
 a4 = 0;
 a5 = 1502.7083; % (m^4*A^-1)
 a6 = -70.8005; % (m^4*A^-1*°C-1)
 a7 = 0;
 nf=a1*exp((a2+a3*T+a4*T^2)/(Ir(i)/A)+(a5+a6*T+a7*T^2)/(Ir(i)/A)^2);
 
 % energy (or voltaje) efficiency of a cell
 ne(i)=Vtn/V(i);
 
 % flow of H2 produced
 Qh2_V(i) = 80.69*Nc*Ir(i)*nf/2/F; %(Nm^3/h)
 Qh2_m(i) = Ne*Nc*Ir(i)*nf/2/F; % (mol/s)
 m_dotH2(i)=Qh2_m(i)*0.001; % (kg/s)
 
%  %compressor model
  gamma = 1.41;
  cpH2=14.31; % kj/kg-K
  Ptank(i)=moles(i)*R*(T+273.15)/Vtank;
  Tout(i)=(T+273.15)*(Ptank(i)/Pe)^((gamma-1)/gamma);
  W_c(i)=(m_dotH2(i)/eta_c)*cpH2*(Tout(i)-(T+273.15)); % (kW)
  Ptot(i) = W_c(i)*1000 + P(i);
%  % Storage tank
  
  moles(i+1) = moles(i) + Qh2_m(i)*1; % number of moles in time i in the tank
  %soc_i = 0.19*Pe_out 

  soc(i) = Ptank(i);
  if soc(i)>= max_charge
  disp('Tank is fully charge:')
  disp(['charging_time: ', num2str(i),' seconds']) 
  Soc = round(soc(i)/max_charge*100,1);
  disp(['State of charge ',num2str(Soc),'%'])
  break
  end
end

%   plot(1:i,W_c(1:i));
%   plot(1:i,Ptank(1:i));
%   plot(1:i,Ir(1:i));
%   plot(1:i,V(1:i));
%   plot(1:i,moles(1:i));
%   plot(m_dotH2(1:i)*1000); % (g/s)
%   plot(1:i,Ptot(1:i));
%   plot(1:i,ne(1:i));
%   plotyy(1:i,Ptot(1:i),1:i,Pr(1:i))




