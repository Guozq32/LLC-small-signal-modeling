%   Description:
% This program is used to calculate the small-signal model of the LLC
% converter based on the time-domain analysis for NP mode.The models for PFM
% TSC are deduced and the calculation of the quiescent-state operating point 
% is a prerequisite for the calculation of this program

% Input voltage Vin is 60V.
% Resonant capacitor Cr is 365nF
% Resonant inductor Lr is 24uH
% Magnetizing inductor of the transformer Lm is 60uH
% Resonant frequency fr is 53.7kHz
% Load resistance R is 40Ω
% turns ratio of the transformer n is 1
% Switching frequency is 65kHz
% Output capacitor Co is 36.2uF
clear
clc
s=tf('s');
Lr=24e-6; Cr=365e-9; Lm=60e-6; Co=36.2e-6; R=40; vin=60; n=1; ts=15.384615e-06;
Ln=Lm/Lr;
Z0=sqrt(Lr/Cr);wr0=1/sqrt(Lr*Cr);
Z1=sqrt((Lr+Lm)/Cr);wr1=1/sqrt((Lr+Lm)*Cr);
fs=1/ts;

In=vin/Z0;
%%%%%Setting the Quiescent Point value%%%%%%%

vo=51.392131;  
ir0= -3.29436729;
vr0= -13.53849728;
ir0N=ir0/In;
vr0N=vr0/vin;
t2= 7.504764e-6;
phi0=wr0*t2;
M=n*vo/vin;

%%%%%%%%%%% From t0 to t2 %%%%%%%%%%%%%
theta0=atan(-ir0N/(vr0N-(1-M)));
r0=sqrt(ir0N^2+(vr0N-(1-M))^2);


g0i=-(vr0N-(1-M))/r0^2;
g0v=ir0N/r0^2;
g0in=-ir0N*M/vin/r0^2;
g0o=n*ir0N/vin/r0^2;

h0i=ir0N/r0;
h0v=(vr0N-(1-M))/r0;
h0in=-(vr0N-(1-M))*M/vin/r0;
h0o=(vr0N-(1-M))*n/vin/r0;

ir2N=r0*sin(phi0+theta0);
vr2N=-r0*cos(phi0+theta0)+(1-M);

k2i=sin(phi0+theta0)*h0i+r0*cos(phi0+theta0)*g0i;
k2v=sin(phi0+theta0)*h0v+r0*cos(phi0+theta0)*g0v;
k2in=sin(phi0+theta0)*h0in+r0*cos(phi0+theta0)*g0in;
k2o=sin(phi0+theta0)*h0o+r0*cos(phi0+theta0)*g0o;
k2m0=r0*cos(phi0+theta0);

l2i=-cos(phi0+theta0)*h0i+r0*sin(phi0+theta0)*g0i;
l2v=-cos(phi0+theta0)*h0v+r0*sin(phi0+theta0)*g0v;
l2in=-cos(phi0+theta0)*h0in+r0*sin(phi0+theta0)*g0in+M/vin;
l2o=-cos(phi0+theta0)*h0o+r0*sin(phi0+theta0)*g0o-n/vin;
l2m0=r0*sin(phi0+theta0);

%%%%%%%%%%%%%%% From t2 to t3 %%%%%%%%%%%%%%%%%%%%%%%%%%%

theta1=pi+atan(-ir2N/(vr2N+1+M));
r1=sqrt(ir2N^2+(vr2N+1+M)^2);
phi1=wr0*(ts/2-phi0/wr0);

g1i=-(vr2N+1+M)/r1^2*k2i+ir2N/r1^2*l2i;
g1v=-(vr2N+1+M)/r1^2*k2v+ir2N/r1^2*l2v;
g1in=-(vr2N+1+M)/r1^2*k2in+ir2N/r1^2*l2in-M*ir2N/vin/r1^2;
g1o=-(vr2N+1+M)/r1^2*k2o+ir2N/r1^2*l2o+n*ir2N/vin/r1^2;
g1m0=-(vr2N+1+M)/r1^2*k2m0+ir2N/r1^2*l2m0;

h1i=ir2N/r1*k2i+(vr2N+1+M)/r1*l2i;
h1v=ir2N/r1*k2v+(vr2N+1+M)/r1*l2v;
h1in=ir2N/r1*k2in+(vr2N+1+M)/r1*l2in-M*(vr2N+1+M)/vin/r1;
h1o=ir2N/r1*k2o+(vr2N+1+M)/r1*l2o+n*(vr2N+1+M)/vin/r1;
h1m0=ir2N/r1*k2m0+(vr2N+1+M)/r1*l2m0;

m1=-1/(r1*cos(phi1+theta1)-sin(phi1+theta1)*h1m0-r1*cos(phi1+theta1)*g1m0);
m1i=m1*(-sin(theta0)*h0i+sin(phi1+theta1)*h1i-r0*cos(theta0)*g0i+r1*cos(phi1+theta1)*g1i);
m1v=m1*(-sin(theta0)*h0v+sin(phi1+theta1)*h1v-r0*cos(theta0)*g0v+r1*cos(phi1+theta1)*g1v);
m1in=m1*(-sin(theta0)*h0in+sin(phi1+theta1)*h1in-r0*cos(theta0)*g0in+r1*cos(phi1+theta1)*g1in+wr0*ts*M/vin/2/Ln);
m1o=m1*(-sin(theta0)*h0o+sin(phi1+theta1)*h1o-r0*cos(theta0)*g0o+r1*cos(phi1+theta1)*g1o-wr0*ts*n/vin/2/Ln);
m1t=m1*wr0/2*(sin(phi1+theta1)*h1m0+r1*cos(phi1+theta1)*g1m0-M/Ln);

ir3N=r1*sin(phi1+theta1);
vr3N=-r1*cos(phi1+theta1)-1-M;

k3i=(sin(phi1+theta1)*h1i+r1*cos(phi1+theta1)*g1i+(r1*cos(phi1+theta1)-sin(phi1+theta1)*h1m0-r1*cos(phi1+theta1)*g1m0)*m1i);
k3v=(sin(phi1+theta1)*h1v+r1*cos(phi1+theta1)*g1v+(r1*cos(phi1+theta1)-sin(phi1+theta1)*h1m0-r1*cos(phi1+theta1)*g1m0)*m1v);
k3in=(sin(phi1+theta1)*h1in+r1*cos(phi1+theta1)*g1in+(r1*cos(phi1+theta1)-sin(phi1+theta1)*h1m0-r1*cos(phi1+theta1)*g1m0)*m1in);
k3o=(sin(phi1+theta1)*h1o+r1*cos(phi1+theta1)*g1o+(r1*cos(phi1+theta1)-sin(phi1+theta1)*h1m0-r1*cos(phi1+theta1)*g1m0)*m1o);
k3t=(sin(phi1+theta1)*h1m0*wr0/2+r1*cos(phi1+theta1)*g1m0*wr0/2+(r1*cos(phi1+theta1)-sin(phi1+theta1)*h1m0-r1*cos(phi1+theta1)*g1m0)*m1t);

l3i=(-cos(phi1+theta1)*h1i+r1*sin(phi1+theta1)*g1i+(r1*sin(phi1+theta1)+cos(phi1+theta1)*h1m0-r1*sin(phi1+theta1)*g1m0)*m1i);
l3v=(-cos(phi1+theta1)*h1v+r1*sin(phi1+theta1)*g1v+(r1*sin(phi1+theta1)+cos(phi1+theta1)*h1m0-r1*sin(phi1+theta1)*g1m0)*m1v);
l3in=(-cos(phi1+theta1)*h1in+r1*sin(phi1+theta1)*g1in+(r1*sin(phi1+theta1)+cos(phi1+theta1)*h1m0-r1*sin(phi1+theta1)*g1m0)*m1in+M/vin);
l3o=(-cos(phi1+theta1)*h1o+r1*sin(phi1+theta1)*g1o+(r1*sin(phi1+theta1)+cos(phi1+theta1)*h1m0-r1*sin(phi1+theta1)*g1m0)*m1o-n/vin);
l3t=(-cos(phi1+theta1)*h1m0*wr0/2+r1*sin(phi1+theta1)*g1m0*wr0/2+(r1*sin(phi1+theta1)+cos(phi1+theta1)*h1m0-r1*sin(phi1+theta1)*g1m0)*m1t);

%%%%%%%%%%%%%%%%%%%% The average output current of the rectifier from t0 to t3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

krec1i=n*In/(wr0*ts)*((cos(theta0)-cos(phi0+theta0))*h0i+(cos(theta1)-cos(phi1+theta1))*h1i+(-r0*sin(theta0)+r0*sin(phi0+theta0))*g0i+(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1i...
    +(-(cos(theta1)-cos(phi1+theta1))*h1m0-(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1m0-r0*sin(phi0+theta0)+r1*sin(phi1+theta1))*m1i);

krec1v=n*In/(wr0*ts)*((cos(theta0)-cos(phi0+theta0))*h0v+(cos(theta1)-cos(phi1+theta1))*h1v+(-r0*sin(theta0)+r0*sin(phi0+theta0))*g0v+(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1v...
    +(-(cos(theta1)-cos(phi1+theta1))*h1m0-(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1m0-r0*sin(phi0+theta0)+r1*sin(phi1+theta1))*m1v);

krec1in=n*In/(wr0*ts)*((cos(theta0)-cos(phi0+theta0))*h0in+(cos(theta1)-cos(phi1+theta1))*h1in+(-r0*sin(theta0)+r0*sin(phi0+theta0))*g0in+(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1in...
    +(-(cos(theta1)-cos(phi1+theta1))*h1m0-(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1m0-r0*sin(phi0+theta0)+r1*sin(phi1+theta1))*m1in)+n*Cr/ts*(r0*cos(theta0)-r0*cos(phi0+theta0)+r1*cos(theta1)-r1*cos(phi1+theta1));

krec1o=n*In/(wr0*ts)*((cos(theta0)-cos(phi0+theta0))*h0o+(cos(theta1)-cos(phi1+theta1))*h1o+(-r0*sin(theta0)+r0*sin(phi0+theta0))*g0o+(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1o...
    +(-(cos(theta1)-cos(phi1+theta1))*h1m0-(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1m0-r0*sin(phi0+theta0)+r1*sin(phi1+theta1))*m1o);

krec1t=n*In/(wr0*ts)*((cos(theta1)-cos(phi1+theta1))*h1m0*wr0/2+(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1m0*wr0/2-1/ts*(r0*cos(theta0)-r0*cos(phi0+theta0)+r1*cos(theta1)-r1*cos(phi1+theta1))...
    +(-(cos(theta1)-cos(phi1+theta1))*h1m0-(-r1*sin(theta1)+r1*sin(phi1+theta1))*g1m0-r0*sin(phi0+theta0)+r1*sin(phi1+theta1))*m1t+r0*sin(phi0+theta0)*wr0/2);

%%%%%%%%%%%%%%%%%% From t3 to t5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta2=pi+atan(-ir3N/(vr3N+(1-M)));
r2=sqrt(ir3N^2+(vr3N+(1-M))^2);
phi2=phi0;

g2i=-(vr3N+(1-M))/r2^2*k3i+ir3N/r2^2*l3i;
g2v=-(vr3N+(1-M))/r2^2*k3v+ir3N/r2^2*l3v;
g2in=-(vr3N+(1-M))/r2^2*k3in+ir3N/r2^2*l3in+ir3N*M/vin/r2^2;
g2o=-(vr3N+(1-M))/r2^2*k3o+ir3N/r2^2*l3o-n*ir3N/vin/r2^2;
g2t=-(vr3N+(1-M))/r2^2*k3t+ir3N/r2^2*l3t;

h2i=ir3N/r2*k3i+(vr3N+(1-M))/r2*l3i;
h2v=ir3N/r2*k3v+(vr3N+(1-M))/r2*l3v;
h2in=ir3N/r2*k3in+(vr3N+(1-M))/r2*l3in+(vr3N+(1-M))*M/vin/r2;
h2o=ir3N/r2*k3o+(vr3N+(1-M))/r2*l3o-(vr3N+(1-M))*n/vin/r2;
h2t=ir3N/r2*k3t+(vr3N+(1-M))/r2*l3t;

ir5N=r2*sin(phi2+theta2);
vr5N=-r2*cos(phi2+theta2)-(1-M);

k5i=sin(phi2+theta2)*h2i+r2*cos(phi2+theta2)*g2i;
k5v=sin(phi2+theta2)*h2v+r2*cos(phi2+theta2)*g2v;
k5in=sin(phi2+theta2)*h2in+r2*cos(phi2+theta2)*g2in;
k5o=sin(phi2+theta2)*h2o+r2*cos(phi2+theta2)*g2o;
k5t=sin(phi2+theta2)*h2t+r2*cos(phi2+theta2)*g2t;
k5m2=r2*cos(phi2+theta2);

l5i=-cos(phi2+theta2)*h2i+r2*sin(phi2+theta2)*g2i;
l5v=-cos(phi2+theta2)*h2v+r2*sin(phi2+theta2)*g2v;
l5in=-cos(phi2+theta2)*h2in+r2*sin(phi2+theta2)*g2in-M/vin;
l5o=-cos(phi2+theta2)*h2o+r2*sin(phi2+theta2)*g2o+n/vin;
l5t=-cos(phi2+theta2)*h2t+r2*sin(phi2+theta2)*g2t;
l5m2=r2*sin(phi2+theta2);

%%%%%%%%%%%%%%%%%%% From t5 to t6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta3=atan(-ir5N/(vr5N-1-M));
r3=sqrt(ir5N^2+(vr5N-1-M)^2);
phi3=wr0*(ts/2-phi2/wr0);

g3i=-(vr5N-1-M)/r3^2*k5i+ir5N/r3^2*l5i;
g3v=-(vr5N-1-M)/r3^2*k5v+ir5N/r3^2*l5v;
g3in=-(vr5N-1-M)/r3^2*k5in+ir5N/r3^2*l5in+M*ir5N/vin/r3^2;
g3o=-(vr5N-1-M)/r3^2*k5o+ir5N/r3^2*l5o-n*ir5N/vin/r3^2;
g3t=-(vr5N-1-M)/r3^2*k5t+ir5N/r3^2*l5t;
g3m2=-(vr5N-1-M)/r3^2*k5m2+ir5N/r3^2*l5m2;

h3i=ir5N/r3*k5i+(vr5N-1-M)/r3*l5i;
h3v=ir5N/r3*k5v+(vr5N-1-M)/r3*l5v;
h3in=ir5N/r3*k5in+(vr5N-1-M)/r3*l5in+M*(vr5N-1-M)/vin/r3;
h3o=ir5N/r3*k5o+(vr5N-1-M)/r3*l5o-n*(vr5N-1-M)/vin/r3;
h3t=ir5N/r3*k5t+(vr5N-1-M)/r3*l5t;
h3m2=ir5N/r3*k5m2+(vr5N-1-M)/r3*l5m2;

m3=-1/(r3*cos(phi3+theta3)-sin(phi3+theta3)*h3m2-r3*cos(phi3+theta3)*g3m2);
m3i=m3*(-sin(theta2)*h2i+sin(phi3+theta3)*h3i-r2*cos(theta2)*g2i+r3*cos(phi3+theta3)*g3i);
m3v=m3*(-sin(theta2)*h2v+sin(phi3+theta3)*h3v-r2*cos(theta2)*g2v+r3*cos(phi3+theta3)*g3v);
m3in=m3*(-sin(theta2)*h2in+sin(phi3+theta3)*h3in-r2*cos(theta2)*g2in+r3*cos(phi3+theta3)*g3in-wr0*ts*M/vin/2/Ln);
m3o=m3*(-sin(theta2)*h2o+sin(phi3+theta3)*h3o-r2*cos(theta2)*g2o+r3*cos(phi3+theta3)*g3o+wr0*ts*n/vin/2/Ln);
m3t=m3*wr0/2*(sin(phi3+theta3)*h3m2+r3*cos(phi3+theta3)*g3m2+M/Ln);

k6i=(sin(phi3+theta3)*h3i+r3*cos(phi3+theta3)*g3i+(r3*cos(phi3+theta3)-sin(phi3+theta3)*h3m2-r3*cos(phi3+theta3)*g3m2)*m3i);
k6v=(sin(phi3+theta3)*h3v+r3*cos(phi3+theta3)*g3v+(r3*cos(phi3+theta3)-sin(phi3+theta3)*h3m2-r3*cos(phi3+theta3)*g3m2)*m3v);
k6in=(sin(phi3+theta3)*h3in+r3*cos(phi3+theta3)*g3in+(r3*cos(phi3+theta3)-sin(phi3+theta3)*h3m2-r3*cos(phi3+theta3)*g3m2)*m3in);
k6o=(sin(phi3+theta3)*h3o+r3*cos(phi3+theta3)*g3o+(r3*cos(phi3+theta3)-sin(phi3+theta3)*h3m2-r3*cos(phi3+theta3)*g3m2)*m3o);
k6t=(sin(phi3+theta3)*(h3t+h3m2*wr0/2)+r3*cos(phi3+theta3)*(g3t+g3m2*wr0/2)+(r3*cos(phi3+theta3)-sin(phi3+theta3)*h3m2-r3*cos(phi3+theta3)*g3m2)*m3t);

l6i=-cos(phi3+theta3)*h3i+r3*sin(phi3+theta3)*g3i+(r3*sin(phi3+theta3)+cos(phi3+theta3)*h3m2-r3*sin(phi3+theta3)*g3m2)*m3i;
l6v=-cos(phi3+theta3)*h3v+r3*sin(phi3+theta3)*g3v+(r3*sin(phi3+theta3)+cos(phi3+theta3)*h3m2-r3*sin(phi3+theta3)*g3m2)*m3v;
l6in=-cos(phi3+theta3)*h3in+r3*sin(phi3+theta3)*g3in+(r3*sin(phi3+theta3)+cos(phi3+theta3)*h3m2-r3*sin(phi3+theta3)*g3m2)*m3in-M/vin;
l6o=-cos(phi3+theta3)*h3o+r3*sin(phi3+theta3)*g3o+(r3*sin(phi3+theta3)+cos(phi3+theta3)*h3m2-r3*sin(phi3+theta3)*g3m2)*m3o+n/vin;
l6t=-cos(phi3+theta3)*(h3t+wr0*h3m2/2)+r3*sin(phi3+theta3)*(g3t+wr0*g3m2/2)+(r3*sin(phi3+theta3)+cos(phi3+theta3)*h3m2-r3*sin(phi3+theta3)*g3m2)*m3t;

%%%%%%%%%%%%%%%%%%%% The average output current of the rectifier from t3 to t6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
krec2i=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2))*h2i+(cos(theta3)-cos(phi3+theta3))*h3i+(-r2*sin(theta2)+r2*sin(phi2+theta2))*g2i+(-r3*sin(theta3)+r3*sin(phi3+theta3))*g3i...
    +(-(cos(theta3)-cos(phi3+theta3))*h3m2-(-r3*sin(theta3)+r3*sin(phi3+theta3))*g3m2-r2*sin(phi2+theta2)+r3*sin(phi3+theta3))*m3i);

krec2v=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2))*h2v+(cos(theta3)-cos(phi3+theta3))*h3v+(-r2*sin(theta2)+r2*sin(phi2+theta2))*g2v+(-r3*sin(theta3)+r3*sin(phi3+theta3))*g3v...
    +(-(cos(theta3)-cos(phi3+theta3))*h3m2-(-r3*sin(theta3)+r3*sin(phi3+theta3))*g3m2-r2*sin(phi2+theta2)+r3*sin(phi3+theta3))*m3v);

krec2in=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2))*h2in+(cos(theta3)-cos(phi3+theta3))*h3in+(-r2*sin(theta2)+r2*sin(phi2+theta2))*g2in+(-r3*sin(theta3)+r3*sin(phi3+theta3))*g3in...
    +(-(cos(theta3)-cos(phi3+theta3))*h3m2-(-r3*sin(theta3)+r3*sin(phi3+theta3))*g3m2-r2*sin(phi2+theta2)+r3*sin(phi3+theta3))*m3in)+n*Cr/ts*(r2*cos(theta2)-r2*cos(phi2+theta2)+r3*cos(theta3)-r3*cos(phi3+theta3));

krec2o=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2))*h2o+(cos(theta3)-cos(phi3+theta3))*h3o+(-r2*sin(theta2)+r2*sin(phi2+theta2))*g2o+(-r3*sin(theta3)+r3*sin(phi3+theta3))*g3o...
    +(-(cos(theta3)-cos(phi3+theta3))*h3m2-(-r3*sin(theta3)+r3*sin(phi3+theta3))*g3m2-r2*sin(phi2+theta2)+r3*sin(phi3+theta3))*m3o);

krec2t=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2))*h2t+(cos(theta3)-cos(phi3+theta3))*(h3t+h3m2*wr0/2)+(-r2*sin(theta2)+r2*sin(phi2+theta2))*g2t+(-r3*sin(theta3)+r3*sin(phi3+theta3))*(g3t+g3m2*wr0/2)...
    +(-(cos(theta3)-cos(phi3+theta3))*h3m2-(-r3*sin(theta3)+r3*sin(phi3+theta3))*g3m2-r2*sin(phi2+theta2)+r3*sin(phi3+theta3))*m3t-1/ts*(r2*cos(theta2)-r2*cos(phi2+theta2)+r3*cos(theta3)-r3*cos(phi3+theta3))+r2*sin(phi2+theta2)*wr0/2);

%%%%%%%%%%%%%%%%%%%% The small signal model for PFM %%%%%%%%%%%%%%%%%%%%
A=[(k6i-1)/ts k6v/ts k6o/ts; l6i/ts (l6v-1)/ts l6o/ts; (krec1i-krec2i)/Co (krec1v-krec2v)/Co (krec1o-krec2o-1/R)/Co;];
B=[k6in/ts k6t/ts; l6in/ts l6t/ts; (krec1in-krec2in)/Co (krec1t-krec2t)/Co;];
C=[0 0 1];
G=C*(s*eye(3)-A)^-1*B
Gts=exp(-ts/2*s)*G(1,2)     %%%%%%% The transfer function from switching period to output voltage 

p=bodeoptions;
p.FreqUnits='Hz';
figure(1),bode(Gts,p);grid on; hold on
xlim([100,1e5]);
%%%%%%%%%%%%%%%%%%%% The small signal model for TSC %%%%%%%%%%%%%%%%%%%%

B1=[k6in/ts;l6in/ts;(krec1in-krec2in)/Co];
B2=[k6t/ts;l6t/ts;(krec1t-krec2t)/Co];
Az=1/(wr0+g2t-m1t-m3t)*[-g0i-g2i+m1i+m3i -g0v-g2v+m1v+m3v -g0o-g2o+m1o+m3o];
Bz=1/(wr0+g2t-m1t-m3t)*[-g0in-g2in+m1in+m3in 2*wr0];
Ac=A+B2*Az;
Bc=[(B1+B2*(-g0in-g2in+m1in+m3in)/(wr0+g2t-m1t-m3t)) 2*wr0/(wr0+g2t-m1t-m3t)*B2];
Gcs=C*(s*eye(3)-Ac)^-1*Bc;
Gtcs=exp(-ts/2*s)*Gcs(1,2);      %%%%%%% The transfer function from switching period to output voltage 

figure(2),bode(Gtcs,p);grid on; hold on
xlim([100,1e5]);
