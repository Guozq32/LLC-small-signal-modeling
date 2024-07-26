%   Description:
% This program is used to calculate the small-signal model of the LLC
% converter based on the time-domain analysis for PO mode. The models for PFM
% TSC are deduced and the calculation of the quiescent-state operating point 
% is a prerequisite for the calculation of this program

% Input voltage Vin is 60V.
% Resonant capacitor Cr is 365nF
% Resonant inductor Lr is 24uH
% Magnetizing inductor of the transformer Lm is 60uH
% Resonant frequency fr is 53.7kHz
% Load resistance R is 40Ω
% turns ratio of the transformer n is 1
% Switching frequency is 43kHz
% Output capacitor Co is 36.2uF
clc;
clear
s=tf('s');
Lr=24e-6; Cr=365e-9; Lm=60e-6; Co=36.2e-6; R=40; vin=60; n=1; ts=23.25581e-06;
Ln=Lm/Lr;
Z0=sqrt(Lr/Cr);wr0=1/sqrt(Lr*Cr);
Z1=sqrt((Lr+Lm)/Cr);wr1=1/sqrt((Lr+Lm)*Cr);
fs=1/ts;

In=vin/Z0;
%%%%%Setting the Quiescent Point value%%%%%%%

vo= 81.5794;
ir0=-7.0514;
vr0=-44.1701;
ir0N=ir0/In;
vr0N=vr0/vin;
t2= 9.8762e-06;
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

m0s=1/(M/Ln-r0*cos(phi0+theta0));
m0i=m0s*((sin(phi0+theta0)-sin(theta0))*h0i+r0*(cos(phi0+theta0)-cos(theta0))*g0i);
m0v=m0s*((sin(phi0+theta0)-sin(theta0))*h0v+r0*(cos(phi0+theta0)-cos(theta0))*g0v);
m0in=m0s*((sin(phi0+theta0)-sin(theta0))*h0in+r0*(cos(phi0+theta0)-cos(theta0))*g0in+M/vin/Ln);
m0o=m0s*((sin(phi0+theta0)-sin(theta0))*h0o+r0*(cos(phi0+theta0)-cos(theta0))*g0o-n/vin/Ln);

ir2N=r0*sin(phi0+theta0);
vr2N=-r0*cos(phi0+theta0)+(1-M);

k2i=sin(phi0+theta0)*h0i+r0*cos(phi0+theta0)*(g0i+m0i);
k2v=sin(phi0+theta0)*h0v+r0*cos(phi0+theta0)*(g0v+m0v);
k2in=sin(phi0+theta0)*h0in+r0*cos(phi0+theta0)*(g0in+m0in);
k2o=sin(phi0+theta0)*h0o+r0*cos(phi0+theta0)*(g0o+m0o);

l2i=-cos(phi0+theta0)*h0i+r0*sin(phi0+theta0)*(g0i+m0i);
l2v=-cos(phi0+theta0)*h0v+r0*sin(phi0+theta0)*(g0v+m0v);
l2in=-cos(phi0+theta0)*h0in+r0*sin(phi0+theta0)*(g0in+m0in)+M/vin;
l2o=-cos(phi0+theta0)*h0o+r0*sin(phi0+theta0)*(g0o+m0o)-n/vin;

%%%%%%%%%%%%%%% From t2 to t3 %%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1=atan(-sqrt(1+Ln)*ir2N/(vr2N-1));
r1=sqrt((1+Ln)*ir2N^2+(vr2N-1)^2);
phi1=wr1*(ts/2-phi0/wr0);

g1i=-sqrt(1+Ln)*(vr2N-1)/r1^2*k2i+sqrt(1+Ln)*ir2N/r1^2*l2i;
g1v=-sqrt(1+Ln)*(vr2N-1)/r1^2*k2v+sqrt(1+Ln)*ir2N/r1^2*l2v;
g1in=-sqrt(1+Ln)*(vr2N-1)/r1^2*k2in+sqrt(1+Ln)*ir2N/r1^2*l2in;
g1o=-sqrt(1+Ln)*(vr2N-1)/r1^2*k2o+sqrt(1+Ln)*ir2N/r1^2*l2o;

h1i=sqrt(1+Ln)*ir2N/r1*k2i+(vr2N-1)/r1*l2i;
h1v=sqrt(1+Ln)*ir2N/r1*k2v+(vr2N-1)/r1*l2v;
h1in=sqrt(1+Ln)*ir2N/r1*k2in+(vr2N-1)/r1*l2in;
h1o=sqrt(1+Ln)*ir2N/r1*k2o+(vr2N-1)/r1*l2o;

m1t=wr1/2;
m1i=-m0i*wr1/wr0;
m1v=-m0v*wr1/wr0;
m1in=-m0in*wr1/wr0;
m1o=-m0o*wr1/wr0;

ir3N=r1/sqrt(1+Ln)*sin(phi1+theta1);
vr3N=-r1*cos(phi1+theta1)+1;

k3i=sin(phi1+theta1)/sqrt(1+Ln)*h1i+r1*cos(phi1+theta1)/sqrt(1+Ln)*(g1i+m1i);
k3v=sin(phi1+theta1)/sqrt(1+Ln)*h1v+r1*cos(phi1+theta1)/sqrt(1+Ln)*(g1v+m1v);
k3in=sin(phi1+theta1)/sqrt(1+Ln)*h1in+r1*cos(phi1+theta1)/sqrt(1+Ln)*(g1in+m1in);
k3o=sin(phi1+theta1)/sqrt(1+Ln)*h1o+r1*cos(phi1+theta1)/sqrt(1+Ln)*(g1o+m1o);
k3t=r1*cos(phi1+theta1)/sqrt(1+Ln)*m1t;

l3i=-cos(phi1+theta1)*h1i+r1*sin(phi1+theta1)*(g1i+m1i);
l3v=-cos(phi1+theta1)*h1v+r1*sin(phi1+theta1)*(g1v+m1v);
l3in=-cos(phi1+theta1)*h1in+r1*sin(phi1+theta1)*(g1in+m1in);
l3o=-cos(phi1+theta1)*h1o+r1*sin(phi1+theta1)*(g1o+m1o);
l3t=r1*sin(phi1+theta1)*m1t;

%%%%%%%%%%%%%%%%%%%% The average output current of the rectifier from t0 to t3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

krec1i=n*In/(wr0*ts)*((cos(theta0)-cos(phi0+theta0)-sin(theta0)*phi0)*h0i...
    +(-r0*sin(theta0)+r0*sin(phi0+theta0)-r0*cos(theta0)*phi0)*g0i+(r0*sin(phi0+theta0)-r0*sin(theta0)-M/Ln*phi0)*m0i);
krec1v=n*In/(wr0*ts)*((cos(theta0)-cos(phi0+theta0)-sin(theta0)*phi0)*h0v...
    +(-r0*sin(theta0)+r0*sin(phi0+theta0)-r0*cos(theta0)*phi0)*g0v+(r0*sin(phi0+theta0)-r0*sin(theta0)-M/Ln*phi0)*m0v);
krec1in=n*In/(wr0*ts)*((cos(theta0)-cos(phi0+theta0)-sin(theta0)*phi0)*h0in...
    +(-r0*sin(theta0)+r0*sin(phi0+theta0)-r0*cos(theta0)*phi0)*g0in+(r0*sin(phi0+theta0)-r0*sin(theta0)-M/Ln*phi0)*m0in+M/(2*vin*Ln)*phi0^2)...
    +n*Cr/ts*(r0*cos(theta0)-r0*cos(phi0+theta0)-r0*sin(theta0)*phi0-M/(2*Ln)*phi0^2);
krec1o=n*In/(wr0*ts)*((cos(theta0)-cos(phi0+theta0)-sin(theta0)*phi0)*h0o...
    +(-r0*sin(theta0)+r0*sin(phi0+theta0)-r0*cos(theta0)*phi0)*g0o+(r0*sin(phi0+theta0)-r0*sin(theta0)-M/Ln*phi0)*m0o-n/(2*vin*Ln)*phi0^2);
krec1t=-n*In/(wr0*ts^2)*(r0*cos(theta0)-r0*cos(phi0+theta0)-r0*sin(theta0)*phi0-M/(2*Ln)*phi0^2);


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

m2s=-1/(r2*cos(phi2+theta2)+M/Ln);
m2i=m2s*((sin(phi2+theta2)-sin(theta2))*h2i+r2*(cos(phi2+theta2)-cos(theta2))*g2i);
m2v=m2s*((sin(phi2+theta2)-sin(theta2))*h2v+r2*(cos(phi2+theta2)-cos(theta2))*g2v);
m2in=m2s*((sin(phi2+theta2)-sin(theta2))*h2in+r2*(cos(phi2+theta2)-cos(theta2))*g2in-M/(vin*Ln)*phi2);
m2o=m2s*((sin(phi2+theta2)-sin(theta2))*h2o+r2*(cos(phi2+theta2)-cos(theta2))*g2o+n/(vin*Ln)*phi2);
m2t=m2s*((sin(phi2+theta2)-sin(theta2))*h2t+r2*(cos(phi2+theta2)-cos(theta2))*g2t);

ir5N=r2*sin(phi2+theta2);
vr5N=-r2*cos(phi2+theta2)-(1-M);

k5i=sin(phi2+theta2)*h2i+r2*cos(phi2+theta2)*(g2i+m2i);
k5v=sin(phi2+theta2)*h2v+r2*cos(phi2+theta2)*(g2v+m2v);
k5in=sin(phi2+theta2)*h2in+r2*cos(phi2+theta2)*(g2in+m2in);
k5o=sin(phi2+theta2)*h2o+r2*cos(phi2+theta2)*(g2o+m2o);
k5t=sin(phi2+theta2)*h2t+r2*cos(phi2+theta2)*(g2t+m2t);

l5i=-cos(phi2+theta2)*h2i+r2*sin(phi2+theta2)*(g2i+m2i);
l5v=-cos(phi2+theta2)*h2v+r2*sin(phi2+theta2)*(g2v+m2v);
l5in=-cos(phi2+theta2)*h2in+r2*sin(phi2+theta2)*(g2in+m2in)-M/vin;
l5o=-cos(phi2+theta2)*h2o+r2*sin(phi2+theta2)*(g2o+m2o)+n/vin;
l5t=-cos(phi2+theta2)*h2t+r2*sin(phi2+theta2)*(g2t+m2t);

%%%%%%%%%%%%%%%%%%% From t5 to t6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta3=pi+atan(-sqrt(1+Ln)*ir5N/(vr5N+1));
r3=sqrt((1+Ln)*ir5N^2+(vr5N+1)^2);
phi3=wr1*(ts/2-phi2/wr0);

g3i=-sqrt(1+Ln)*(vr5N+1)/r3^2*k5i+sqrt(1+Ln)*ir5N/r3^2*l5i;
g3v=-sqrt(1+Ln)*(vr5N+1)/r3^2*k5v+sqrt(1+Ln)*ir5N/r3^2*l5v;
g3in=-sqrt(1+Ln)*(vr5N+1)/r3^2*k5in+sqrt(1+Ln)*ir5N/r3^2*l5in;
g3o=-sqrt(1+Ln)*(vr5N+1)/r3^2*k5o+sqrt(1+Ln)*ir5N/r3^2*l5o;
g3t=-sqrt(1+Ln)*(vr5N+1)/r3^2*k5t+sqrt(1+Ln)*ir5N/r3^2*l5t;

h3i=sqrt(1+Ln)*ir5N/r3*k5i+(vr5N+1)/r3*l5i;
h3v=sqrt(1+Ln)*ir5N/r3*k5v+(vr5N+1)/r3*l5v;
h3in=sqrt(1+Ln)*ir5N/r3*k5in+(vr5N+1)/r3*l5in;
h3o=sqrt(1+Ln)*ir5N/r3*k5o+(vr5N+1)/r3*l5o;
h3t=sqrt(1+Ln)*ir5N/r3*k5t+(vr5N+1)/r3*l5t;

m3t=(wr1/2-wr1/wr0*m2t);
m3i=-wr1/wr0*m2i;
m3v=-wr1/wr0*m2v;
m3in=-wr1/wr0*m2in;
m3o=-wr1/wr0*m2o;

k6i=sin(phi3+theta3)/sqrt(1+Ln)*h3i+r3*cos(phi3+theta3)/sqrt(1+Ln)*(g3i+m3i);
k6v=sin(phi3+theta3)/sqrt(1+Ln)*h3v+r3*cos(phi3+theta3)/sqrt(1+Ln)*(g3v+m3v);
k6in=sin(phi3+theta3)/sqrt(1+Ln)*h3in+r3*cos(phi3+theta3)/sqrt(1+Ln)*(g3in+m3in);
k6o=sin(phi3+theta3)/sqrt(1+Ln)*h3o+r3*cos(phi3+theta3)/sqrt(1+Ln)*(g3o+m3o);
k6t=sin(phi3+theta3)/sqrt(1+Ln)*h3t+r3*cos(phi3+theta3)/sqrt(1+Ln)*(g3t+m3t);

l6i=-cos(phi3+theta3)*h3i+r3*sin(phi3+theta3)*(g3i+m3i);
l6v=-cos(phi3+theta3)*h3v+r3*sin(phi3+theta3)*(g3v+m3v);
l6in=-cos(phi3+theta3)*h3in+r3*sin(phi3+theta3)*(g3in+m3in);
l6o=-cos(phi3+theta3)*h3o+r3*sin(phi3+theta3)*(g3o+m3o);
l6t=-cos(phi3+theta3)*h3t+r3*sin(phi3+theta3)*(g3t+m3t);

%%%%%%%%%%%%%%%%%%%% The average output current of the rectifier from t3 to t6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
krec2i=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2)-sin(theta2)*phi2)*h2i...
    +(-r2*sin(theta2)+r2*sin(phi2+theta2)-r2*cos(theta2)*phi2)*g2i+(r2*sin(phi2+theta2)-r2*sin(theta2)+M/Ln*phi2)*m2i);
krec2v=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2)-sin(theta2)*phi2)*h2v...
    +(-r2*sin(theta2)+r2*sin(phi2+theta2)-r2*cos(theta2)*phi2)*g2v+(r2*sin(phi2+theta2)-r2*sin(theta2)+M/Ln*phi2)*m2v);
krec2in=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2)-sin(theta2)*phi2)*h2in...
    +(-r2*sin(theta2)+r2*sin(phi2+theta2)-r2*cos(theta2)*phi2)*g2in+(r2*sin(phi2+theta2)-r2*sin(theta2)+M/Ln*phi2)*m2in-M/(2*vin*Ln)*phi2^2)...
    +n*Cr/ts*(r2*cos(theta2)-r2*cos(phi2+theta2)-r2*sin(theta2)*phi2+M/(2*Ln)*phi0^2);
krec2o=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2)-sin(theta2)*phi2)*h2o...
    +(-r2*sin(theta2)+r2*sin(phi2+theta2)-r2*cos(theta2)*phi2)*g2o+(r2*sin(phi2+theta2)-r2*sin(theta2)+M/Ln*phi2)*m2o+n/(2*vin*Ln)*phi2^2);
krec2t=n*In/(wr0*ts)*((cos(theta2)-cos(phi2+theta2)-sin(theta2)*phi2)*h2t...
    +(-r2*sin(theta2)+r2*sin(phi2+theta2)-r2*cos(theta2)*phi2)*g2t+(r2*sin(phi2+theta2)-r2*sin(theta2)+M/Ln*phi2)*m2t)...
    -n*In/(wr0*ts^2)*(r2*cos(theta2)-r2*cos(phi2+theta2)-r2*sin(theta2)*phi2+M/(2*Ln)*phi0^2);

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

Az=1/(wr0+g2t)*[-g0i-g2i -g0v-g2v -g0o-g2o];
Bz=1/(wr0+g2t)*[-g0in-g2in 2*wr0];
Ac=A+B2*Az;
Bc=[(B1+B2*(-g0in-g2in)/(wr0+g2t)) 2*wr0/(wr0+g2t)*B2];
Gcs=C*(s*eye(3)-Ac)^-1*Bc;% 对控制时间tcs的电流控制传函
Gtcs=exp(-ts/2*s)*Gcs(1,2);      %%%%%%% The transfer function from switching period to output voltage 

figure(2),bode(Gtcs,p);grid on; hold on
xlim([100,1e5]);
G=exp(-ts/2*s)*C*(s*eye(3)-A)^-1*B; %输出扰动对周期传函

