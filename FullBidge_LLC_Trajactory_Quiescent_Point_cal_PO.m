%   Description:
% This program is used to calculate the quiescent point of the LLC
% converter based on the time-domain analysis for PO mode. Given
% input voltage, resonance parameters, switching period, and load,
% the output voltage can be calculated. In addition, the time-domain 
% expressions at quisient state can also be obtained.

% Input voltage Vin is 60V.
% Resonant capacitor Cr is 365nF
% Resonant inductor Lr is 24uH
% Magnetizing inductor of the transformer Lm is 60uH
% Resonant frequency fr is 53.7kHz
% Load resistance R is 40Ω
% turns ratio of the transformer n is 1
% Switching frequency is 43kHz
% Output capacitor C is 36.2uF
clc;
clear
Lr=24e-6;Cr=365e-9;Lm=60e-6;C=36.2e-6;R=40;Vin=60;n=1; Ts=23.25581e-06;
Ln=Lm/Lr;
Z0=sqrt(Lr/Cr);wr0=1/sqrt(Lr*Cr);
Z1=sqrt((Lr+Lm)/Cr);wr1=1/sqrt((Lr+Lm)*Cr);

In=Vin/Z0;

%%%%% Calculate the output voltage through the switching cycle  %%%%%
%%%%% Setting the Quiescent Point initial value for iteration %%%%%%%
Vo=60; 
Ir0=-7;
Vcr0=-45;
Ir0N=Ir0/In;
Vcr0N=Vcr0/Vin;
Ir2=5;
Vcr2=40;
Ir2N=Ir2/In;
Vcr2N=Vcr2/Vin;

t2=10e-6;
tr1=12e-7;
M=n*Vo/Vin;

eps1=1e-7;
eps2=1e-7;

lemda=0.005;

k=0;


r0=sqrt(Ir0N^2+(Vcr0N-(1-M))^2);
the0=atan(-Ir0N/(Vcr0N-(1-M)));
r1=sqrt((1+Ln)*Ir2N^2+(Vcr2N-1)^2);
the1=atan(-sqrt(1+Ln)*Ir2N/(Vcr2N-1));
phi0=t2*wr0;
phi1=tr1*wr1;                  

x=[r0,phi0,the0,r1,phi1,the1,M]';

%%%%%%%%%%%%%%%%%%%%%%
for k=1:5000
 b=F(x);
 norm_b=norm(b,inf);
 if norm_b<eps1
     break;
 end
 A=Jac(x);
 dx=-inv(A)*b;
 x=x+lemda*dx;
 norm_dx=norm(dx,inf);
 if norm_dx<eps2
     break;
 end
end
k
r0=x(1)
phi0=x(2)
the0=x(3)
r1=x(4)
phi1=x(5)
the1=x(6)
M=x(7)


Ir0N=r0*sin(the0)
Vcr0N=-r0*cos(the0)+(1-M)
Ir3N=r1/(sqrt(1+Ln))*sin(phi1+the1)
Vcr3N=-r1*cos(phi1+the1)+1
Ir2N=r0*sin(phi0+the0)
Vcr2N=-r0*cos(phi0+the0)+1-M
t2=phi0/wr0
Vo=M*Vin/n

Ir0=In*Ir0N
Vcr0=Vin*Vcr0N
Ir2=In*Ir2N
Vcr2=Vin*Vcr2N
b

Jac(x)

    
function y=F(x)
y=zeros(size(x));
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);

Lr=24e-6;Cr=365e-9;Lm=60e-6;C=36.2e-6;R=40;Vin=60;n=1; Ts=23.25581e-06;
Ln=Lm/Lr;
Z0=sqrt(Lr/Cr);wr0=1/sqrt(Lr*Cr);
Z1=sqrt((Lr+Lm)/Cr);wr1=1/sqrt((Lr+Lm)*Cr);
In=Vin/Z0;
IoN=2*n^2*R*In/(Ts*wr0*Vin);

y(1)=r0*sin(phi0+the0)-r1/(sqrt(1+Ln))*sin(the1);
y(2)=-r0*cos(phi0+the0)-M+r1*cos(the1);
y(3)=r1/(sqrt(1+Ln))*sin(phi1+the1)+r0*sin(the0);
y(4)=-r1*cos(phi1+the1)-r0*cos(the0)+2-M;
y(5)=r0*sin(phi0+the0)-r0*sin(the0)-M/Ln*phi0;
y(6)=M-IoN*(r0*cos(the0)-r0*cos(phi0+the0)-r0*sin(the0)*phi0-M/(2*Ln)*phi0^2);
y(7)=phi0/wr0+phi1/wr1-Ts/2;

end

%%% Jacobian matrix %%%%%
function Mx=Jac(x)
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);


Lr=24e-6;Cr=365e-9;Lm=60e-6;C=36.2e-6;R=40;Vin=60;n=1; Ts=23.25581e-06;
Ln=Lm/Lr;
Z0=sqrt(Lr/Cr);wr0=1/sqrt(Lr*Cr);
Z1=sqrt((Lr+Lm)/Cr);wr1=1/sqrt((Lr+Lm)*Cr);
In=Vin/Z0;
IoN=2*n^2*R*In/(Ts*wr0*Vin);



Mx=[sin(phi0+the0) r0*cos(phi0+the0) r0*cos(phi0+the0) -1/sqrt(1+Ln)*sin(the1) 0 r1/sqrt(1+Ln)*cos(the1) 0;...
    -cos(phi0+the0) r0*sin(phi0+the0) r0*sin(phi0+the0) cos(the1) 0 -r1*sin(the1) -1;...
    sin(the0) 0 r0*cos(the0) 1/sqrt(1+Ln)*sin(phi1+the1) r1/sqrt(1+Ln)*cos(phi1+the1) r1/sqrt(1+Ln)*cos(phi1+the1) 0;...
    -cos(the0) 0 r0*sin(the0) -cos(phi1+the1) r1*sin(phi1+the1) r1*sin(phi1+the1) -1;...
    sin(phi0+the0)-sin(the0) r0*cos(phi0+the0)-M/Ln r0*cos(phi0+the0)-r0*cos(the0) 0 0 0 -phi0/Ln;...
    -IoN*(cos(the0)-cos(phi0+the0)-sin(the0)*phi0) -IoN*(r0*sin(phi0+the0)-r0*sin(the0)-M/Ln*phi0) -IoN*(-r0*sin(the0)+r0*sin(phi0+the0)-r0*cos(the0)*phi0) 0 0 0 1-IoN*(-1/(2*Ln)*phi0^2);...
    0 1/wr0 0 0 1/wr1 0 0]; 

end












