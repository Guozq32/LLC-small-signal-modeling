%   Description:
% This program is used to calculate the quiescent point of the LLC
% converter based on the time-domain analysis for NP mode. Given
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
% Switching frequency is 65kHz
% Output capacitor C is 36.2uF
clear
clc
Lr=24e-6;Cr=365e-9;Lm=60e-6;C=36.2e-6;R=40;Vin=60;n=1; Ts=15.384615e-06;
Ln=Lm/Lr;
Z0=sqrt(Lr/Cr);
wr0=1/sqrt(Lr*Cr);

In=Vin/Z0;
%%%%% Calculate the output voltage through the switching cycle  %%%%%
%%%%%Setting the Quiescent_Point initial value for iteration%%%%%%%
Vo=50; 
Ir0=-3;
Vcr0=-14; 
Ir0N=Ir0/In;
Vcr0N=Vcr0/Vin;
Ir2=3;
Vcr2=14;
Ir2N=Ir2/In;
Vcr2N=Vcr2/Vin;

t2=7e-6;
tr1=9e-7;
M=n*Vo/Vin;

eps1=1e-7;
eps2=1e-7;

lemda=0.5;

k=0;

r0=sqrt(Ir0N^2+(Vcr0N-(1-M))^2);
the0=atan(-Ir0N/(Vcr0N-(1-M)));
r1=sqrt(Ir2N^2+(Vcr2N+1+M)^2);
the1=pi+atan(-Ir2N/(Vcr2N+1+M));
phi0=t2*wr0;
phi1=tr1*wr0;

x=[r0,phi0,the0,r1,phi1,the1,M]';

%%%%%%%%%%%%%%%%%%%%%%
for k=1:500
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
Ir3N=r1*sin(phi1+the1)
Vcr3N=-r1*cos(phi1+the1)-1-M
Ir2N=r0*sin(phi0+the0)
Vcr2N=-r0*cos(phi0+the0)+1-M
t2=phi0/wr0
Vo=M*Vin/n

Ir0=In*Ir0N
Vcr0=Vin*Vcr0N
b


    
function y=F(x)
y=zeros(size(x));
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);

Lr=24e-6;Cr=365e-9;Lm=60e-6;C=36.2e-6;R=40;Vin=60;n=1; Ts=15.384615e-06;
Ln=Lm/Lr;
Z0=sqrt(Lr/Cr);wr0=1/sqrt(Lr*Cr);
%Z1=sqrt((Lr+Lm)/Cr);wr1=1/sqrt((Lr+Lm)*Cr);
In=Vin/Z0;
IoN=2*n^2*R*In/(Ts*wr0*Vin);

y(1)=r0*sin(the0)+M*wr0*Ts/(4*Ln);
y(2)=r0*sin(phi0+the0)-r1*sin(the1);
y(3)=-r0*cos(phi0+the0)+r1*cos(the1)+2;
y(4)=r1*sin(phi1+the1)+r0*sin(the0);
y(5)=-r1*cos(phi1+the1)-r0*cos(the0)-2*M;
y(6)=M-IoN*(r0*cos(the0)-r0*cos(phi0+the0)+r1*cos(the1)-r1*cos(phi1+the1));
y(7)=phi0/wr0+phi1/wr0-Ts/2;

end

function Mx=Jac(x)
r0=x(1);
phi0=x(2);
the0=x(3);
r1=x(4);
phi1=x(5);
the1=x(6);
M=x(7);

Lr=24e-6;Cr=365e-9;Lm=60e-6;C=36.2e-6;R=40;Vin=60;n=1; Ts=15.384615e-06;
Ln=Lm/Lr;
Z0=sqrt(Lr/Cr);wr0=1/sqrt(Lr*Cr);
%Z1=sqrt((Lr+Lm)/Cr);wr1=1/sqrt((Lr+Lm)*Cr);
In=Vin/Z0;
IoN=2*n^2*R*In/(Ts*wr0*Vin);

Mx=[sin(the0) 0 r0*cos(the0) 0 0 0 Ts*wr0/(4*Ln);...
    sin(phi0+the0) r0*cos(phi0+the0) r0*cos(phi0+the0) -sin(the1) 0 -r1*cos(the1) 0;...
    -cos(phi0+the0) r0*sin(phi0+the0) r0*sin(phi0+the0) cos(the1) 0 -r1*sin(the1) 0;...
    sin(the0) 0 r0*cos(the0) sin(phi1+the1) r1*cos(phi1+the1) r1*cos(phi1+the1) 0;...
    -cos(the0) 0 r0*sin(the0) -cos(phi1+the1) r1*sin(phi1+the1) r1*sin(phi1+the1) -2;...
    -IoN*(cos(the0)-cos(phi0+the0)) -IoN*r0*sin(phi0+the0) -IoN*(-r0*sin(the0)+r0*sin(phi0+the0)) -IoN*(cos(the1)-cos(phi1+the1)) -IoN*r1*sin(phi1+the1) -IoN*(-r1*sin(the0)+r1*sin(phi1+the1)) 1;...
    0 1/wr0 0 0 1/wr0 0 0]; 






% Mx=[r0I0*sin(angel)+r0*cos(angel)*theI0 r0V0*sin(angel)+r0*cos(angel)*theV0 -1 0 r0*cos(angel) r0Vout*sin(angel)+r0*cos(angel)*theVout;...
%     -r0I0*cos(angel)+r0*sin(angel)*theI0 -r0V0*cos(angel)+r0*sin(angel)*theV0 0 -1 r0*sin(angel) -r0Vout*cos(angel)+r0*sin(angel)*theVout-1;...
%     1 0 r1I2/sqrt(1+Ln)*sin(angel2)+r1/sqrt(1+Ln)*cos(angel2)*theI2 r1V2/sqrt(1+Ln)*sin(angel2)+r1/sqrt(1+Ln)*cos(angel2)*theV2 -r1/sqrt(1+Ln)*cos(angel2)*wr1/wr0 0;...
%     0 1 -r1I2*cos(angel2)+r1*sin(angel2)*theI2 -r1V2*cos(angel2)+r1*sin(angel2)*theV2 -r1*sin(angel2)*wr1/wr0 0;...
%     r0I0*sin(angel)+r0*cos(angel)*theI0-r0I0*sin(theta0)-r0*cos(theta0)*theI0 r0V0*sin(angel)+r0*cos(angel)*theV0-r0V0*sin(theta0)-r0*cos(theta0)*theV0 0 0 r0*cos(angel)-M/Ln r0Vout*sin(angel)+r0*cos(angel)*theVout-r0Vout*sin(theta0)-r0*cos(theta0)*theVout-wr0t2/Ln;...
%     -Nnorm*(r0I0*cos(theta0)/wr0-r0*sin(theta0)*theI0/wr0-r0I0*cos(angel)/wr0+r0*sin(angel)*theI0/wr0-r0I0*sin(theta0)/wr0*wr0t2-r0*cos(theta0)*theI0/wr0*wr0t2) ...
%     -Nnorm*(r0V0*cos(theta0)/wr0-r0*sin(theta0)*theV0/wr0-r0V0*cos(angel)/wr0+r0*sin(angel)*theV0/wr0-r0V0*sin(theta0)/wr0*wr0t2-r0*cos(theta0)*theV0/wr0*wr0t2) ...
%     0 0 -Nnorm*(r0*sin(angel)/wr0-r0*sin(theta0)/wr0-M*wr0t2/(wr0*Ln)) 1-Nnorm*(r0Vout/wr0*cos(theta0)+r0/wr0*sin(theta0)*theVout-r0Vout/wr0*cos(angel)+r0/wr0*sin(angel)*theVout-r0Vout*sin(theta0)/wr0*wr0t2-r0*cos(theta0)*theVout/wr0*wr0t2-1/(2*wr0*Ln)*wr0t2^2)];
%     

% M=[Ir0N/r0*sin(angel)+r0*cos(angel)*(Vcr0N-(1-M))/r0^2 (Vcr0N-(1-M))/r0*sin(angel)+r0*cos(angel)*(-Ir0N)/r0^2 -1 0 r0*wr0*cos(angel) (Vcr0N-(1-M))*n/Vin/r0*sin(angel)+r0*cos(angel)*(-Ir0N)*n/Vin/r0^2;...
%     -Ir0N/r0*cos(angel)+r0*sin(angel)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)/(Vcr0N-(1-M)) -(Vcr0N-(1-M))/r0*cos(angel)+r0*sin(angel)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)/(Vcr0N-(1-M))^2 0 -1 -r0*wr0*sin(angel) -(Vcr0N-(1-M))*n/Vin/r0*cos(angel)+r0*wr0*sin(angel)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)*n/Vin/(Vcr0N-(1-M))^2-n/Vin;...
%     1 0 Ir2N/r1*sin(wr0*(ts/2-t2)+theta1)+r1/sqrt(1+Ln)*cos(wr0*(ts/2-t2)+theta1)*1/(1+(sqrt(1+Ln)*Ir0N/(Vcr2N-1))^2)*sqrt(1+Ln)/(Vcr2N-1) (Vcr2N-1)/(sqrt(1+Ln)*r1)*sin(wr0*(ts/2-t2)+theta1)+r1/sqrt(1+Ln)*cos(wr0*(ts/2-t2)+theta1)*1/(1+(sqrt(1+Ln)*Ir0N/(Vcr2N-1))^2)*(-sqrt(1+Ln)*Ir2N)/(Vcr2N-1)^2 -r1/(sqrt(1+Ln))*wr1*cos(wr0*(ts/2-t2)+theta1) 0;...
%     0 1 -(sqrt(1+Ln)*Ir2N)/r1*cos(wr0*(ts/2-t2)+theta1)+r1*sin(wr0*(ts/2-t2)+theta1)*1/(1+(sqrt(1+Ln)*Ir2N/(Vcr2N-1))^2)*sqrt(1+Ln)/(Vcr2N-1) -(Vcr2N-1)/r1*cos(wr0*(ts/2-t2)+theta1)+r1*sin(wr0*(ts/2-t2)+theta1)*1/(1+(sqrt(1+Ln)*Ir2N/(Vcr2N-1))^2)*(-sqrt(1+Ln)*Ir2N)/(Vcr2N-1)^2 r1*wr1*cos(wr1*t2+theta1) 0;...
%     Ir0N/r0*sin(angel)+r0*cos(angel)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)/(Vcr0N-(1-M))-Ir0N/r0*sin(theta0)-r0*cos(theta0)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)/(Vcr0N-(1-M)) (Vcr0N-(1-M))/r0*sin(angel)+r0*cos(angel)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)/(Vcr0N-(1-M))^2-(Vcr0N-(1-M))/r0*sin(theta0)-r0*cos(theta0)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)/(Vcr0N-(1-M))^2 0 0 r0*wr0*cos(angel)-M/Ln*wr0 (Vcr0N-(1-M))*n/Vin/r0*sin(angel)-(Vcr0N-(1-M))*n/Vin/r0*sin(theta0)-n/(Vin*Ln)*wr0*t2;...
%     -Nnorm*((Ir0N/r0*cos(theta0)-r0*sin(theta0)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)/(Vcr0N-(1-M))-Ir0N/r0*cos(angel)+r0*sin(angel)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)/(Vcr0N-(1-M)))/wr0-Ir0N/r0*sin(theta0)+r0*cos(theta0)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)/(Vcr0N-(1-M)))  -Nnorm*(((Vcr0N-(1-M))/r0*cos(theta0)-r0*sin(theta0)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)/(Vcr0N-(1-M))^2-(Vcr0N-(1-M))/r0*cos(angel)+r0*sin(angel)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)/(Vcr0N-(1-M))^2)/wr0-(Vcr0N-(1-M))/r0*sin(theta0)+r0*cos(theta0)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)/(Vcr0N-(1-M))^2) ...
%     0 0 -Nnorm*(r0*sin(angel)-r0*sin(theta0)-M/Ln*wr0*t2) 1-Nnorm*((Vcr0N-(1-M))*n/Vin/r0*cos(theta0)-r0*wr0*sin(theta0)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)*n/Vin/(Vcr0N-(1-M))^2-(Vcr0N-(1-M))*n/Vin/r0*cos(angel)+r0*wr0*sin(angel)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)*n/Vin/(Vcr0N-(1-M))^2)/wr0-Nnorm*((Vcr0N-(1-M))*n/Vin/r0*sin(theta0)+r0*cos(theta0)*1/(1+(Ir0N/(Vcr0N-(1-M)))^2)*(-Ir0N)*n/Vin/(Vcr0N-(1-M))^2)-Nnorm*n/Vin/(2*Ln)*wr0*t2^2];

end












