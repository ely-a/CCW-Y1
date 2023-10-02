function K=getK(Vel,Kh,Kalfa,c,xf,rho)
%This function returns the stiffness matrix (K) of the aeroelastic 
%equations of motion in 2 degrees of freedom (M xdot = K x). 
%
%[K]=getK(Vel,Kh,Kalfa,c,xf,rho)
%
%outputs: K     -   stiffness matrix of dimensions 8x8
%inputs:  Vel   -   Airflow velocity (m/s)
%         Kh    -   sectional stiffness in heave (N/m)
%         Kalfa -   sectional torsional stiffness (Nm/rad)
%         c     -   airfoil chord length (m)
%         xf    -   sectional flexural centre position from l.e. (% chord)
%         rho   -   density of air (kg/m^3)


%define wagner function approximation constants
psi1=0.165;
psi2=0.335;
eps1=0.0455;
eps2=0.3;
phi0=1-psi1-psi2;
phidot0=1+psi1*eps1*Vel*2/c+psi2*eps2*Vel*2/c;

b=c/2;
ec=xf*c-c/4;
%define damping submatrix
C=[phi0,c/4+phi0*(3*c/4-xf);-ec*phi0,(3*c/4-xf*c)*(c/4-ec*phi0)];

C=C*rho*pi*Vel*c;
%define stiffness submatrix
K(1,1)=Kh-pi*rho*Vel*c*phidot0;
K(1,2)=pi*rho*Vel*c*(Vel*phi0-(3/4-xf)*c*phidot0);
K(2,1)=pi*rho*Vel*ec*c*phidot0;
K(2,2)=Kalfa-pi*rho*Vel*ec*c*(Vel*phi0-(3/4-xf)*c*phidot0);

%define aerodynamic submatix
W(1,1:2)=[-psi1*eps1^2/b,-psi2*eps2^2/b];
W(1,3)=psi1*eps1*(1-eps1*(1-2*ec/c));
W(1,4)=psi2*eps2*(1-eps2*(1-2*ec/c));
W(2,1:2)=-ec*W(1,1:2);

W=W*2*pi*rho*Vel^3;
%define Wo submatrix
Wo=zeros(4,6);
Wo([1,2],1)=1;
Wo([3,4],2)=1;
Wo(1,3)=-eps1*Vel*2/c;
Wo(2,4)=-eps2*Vel*2/c;
Wo(3,5)=-eps1*Vel*2/c;
Wo(4,6)=-eps2*Vel*2/c;

K=[C,K,W;eye(2),zeros(2,6);zeros(4,2),Wo];
return