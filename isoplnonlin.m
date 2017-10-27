function [q_e, k_e]=isoplnonlin(E, nu, th, x, d)
% ISOPLNONLIN calculates and returns the internal force vector and tangent
% stiffness matrix for 0an non-linear 2D planar element.

C=E/(1-nu^2)*[1 nu 0                          
              nu 1 0
              0 0 (1-nu)/2];

A1=[1 0 0 0
    0 0 0 1
    0 1 1 0];

O2=zeros(2,2);

k_e=zeros(8,8);

a=1/sqrt(3);

xi=[-a -a  a a];
mu=[-a  a -a a];

for i=1:4

Dn=1/4*[-(1-mu(i))  (1-mu(i)) (1+mu(i)) -(1+mu(i))
        -(1-xi(i)) -(1+xi(i)) (1+xi(i))  (1-xi(i))];

J=Dn*[x(1) x(2)
      x(3) x(4)
      x(5) x(6)
      x(7) x(8)];

Ga=inv(J);      
      
A2=[Ga O2
    O2 Ga];

A3=1/4*[-(1-mu(i))   0     (1-mu(i))  0     (1+mu(i))  0    -(1+mu(i))  0
        -(1-xi(i))   0    -(1+xi(i))  0     (1+xi(i))  0     (1-xi(i))  0
          0     -(1-mu(i))  0     (1-mu(i))  0     (1+mu(i))  0    -(1+mu(i))
          0     -(1-xi(i))  0    -(1+xi(i))  0     (1+xi(i))  0     (1-xi(i))];
      
B=A1*A2*A3;


k_e=k_e+B'*C*B*det(J);

end;

k_e=k_e*th;