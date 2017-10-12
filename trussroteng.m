% non-linear bar
% green's strain

function [q,K]=trussroteng(x,p,E,Ao);

x21=x(3)-x(1);
y21=x(4)-x(2);
u21=p(3)-p(1);
v21=p(4)-p(2);

Lo=sqrt(x21^2+y21^2);
Ln=sqrt((x21+u21)^2+(y21+v21)^2);

eps=(Ln-Lo)/Lo;

sig=E*eps;

b=1/(Lo*Ln)*[-(x21+u21) -(y21+v21) (x21+u21) (y21+v21)];
c=[-(x21+u21) -(y21+v21) (x21+u21) (y21+v21)];
Ka=1/(Lo*Ln)*[1  0 -1  0
           0  1  0 -1
          -1  0  1  0
           0 -1  0  1];
A=[1  0 -1  0
           0  1  0 -1
          -1  0  1  0
           0 -1  0  1];
q=Ao*Lo*sig*b';

K=(Ao*E/(Lo*Ln^2)-sig*E/Ln^3)*c'*c+Ao*sig*A/Ln;

