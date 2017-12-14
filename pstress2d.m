% 2D plasticity

% calculate the new stresses (sn) and new equivalent 
% plastic strain (epsn) from the old stresses (so), 
% the old equivalent plastic strain (epso)
% and the strain increment (de);

function [sn,epsn,Ct]=pstress2d(so,epso,de)


global E nu H yield

Ce=E/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

sb=so+Ce*de;
sy=yield+H*epso;
fb=sqrt(sb(1)^2+sb(2)^2-sb(1)*sb(2)+3*sb(3)^2)-sy;

if fb<=0
   sn=sb;
   epsn=epso;
   Ct=Ce;   
   iter=0;
end

if fb>0
   ab=1/2/sqrt(sb(1)^2+sb(2)^2-sb(1)*sb(2)+3*sb(3)^2)*[2*sb(1)-sb(2)
                                                       2*sb(2)-sb(1)
                                                       6*sb(3)];
   dl=fb/(ab'*Ce*ab+H);
   s=sb-dl*Ce*ab;
   eps=epso+dl;
   sy=yield+H*eps;
   f=sqrt(s(1)^2+s(2)^2-s(1)*s(2)+3*s(3)^2)-sy;
   iter=0;   
   test=0;
   while f>1E-5 | test==0
      test=1;
      iter=iter+1;
      a=1/2/sqrt(s(1)^2+s(2)^2-s(1)*s(2)+3*s(3)^2)*[2*s(1)-s(2)
                                                    2*s(2)-s(1)
                                                    6*s(3)];
      r=s-(sb-dl*Ce*a);
      
      ad=3/4/(s(1)^2+s(2)^2-s(1)*s(2)+3*s(3)^2)^(3/2)* [s(2)^2+4*s(3)^2       -s(1)*s(2)-2*s(3)^2   -2*s(3)*(2*s(1)-s(2))  
                                                       -s(1)*s(2)-2*s(3)^2     s(1)^2+4*s(3)^2       2*s(3)*(-2*s(2)+s(1))
                                                        -2*s(3)*(2*s(1)-s(2))  2*s(3)*(-2*s(2)+s(1)) 4*(s(1)^2+s(2)^2-s(1)*s(2))];
      Q=(eye(3)+dl*Ce*ad);
      lp=(f-a'*inv(Q)*r)/(a'*inv(Q)*Ce*a+H);
      s=s-inv(Q)*r-lp*inv(Q)*Ce*a;
      dl=dl+lp;
      eps=epso+dl;
      sy=yield+H*eps;
      f=sqrt(s(1)^2+s(2)^2-s(1)*s(2)+3*s(3)^2)-sy;
   end 
   iter;
   sn=s;
   epsn=eps;
                                          
   R=inv(Q)*Ce;
   Ct=R*(eye(3)-a*a'*R/(a'*R*a+H));
   
end   


