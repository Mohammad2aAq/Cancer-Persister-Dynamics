function [a,B,C] = pm_coefficients(p,kind)
% Equilibrium polynomial a*r^2+B(w)*r+C(w)=0, r=x_R/K.
% w=(d_P+m/k_P)/(1-mu). B and C use ASCENDING powers of w.
A=p.A; rR=p.rmax*exp(-p.g*p.uR);
D=p.dP+p.m/p.kP; U=(p.rmax-p.d-p.m/p.k)/p.rmax;
cR=rR-p.d-p.m/(p.k+p.b*p.uR);
M22=det(A([1,3],[1,3])); M23=det(A([1,3],[1,2]));
switch kind
    case 'CRP'
        a=-rR*det(A(2:3,2:3));
        B=[A(3,3)*cR-rR*A(2,3)+A(3,2)*D,...
            rR*A(2,3)/p.rP-A(3,2)];
        C=conv([-D,1],[1,-1/p.rP]);
    case 'CAT'
        sn=[A(3,3)*U-A(1,3),A(1,3)/p.rP];
        pn=[A(1,1)-A(3,1)*U,-A(1,1)/p.rP];
        a=-rR*det(A);
        B=-rR*(A(2,1)*sn+A(2,3)*pn)-M23*[-D,1];
        B(1)=B(1)+M22*cR;
        C=conv([-D,1],pn);
    otherwise
        error('Equilibrium kind must be CRP or CAT.');
end
end
