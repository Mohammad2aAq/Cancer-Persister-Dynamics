function [f,J,birthFactor] = pm_model(x,p)
% Original manuscript equations in cell counts x=[x_S;x_R;x_P].
% No clipping of populations or mutation flux is performed.
x=x(:); z=x/p.K;
rR=p.rmax*exp(-p.g*p.uR);
cS=p.rmax-p.d-p.m/p.k;
cR=rR-p.d-p.m/(p.k+p.b*p.uR);
D=p.dP+p.m/p.kP;
birthFactor=1-p.A(3,:)*z;
gS=cS-p.rmax*p.A(1,:)*z;
gR=cR-rR*p.A(2,:)*z;
t=(1-p.mu)*p.rP; v=p.mu*p.rP;
f=[x(1)*gS; x(2)*gR+v*x(3)*birthFactor;...
    x(3)*(t*birthFactor-D)];
if nargout>1
    J=[-p.rmax*z(1)*p.A(1,:);...
        -rR*z(2)*p.A(2,:)-v*z(3)*p.A(3,:);...
        -t*z(3)*p.A(3,:)];
    J(1,1)=J(1,1)+gS;
    J(2,2)=J(2,2)+gR;
    J(2,3)=J(2,3)+v*birthFactor;
    J(3,3)=J(3,3)+t*birthFactor-D;
end
end
