function eq = pm_equilibria(p)
% Feasible equilibria, including unstable equilibria, of the original ODE.
% z is normalized by K; eigenvalues refer to the full three-population system.
% These formulas cover the nonsingular, positive-mortality figure scenarios.
% Continuum cases and mu=0,1 require the separate analyses in the manuscript.
A=p.A; rR=p.rmax*exp(-p.g*p.uR); D=p.dP+p.m/p.kP;
assert(D>0 && p.mu>0 && p.mu<1,'Require D>0 and 0<mu<1.');
U=(p.rmax-p.d-p.m/p.k)/p.rmax;
W=(rR-p.d-p.m/(p.k+p.b*p.uR))/rR;
w=D/(1-p.mu); V=1-w/p.rP;
M22=det(A([1,3],[1,3])); M21=det(A([1,3],[2,3]));
M23=det(A([1,3],[1,2])); M33=det(A(1:2,1:2));
assert(abs(M22)>1e-12 && abs(M33)>1e-12,...
    'Singular minor: these isolated-equilibrium formulas are not applicable.');
eq=struct('name',{},'z',{},'eigenvalues',{},'stable',{},'residual',{},'beta',{});
eq=append(eq,'TE',zeros(3,1),p);
if U/A(1,1)>1e-10, eq=append(eq,'FSE',[U/A(1,1);0;0],p); end
if W/A(2,2)>1e-10, eq=append(eq,'FRE',[0;W/A(2,2);0],p); end
sr=A(1:2,1:2)\[U;W];
if all(sr>1e-10), eq=append(eq,'CSR',[sr;0],p); end
for kind={'CRP','CAT'}
    [a,B,C]=pm_coefficients(p,kind{1});
    c=[a,polyval(fliplr(B),w),polyval(fliplr(C),w)];
    first=find(abs(c)>1e-13,1);
    assert(~isempty(first),'Equilibrium polynomial vanishes identically.');
    rr=roots(c(first:end)); rr=sort(real(rr(abs(imag(rr))<1e-8)));
    for j=1:numel(rr)
        r=rr(j);
        if strcmp(kind{1},'CRP')
            z=[0;r;(V-A(3,2)*r)/A(3,3)];
            feasible=all(z(2:3)>1e-10);
        else
            z=[(A(3,3)*U-A(1,3)*V-M21*r)/M22;r;...
                (A(1,1)*V-A(3,1)*U-M23*r)/M22];
            feasible=all(z>1e-10);
        end
        if feasible, eq=append(eq,sprintf('%s%d',kind{1},j),z,p); end
    end
end
end

function eq=append(eq,name,z,p)
[f,J]=pm_model(p.K*z,p); ev=eig(J);
eq(end+1)=struct('name',name,'z',z,'eigenvalues',ev,...
    'stable',max(real(ev))<-1e-9,'residual',norm(f/p.K,inf),...
    'beta',real(poly(J)));
end
