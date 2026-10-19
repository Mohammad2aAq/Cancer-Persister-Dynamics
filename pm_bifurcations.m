function result = pm_bifurcations()
% Locate the bifurcation values plotted in Figures 2, 3, 4 and 6.
% These brackets and parameter scenarios are those of the paper; this is not
% a global continuation package for arbitrary edited parameters.
brackets=[.474,.480;.55,.565]; spans=[.006,.025];
hopf=struct([]);
for j=1:2
    h=fzero(@hopf_condition,brackets(j,:),optimset('TolX',1e-13));
    p=pm_parameters(.45,h); e=cat_equilibrium(p);
    ev=e.eigenvalues; [omega,k]=max(imag(ev));
    assert(omega>0 && abs(real(ev(k)))<1e-9 && e.beta(2)>0);
    horizontal=linspace(h-spans(j),h+spans(j),401); rr=zeros(size(horizontal));
    for i=1:numel(horizontal)
        ee=cat_equilibrium(pm_parameters(.45,horizontal(i)));
        [im,k]=max(imag(ee.eigenvalues)); assert(im>0);
        rr(i)=real(ee.eigenvalues(k));
    end
    step=1e-6;
    slope=-(hopf_condition(h+step)-hopf_condition(h-step))/(2*step)...
        /(2*(e.beta(2)^2+e.beta(3)));
    hopf(j)=struct('mu',h,'horizontal',horizontal,'realpart',rr,...
        'z',e.z,'omega',omega,'slope',slope);
end

folds=struct([]); kinds={'CRP','CAT'}; doses=[.2,.35];
for j=1:2
    p=pm_parameters(doses(j)); A=p.A; D=p.dP+p.m/p.kP;
    U=(p.rmax-p.d-p.m/p.k)/p.rmax;
    M22=det(A([1,3],[1,3])); M21=det(A([1,3],[2,3]));
    M23=det(A([1,3],[1,2]));
    [a,B,C]=pm_coefficients(p,kinds{j});
    discriminant=conv(B,B)-4*a*C;
    ww=roots(fliplr(discriminant)); count=0;
    for w=ww(:)'
        if abs(imag(w))>1e-9 || real(w)<=0, continue; end
        w=real(w); mu=1-D/w; V=1-w/p.rP;
        r=-polyval(fliplr(B),w)/(2*a);
        if j==1
            z=[0;r;(V-A(3,2)*r)/A(3,3)]; feasible=all(z(2:3)>0);
        else
            z=[(A(3,3)*U-A(1,3)*V-M21*r)/M22;r;...
                (A(1,1)*V-A(3,1)*U-M23*r)/M22]; feasible=all(z>0);
        end
        if mu<=0 || mu>=1 || ~feasible, continue; end
        count=count+1;
        folds(j)=struct('mu',mu,'m',doses(j),'kind',kinds{j},'z',z,'branches',[]);
    end
    assert(count==1,'Expected one physical fold in this paper scenario.');
end

% Sensitive-invasion boundaries on CRP, as a quadratic in w.
p=pm_parameters(); A=p.A; rR=p.rmax*exp(-p.g*p.uR);
D=p.dP+p.m/p.kP; U=(p.rmax-p.d-p.m/p.k)/p.rmax;
cR=rR-p.d-p.m/(p.k+p.b*p.uR); T=det(A([1,3],[2,3]));
Rn=[A(3,3)*U-A(1,3),A(1,3)/p.rP];
Pn=[A(1,2)-A(3,2)*U,-A(1,2)/p.rP];
H=conv(Rn,[cR*T,0]-rR*(A(2,2)*Rn+A(2,3)*Pn))+T*conv([-D,1],Pn);
ww=roots(fliplr(H)); collisions=[];
for w=ww(:)'
    if abs(imag(w))>1e-9 || real(w)<=0, continue; end
    w=real(w); mu=1-D/w;
    if mu>0 && mu<1 && polyval(fliplr(Rn),w)/T>0 && polyval(fliplr(Pn),w)/T>0
        collisions(end+1)=mu; %#ok<AGROW>
    end
end
assert(numel(collisions)==2);
muFRE=1-D/(p.rP*(1-A(3,2)*cR/(rR*A(2,2))));
muFSE=1-D/(p.rP*(1-A(3,1)*U/A(1,1)));
result=struct('hopf',hopf,'folds',folds,...
    'events',sort([collisions,[hopf.mu],muFRE,muFSE]));
end

function value=hopf_condition(mu)
e=cat_equilibrium(pm_parameters(.45,mu));
value=e.beta(2)*e.beta(3)-e.beta(4);
end

function e=cat_equilibrium(p)
ee=pm_equilibria(p); ee=ee(startsWith({ee.name},'CAT'));
assert(numel(ee)==1,'Expected a unique feasible CAT branch near this Hopf point.');
e=ee(1);
end
