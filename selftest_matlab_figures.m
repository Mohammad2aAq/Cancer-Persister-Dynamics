function selftest_matlab_figures()
% Numerical checks to run in MATLAB before performing modified experiments.
% Checks the retained model, equilibrium and bifurcation code; does not plot figures.
root=fileparts(mfilename('fullpath'));
paper=load(fullfile(root,'data','paper_data.mat'),'d'); d=paper.d;
required={'equilibrium_names','maps','branches','events','trajectories','hopf','folds','csr'};
assert(isequal(sort(fieldnames(d)),sort(required(:))),'Unexpected reference-data fields.');
p=pm_parameters(); x=[2000;1000;1500]; [f,J]=pm_model(x,p);
Jfd=zeros(3); step=1e-3;
for j=1:3
    delta=zeros(3,1); delta(j)=step;
    Jfd(:,j)=(pm_model(x+delta,p)-pm_model(x-delta,p))/(2*step);
end
assert(norm(J-Jfd,inf)<1e-8,'Jacobian finite-difference check failed.');
growth=1-p.A*x/p.K; rR=p.rmax*exp(-p.g*p.uR);
total=p.rmax*x(1)*growth(1)+rR*x(2)*growth(2)+p.rP*x(3)*growth(3)...
    -(p.d+p.m/p.k)*x(1)-(p.d+p.m/(p.k+p.b*p.uR))*x(2)...
    -(p.dP+p.m/p.kP)*x(3);
assert(abs(sum(f)-total)<1e-10,'Mutation flux did not cancel from total population.');
% Deliberately verify the original model's positivity limitation; never repair it.
q=pm_parameters(); f=pm_model(q.K*[4;0;.1],q);
assert(f(2)<0,'The original mutation flux appears to have been altered.');

mus=[.10,.30,.45,.48,.50,.65,.80];
actual=pm_branch_rows(mus,.45);
names=d.equilibrium_names;
for i=1:size(actual,1)
    rr=d.branches(abs(d.branches(:,1)-actual(i,1))<1e-10 & d.branches(:,2)==actual(i,2),:);
    assert(size(rr,1)==1,'Missing reference branch.');
    assert(max(abs(rr(3:6)-actual(i,3:6)))<1e-9,...
        ['Equilibrium comparison failed: ' names{actual(i,2)}]);
    assert(rr(7)==actual(i,7) && actual(i,8)<1e-10);
end
critical=pm_bifurcations();
assert(max(abs(critical.events(:)-d.events(:)))<1e-9,'Bifurcation locations disagree.');
for j=1:2
    assert(abs(critical.folds(j).mu-d.folds(j).mu)<1e-9);
    assert(norm(critical.folds(j).z(:)-d.folds(j).z(:),inf)<1e-8);
    q=pm_parameters(critical.folds(j).m,critical.folds(j).mu);
    [f,J]=pm_model(q.K*critical.folds(j).z,q);
    assert(norm(f/q.K,inf)<1e-10 && min(abs(eig(J)))<1e-8);
    assert(abs(critical.hopf(j).slope-d.hopf(j).slope)<1e-6);
end

% Compare one adaptive MATLAB trajectory with its independently saved LSODA result.
trajectory=pm_trajectory(pm_parameters(.45,.45),[1000,1000,1000],12000,6001);
err=abs(trajectory.populations-d.trajectories(1).populations)/p.K;
assert(max(err(:))<1e-5,'Adaptive trajectory comparison failed.');

fprintf('PASS: model, Jacobian, equilibria, bifurcations, and adaptive trajectory checks.\n');
end
