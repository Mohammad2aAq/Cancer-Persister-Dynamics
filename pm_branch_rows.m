function rows = pm_branch_rows(mus,m,scenario)
% Columns: mu, equilibrium code, xS/K, xR/K, xP/K, max Re(lambda), stable, residual.
if nargin<3, scenario='main'; end
names={'TE','FSE','FRE','CSR','CRP1','CRP2','CAT1','CAT2'};
rows=zeros(numel(mus)*8,8); n=0;
for mu=mus(:)'
    eq=pm_equilibria(pm_parameters(m,mu,scenario));
    for e=eq
        n=n+1; code=find(strcmp(names,e.name));
        rows(n,:)=[mu,code,e.z',max(real(e.eigenvalues)),double(e.stable),e.residual];
    end
end
rows=rows(1:n,:);
end
