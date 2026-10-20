function d = recompute_manuscript_data()
% Recompute Figures 1--7 of our manuscript.
d.equilibrium_names={'TE','FSE','FRE','CSR','CRP1','CRP2','CAT1','CAT2'};
raw=cell(1,2); horizontal={linspace(.005,.995,121),linspace(0,10,121)};
mvals=linspace(.005,.995,121);
for panel=1:2
    grid=cell(121,121);
    for i=1:121
        for j=1:121
            m=mvals(i); x=horizontal{panel}(j);
            if panel==1, p=pm_parameters(m,x);
            else, p=pm_parameters(m,(1-m)/(1+m)); p.uR=x;
            end
            eq=pm_equilibria(p); tags={};
            for e=eq
                if e.stable
                    tag=regexprep(e.name,'[12]$',''); tags{end+1}=tag; %#ok<AGROW>
                end
            end
            grid{i,j}=strjoin(sort(unique(tags)),' + ');
        end
    end
    raw{panel}=grid; fprintf('Stability map %d/2 calculated.\n',panel);
end
allsets=unique([raw{1}(:);raw{2}(:)]);
allsets=allsets(~cellfun(@isempty,allsets)); classes=[{''};sort(allsets)];
labels=classes; labels{1}='No stable equilibrium';
for panel=1:2
    [~,classification]=ismember(raw{panel},classes);
    d.maps(panel)=struct('horizontal',horizontal{panel},'m',mvals,...
        'classes',{labels},'classification',classification);
end
d.branches=pm_branch_rows(linspace(.01,.99,1961),.45);
critical=pm_bifurcations(); d.events=critical.events; d.hopf=critical.hopf;
d.folds=critical.folds;
for j=1:2
    f=d.folds(j); rows=pm_branch_rows(linspace(f.mu-.018,f.mu+.018,500),f.m);
    if j==1, keep=ismember(rows(:,2),[5,6]); else, keep=ismember(rows(:,2),[7,8]); end
    d.folds(j).branches=rows(keep,:);
end
mus=[.45,.5,.65]; stops=[12000,50000,12000];
for j=1:3
    d.trajectories(j)=pm_trajectory(pm_parameters(.45,mus(j)),[1000,1000,1000],stops(j),6001);
end
p=pm_parameters(.3,.9,'weak'); initials=[6500,20,20;20,4000,20;1000,1000,1000];
for j=1:3, d.csr.trajectories(j)=pm_trajectory(p,initials(j,:),6000); end
eq=pm_equilibria(p); eq=eq(strcmp({eq.name},'CSR'));
assert(numel(eq)==1); d.csr.equilibrium=eq.z;
d.csr.threshold=1-(p.dP+p.m/p.kP)/(p.rP*(1-p.A(3,:)*eq.z));
d.csr.branches=pm_branch_rows(linspace(.55,.9,351),.3,'weak');
fprintf('Equilibria, bifurcations, and trajectories calculated.\n');
end
