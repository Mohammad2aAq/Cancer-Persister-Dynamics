function p = pm_parameters(m,mu,scenario)
% Parameters used in the manuscript; populations are cell counts.
if nargin<1, m=0.45; end
if nargin<2, mu=0.5; end
if nargin<3, scenario='main'; end
p=struct('rmax',.45,'rP',.2,'K',1e4,'d',.01,'dP',.001,...
    'b',10,'k',2,'kP',10,'g',.8,'uR',3.3,'m',m,'mu',mu,...
    'A',[1,1.6,.2;1.1,1,1.7;.3,1.9,1]);
switch scenario
    case 'main'
    case 'weak'
        p.A=[1,.6,.2;.1,1,.7;.3,.9,1];
    otherwise
        error('Unknown parameter scenario: %s',scenario);
end
end
