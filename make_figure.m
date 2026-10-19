function fig = make_figure(number,d,outdir)
% MAKE_FIGURE(1), ..., MAKE_FIGURE(7): final-manuscript figure numbers.
% Figure 8 is supplied separately by the author.
% Optional d supplies newly calculated data; outdir selects the export folder.
% Parameter values and quantities follow the manuscript; styling may differ.
% Maps classify strictly stable equilibrium sets on the archived interior grid.
root=fileparts(mfilename('fullpath'));
if nargin<2 || isempty(d)
    loaded=load(fullfile(root,'data','paper_data.mat'),'d'); d=loaded.d;
end
if nargin<3, outdir=fullfile(root,'figures_paper'); end
assert(isscalar(number) && ismember(number,1:7),...
    'Use figure numbers 1 to 7; Figure 8 has separate author-supplied code.');
names={'figure_01_stability_mu_m','figure_02_stability_uR_m',...
    'figure_03_bifurcation_branches','figure_04_hopf_crossings',...
    'figure_05_saddle_nodes','figure_06_trajectories','figure_07_corrected_csr'};
heights=[5.8,5.8,6.5,3.7,3.8,7.2,3.7];
col=[33,102,172;178,24,43;118,42,131]/255;
fig=figure('Color','w','Name',names{number},'NumberTitle','off',...
    'Units','inches','Position',[.5,.5,8.5,heights(number)]);
switch number
    case {1,2}
        cmap=[255,255,255;33,102,172;178,24,43;118,42,131;...
            27,120,55;224,130,20;67,147,195;153,112,171;...
            128,205,193;223,194,125;197,27,125]/255;
        s=d.maps(number); classes=s.classes; n=numel(classes);
        assert(n<=size(cmap,1),'More stable sets than available map colors.');
        ax=axes(fig,'Position',[.13,.30,.80,.62]);
        imagesc(ax,s.horizontal,s.m,s.classification,[.5,n+.5]);
        set(ax,'YDir','normal'); colormap(ax,cmap(1:n,:));
        if number==1
            lab(ax,'$\mu$','$m$','$u_R=3.3$');
        else
            lab(ax,'$u_R$','$m$','$\mu=(1-m)/(1+m)$');
        end
        hold(ax,'on'); handles=gobjects(1,n);
        for j=1:n
            handles(j)=plot(ax,nan,nan,'s','MarkerFaceColor',cmap(j,:),...
                'MarkerEdgeColor',[.5,.5,.5],'MarkerSize',9);
        end
        classes{1}='No stable equilibrium';
        lg=legend(ax,handles,classes,'NumColumns',3,'Interpreter','none','FontSize',10);
        lg.Units='normalized'; lg.Position=[.08,.02,.84,.19];
    case 3
        r=d.branches;
        r=r(~ismember(r(:,2),[1,3]) & ~(r(:,2)==2 & r(:,1)<.68),:);
        labels={'$x_S/K$','$x_R/K$','$x_P/K$','Largest real part of eigenvalues'};
        for j=1:4
            ax=subplot(2,2,j,'Parent',fig); hold(ax,'on');
            branchplot(ax,r,j+2,col,2);
            for mu=d.events(:)', xline(ax,mu,':','Color',[.6,.6,.6],'LineWidth',.6); end
            lab(ax,'$\mu$',labels{j},'');
            if j==4, ylim(ax,[-.015,.015]); yline(ax,0,'k-','LineWidth',.5); end
        end
        sgtitle(fig,'Stable equilibria: blue; unstable equilibria: red','FontSize',12);
    case 4
        for j=1:2
            ax=subplot(1,2,j,'Parent',fig); hold(ax,'on'); s=d.hopf(j);
            plot(ax,s.horizontal,s.realpart,'Color',col(1,:),'LineWidth',1.2);
            yline(ax,0,'k-','LineWidth',.6); xline(ax,s.mu,'--','Color',col(2,:));
            lab(ax,'$\mu$','Real part of the critical pair',sprintf('$\\mu=%.9f$',s.mu));
        end
    case 5
        for j=1:2
            ax=subplot(1,2,j,'Parent',fig); hold(ax,'on'); s=d.folds(j);
            branchplot(ax,s.branches,5,col,3);
            plot(ax,s.mu,s.z(3),'ko','MarkerFaceColor','k','MarkerSize',5);
            xline(ax,s.mu,'--','Color',[.5,.5,.5],'LineWidth',.8);
            lab(ax,'$\mu$','$x_P/K$',sprintf('%s: $m=%.2f$, $\\mu=%.9f$',s.kind,s.m,s.mu));
        end
        sgtitle(fig,'Saddle-node bifurcations; blue stable, red unstable','FontSize',12);
    case 6
        for j=1:3
            ax=subplot(3,1,j,'Parent',fig); s=d.trajectories(j);
            populationplot(ax,s,col,j==1);
            lab(ax,'Time','Population / $K$',sprintf('$m=0.45$, $\\mu=%.2f$',s.mu));
        end
    case 7
        ax=subplot(1,2,1,'Parent',fig); hold(ax,'on'); h=gobjects(1,4); labels=cell(1,4);
        for j=1:3
            s=d.csr.trajectories(j);
            h(j)=plot(ax,s.populations(:,1)/s.K,s.populations(:,2)/s.K,'Color',col(j,:));
            plot(ax,s.initial(1)/s.K,s.initial(2)/s.K,'o','Color',col(j,:),'MarkerSize',3);
            labels{j}=sprintf('(%g, %g, %g)',s.initial);
        end
        z=d.csr.equilibrium; h(4)=plot(ax,z(1),z(2),'k*','MarkerSize',10); labels{4}='Stable CSR';
        legend(ax,h,labels,'Location','best','FontSize',9);
        lab(ax,'$x_S/K$','$x_R/K$','$m=0.30$, $\mu=0.90$');
        ax=subplot(1,2,2,'Parent',fig); hold(ax,'on'); r=d.csr.branches;
        r=r(ismember(r(:,2),[4,7,8]),:); branchplot(ax,r,5,col,3);
        h=xline(ax,d.csr.threshold,'--','Color',col(2,:));
        legend(ax,h,sprintf('$\\mu=%.9f$',d.csr.threshold),'Interpreter','latex','FontSize',10);
        lab(ax,'$\mu$','$x_P/K$','Persister invasion at CSR');
end
drawnow;
if ~exist(outdir,'dir'), mkdir(outdir); end
base=fullfile(outdir,names{number});
set(fig,'PaperUnits','inches','PaperSize',[8.5,heights(number)],...
    'PaperPosition',[0,0,8.5,heights(number)]);
print(fig,[base '.pdf'],'-dpdf','-painters');
print(fig,[base '.png'],'-dpng','-r180');
savefig(fig,[base '.fig']);
fprintf('Figure %d: %s\n',number,base);
end

function branchplot(ax,r,column,col,markerSize)
for stable=[false,true]
    take=logical(r(:,7))==stable; color=col(2-double(stable),:);
    plot(ax,r(take,1),r(take,column),'.','Color',color,'MarkerSize',markerSize);
end
end

function populationplot(ax,s,col,showLegend)
hold(ax,'on'); h=gobjects(1,3);
for j=1:3, h(j)=plot(ax,s.time,s.populations(:,j)/s.K,'Color',col(j,:),'LineWidth',1); end
if showLegend
    legend(ax,h,{'$x_S/K$','$x_R/K$','$x_P/K$'},'Interpreter','latex','Location','best','FontSize',10);
end
end

function lab(ax,x,y,t)
set(ax,'FontSize',11,'Box','off','TickLabelInterpreter','latex');
xlabel(ax,x,'Interpreter','latex','FontSize',13);
ylabel(ax,y,'Interpreter','latex','FontSize',12);
title(ax,t,'Interpreter','latex','FontSize',12);
end

