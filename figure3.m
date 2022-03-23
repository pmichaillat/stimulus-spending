%
% Produce figure 3
%

close all;clear;

%% Calibration of formulas (23) and (24) <> see section 4

eta=0.6;
gy=@(gc)gc./(1+gc);
gystar=0.165;
gcstar=gystar./(1-gystar);
ustar=0.06;
u0=0.09;
tau0=0.017;
coeff=gystar.*eta./(1-eta).*tau0./u0;
M=[0:0.01:3]; % range of unemployment multipliers
m=M.*(1-u0)./(1-coeff.*M);
z1=(1-gystar).*gystar./ustar;
z0=1./((1-eta).*(1-ustar)^2);
du=(u0-ustar)./ustar;

%% Compute stimulus spending and unemployment rate from formulas (23) and (24) 

% Stimulus spending and unemployment rate with epsilon = 1, for a range of unemployment multipliers

epsilon=1;
dgc=epsilon.*z0.*m./(1+z1.*z0.*epsilon.*m.^2).*du; %formula (23)
gc=gcstar.*(1+dgc);
gy1=gy(gc);
u1=ustar+du.*ustar./(1+z1.*z0.*epsilon.*m.^2); %formula (24)

% Stimulus spending and unemployment rate with epsilon = 2, for a range of unemployment multipliers

epsilon=2;
dgc=epsilon.*z0.*m./(1+z1.*z0.*epsilon.*m.^2).*du; %formula (23)
gc=gcstar.*(1+dgc);
gy2=gy(gc);
u2=ustar+du.*ustar./(1+z1.*z0.*epsilon.*m.^2); %formula (24)

% Stimulus spending and unemployment rate with epsilon = 0.5, for a range of unemployment multipliers

epsilon=0.5;
dgc=epsilon.*z0.*m./(1+z1.*z0.*epsilon.*m.^2).*du; %formula (23)
gc=gcstar.*(1+dgc);
gy05=gy(gc);
u05=ustar+du.*ustar./(1+z1.*z0.*epsilon.*m.^2); %formula (24)

%% Plot figure 3

% Formatting

set(groot,'DefaultFigureUnits', 'inches')
set(groot,'DefaultFigurePosition', [0,0,7.7778,5.8333]);
set(groot,'DefaultFigurePaperPosition', [0, 0, 7.7779,5.8334]);
set(groot,'DefaultFigurePaperSize', [7.7779,5.8334]);
set(groot,'DefaultAxesFontName', 'Times New Roman')
set(groot,'DefaultAxesFontSize', 20)
set(groot,'DefaultLineLineWidth', 3)
set(groot,'DefaultAxesLineWidth', 1)
set(groot,'DefaultAxesXColor','k')
set(groot,'DefaultAxesYColor','k')
set(groot,'DefaultAxesYGrid','on')
set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesTickLength',[0 0])
co = [33,113,181
			8,48,107
			107,174,214];
co = co./255;
set(groot,'DefaultAxesColorOrder',co)

% Figure 3A : optimal stimulus spending

y1=(gy1-gystar).*100;
y2=(gy2-gystar).*100;
y3=(gy05-gystar).*100;

figure(1)
clf
hold on
plot(M,y1,'-','Color',co(1,:))
plot(M,y2,'-','Color',co(2,:))
plot(M,y3,'-','Color',co(3,:))
set(gca,'YLim',[0,6],'YTick',[0:2:6],'YTickLabel',['0%';'2%';'4%';'6%'])
set(gca,'XLim',[0,2],'XTick',[0:0.5:2])
xlabel('Unemployment multiplier')
ylabel('Optimal stimulus spending (% of GDP)')
print('-dpdf',['figure3A.pdf'])

% Figure 3B: unemployment rate

y1=u1.*100;
y2=u2.*100;
y3=u05.*100;

figure(2)
clf
hold on
plot(M,y1,'-','Color',co(1,:))
plot(M,y2,'-','Color',co(2,:))
plot(M,y3,'-','Color',co(3,:))
set(gca,'YLim',[6,9],'YTick',[6:1:9],'YTickLabel',['6%';'7%';'8%';'9%'])
set(gca,'XLim',[0,2],'XTick',[0:0.5:2])
xlabel('Unemployment multiplier')
ylabel('Resulting unemployment rate')
print('-dpdf',['figure3B.pdf'])