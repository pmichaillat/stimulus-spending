%
% Produce figure 4
%

close all;clear;

%% Calibration of the matching model with land to US data <> see sections 2.2 and 2.4, and online appendix A

epsilon=1;
Mbar=0.5;
eta=0.6;
s=0.028;
ubar=0.061;
xbar=0.43;
lbar=(1-ubar);
Ybar=lbar;
GCbar=0.197;
omega=s.*(lbar./ubar).*xbar.^(eta-1);
gamma=1./(1+1./GCbar);

GY=@(gc)gc./(1+gc);
CY=@(gc)1-GY(gc);
GYbar=GY(GCbar);
CYbar=CY(GCbar);
GC=@(GY)GY./(1-GY);
Gbar=GYbar.*Ybar;
Cbar=CYbar.*Ybar;
r=(Mbar.*epsilon.*CYbar)./(1-Mbar.*GYbar);
mbar=Mbar.*(1-ubar)./(1-GYbar.*Mbar);

q=@(x)omega.*x.^(-eta);
f=@(x)omega.*x.^(1-eta);
u=@(x)s./(s+f(x));
Y=@(x)1-u(x);
taubar=(1-eta).*ubar./eta;
rho=q(xbar)./s.*taubar./(1+taubar);
tau=@(x)s.*rho./(q(x)-s.*rho);
dlnydlnx=@(x)(1-eta).*u(x)-eta.*tau(x);

% General CES utility function <> the derivatives of U are computed in online appendix A

scalar=(1-gamma).^(1-gamma).*gamma.^gamma;
if epsilon==1
	U=@(c,g)c.^(1-gamma).*g.^(gamma)./scalar;
else
	U=@(c,g)((1-gamma).^(1./epsilon).*c.^((epsilon-1)./epsilon)+gamma.^(1./epsilon).*g.^((epsilon-1)./epsilon)).^(epsilon./(epsilon-1));
end
dUdc=@(gc)((1-gamma).*U(1,gc)).^(1./epsilon);
dUdcbar=dUdc(GCbar);
dUdg=@(gc)(gamma.*U(1./gc,1)).^(1./epsilon);
MRS=@(gc)gamma.^(1./epsilon)./(1-gamma).^(1./epsilon).*gc.^(-1./epsilon);
dlnUdlnc=@(gc)(1-gamma).^(1./epsilon).*(U(1,gc)).^((1-epsilon)./epsilon);
dlnUdlng=@(gc)gamma.^(1./epsilon).*(U(1./gc,1)).^((1-epsilon)./epsilon);
dlnUcdlnc=@(gc)(dlnUdlnc(gc)-1)./epsilon;
dlnUcdlng=@(gc)dlnUdlng(gc)./epsilon;

p0=dUdcbar.^r./(1+taubar);
p=@(G)p0.*dUdc(G./(Ybar-G)).^(1-r);
dlnpdlng=@(G)(1-r).*(dlnUcdlng(G./(Ybar-G))-dlnUcdlnc(G./(Ybar-G)).*(G./(Ybar-G)));% equation (A4)
dlncdlnx=@(G,x)eta.*tau(x)./dlnUcdlnc(G./(Y(x)-G));% equation (A6)
dlncdlng=@(G,x)(dlnpdlng(G)-dlnUcdlng(G./(Y(x)-G)))./dlnUcdlnc(G./(Y(x)-G));% equation (A7)
dlnxdlng=@(G,x)(G./Y(x)+(1-G./Y(x)).*dlncdlng(G,x))./(dlnydlnx(x)-(1-G./Y(x)).*dlncdlnx(G,x));% equation (A9)
m=@(G,x)(1-eta).*u(x).*(1-u(x)).*dlnxdlng(G,x).*Y(x)./G;% equation (A11)
M=@(G,x)m(G,x)./(1-u(x)+G./Y(x).*eta.*tau(x)./(1-eta)./u(x).*m(G,x));% equation (A11)

findx=@(G,x,alpha)abs(dUdc(G./(Y(x)-G))-((1+tau(x)).*p(G)./alpha)); % gap between aggregate supply & aggregate demand
optimal=@(G,x)abs(1-MRS(G./(Y(x)-G))-dlnydlnx(x).*dlnxdlng(G,x).*Y(x)./G); % first-order condition (18)
z1=GYbar.*CYbar./ubar;
z0=1./((1-eta).*(1-ubar)^2);
suffstat=@(G,x)epsilon.*z0.*m(G,x)./(1+z1.*z0.*epsilon.*m(G,x).^2);% sufficient-statistic formula (23) 

%% Business-cycle simulations under aggregate-demand shocks <> see section 5

% First step: computing equilibrium public expenditure G and equilibrium tightness x for a range of aggregate demands
% The computation is repeated for three public-expenditure policies: 
% (1) G/Y=16.5%
% (2) G/Y from sufficient-statistic formula (23) 
% (3) optimal G/Y

j=0;
ALPHA=[0.97:0.0025:1.03]; %range of aggregate demand
x0=[0.001:0.0002:2]; %grid to search for equilibrium tightness x
GY0=[0.07:0.00005:0.27]; %grid to search for equilibrium public expenditure G/Y
[GY1,x1]=ndgrid(GY0,x0);

for alpha=ALPHA
	
	j=j+1;
	
	% First policy: G/Y=16.5%

	G0=GYbar.*Y(x0);	% G such that G/Y=16.5%
	eva=findx(G0,x0,alpha);
	[val,ind]=min(eva); % finds x such aggregate supply = aggregate demand
	xconstant(j)=x0(ind); % equilibrium x 
	Gconstant(j)=G0(ind); % equilibrium G 

  % Second policy: G/Y from formula (23) 

	du=(u(xconstant(j))-ubar)./ubar;
	suffstat0=suffstat(Gconstant(j),xconstant(j)).*du; % computes G/Y with formula (23) from initial steady state
	G0=GY((1+suffstat0).*GCbar).*Y(x0); % corresponding G
	eva=findx(G0,x0,alpha);
	[val,ind]=min(eva); % finds x such aggregate supply = aggregate demand
	xsuffstat(j)=x0(ind); % equilibrium x 
	Gsuffstat(j)=G0(ind);  % equilibrium G 
	
	% Third policy: optimal G/Y

	G1=GY1.*Y(x1);
	eva=findx(G1,x1,alpha);
	[val,ind]=min(eva,[],2); % finds x such aggregate supply = aggregate demand, for all G/Y in GY0
	x2=x0(ind);
	G2=GY0.*Y(x2); 
	
	eva=optimal(G2,x2);
	[val,ind]=min(eva); % among all G/Y in GY0, finds G/Y satisfying condition (18) 
	xoptimal(j)=x2(ind); % equilibrium x 
	Goptimal(j)=G2(ind);  % equilibrium G 

end

clear GY1 x1 G1 G0 x2 G2 Y2 c2 g2

% Second step: computing all other equilibrium variables
% Equilibrium with G/Y=16.5%

Yconstant=Y(xconstant);
GYconstant=Gconstant./Yconstant;
uconstant=u(xconstant);
Mconstant=M(Gconstant,xconstant);

% Equilibrium with G/Y from formula (23) 

Ysuffstat=Y(xsuffstat);
GYsuffstat=Gsuffstat./Ysuffstat;
usuffstat=u(xsuffstat);
Msuffstat=M(Gsuffstat,xsuffstat);

% Equilibrium with optimal G/Y

Yoptimal=Y(xoptimal);
GYoptimal=Goptimal./Yoptimal;
uoptimal=u(xoptimal);
Moptimal=M(Goptimal,xoptimal);

%% Plot figure 4

% Formatting

set(groot,'DefaultFigureUnits', 'inches')
set(groot,'DefaultFigurePosition', [0,0,7.7778,5.8333]);
set(groot,'DefaultFigurePaperPosition', [0, 0, 7.7779,5.8334]);
set(groot,'DefaultFigurePaperSize', [7.7779,5.8334]);
set(groot,'DefaultAxesFontName', 'Times New Roman')
set(groot,'DefaultAxesFontSize', 25)
set(groot,'DefaultLineLineWidth', 4)
set(groot,'DefaultAxesLineWidth', 1)
set(groot,'DefaultAxesXColor','k')
set(groot,'DefaultAxesYColor','k')
set(groot,'DefaultAxesYGrid','on')
set(groot,'DefaultAxesXGrid','off')
set(groot,'DefaultAxesTickLength',[0 0])
co = [217,95,2
		  8,81,156
			107,174,214];
co = co./255;
set(groot,'DefaultAxesColorOrder',co)

xmi=min(ALPHA);xma=max(ALPHA);
fign=0;

% Figure 4A: unemployment rate

x1=ALPHA;
y1=uoptimal.*100;
y2=uconstant.*100;

fign=fign+1;
figure(fign)
clf
hold on
plot(x1,y2,'-.','Color',co(1,:))
plot(x1,y1,'-','Color',co(2,:))
xlabel('Aggregate demand')
set(gca,'XLim',[xmi,xma],'XTick',[0.97,1,1.03])
set(gca,'YLim',[4,12],'YTick',[4:2:12],'YTickLabel',[' 4%';' 6%';' 8%';'10%';'12%'])
ylabel('Unemployment rate')
print('-dpdf','figure4A.pdf')

% Figure 4B: unemployment multiplier

x1=ALPHA;
y1=Moptimal;
y2=Mconstant;

fign=fign+1;
figure(fign)
clf
hold on
plot(x1,y2,'-.','Color',co(1,:))
plot(x1,y1,'-','Color',co(2,:))
xlabel('Aggregate demand')
set(gca,'XLim',[xmi,xma],'XTick',[0.97,1,1.03])
set(gca,'YLim',[0,1.5],'YTick',[0:0.5:1.5])
ylabel('Unemployment multiplier')
print('-dpdf','figure4B.pdf')

% Figure 4C: optimal stimulus spending

x1=ALPHA;
y1=GYoptimal.*100;
y2=GYconstant.*100;

fign=fign+1;
figure(fign)
clf
hold on
plot(x1,y2,'-.','Color',co(1,:))
plot(x1,y1,'-','Color',co(2,:))
xlabel('Aggregate demand')
set(gca,'XLim',[xmi,xma],'XTick',[0.97,1,1.03])
set(gca,'YLim',[13,21],'YTick',[13:2:21],'YTickLabel',['13%';'15%';'17%';'19%';'21%'])
ylabel('Public expenditure (% of GDP)')
print('-dpdf','figure4C.pdf')

% Figure 4D: accuracy of sufficient-statistic formula

x1=ALPHA;
y2=GYsuffstat.*100;
y1=GYoptimal.*100;

fign=fign+1;
figure(fign)
clf
hold on
plot(x1,y2,'-.','Color',co(3,:))
plot(x1,y1,'-','Color',co(2,:))
xlabel('Aggregate demand')
set(gca,'XLim',[xmi,xma],'XTick',[0.97,1,1.03])
set(gca,'YLim',[13,21],'YTick',[13:2:21],'YTickLabel',['13%';'15%';'17%';'19%';'21%'])
ylabel('Public expenditure (% of GDP)')
print('-dpdf','figure4D.pdf')