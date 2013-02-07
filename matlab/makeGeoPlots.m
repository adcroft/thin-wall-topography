function [] = makeGeoPlots(reg)

switch reg
 case 'test'
  x=[-80 -60];y=[-75 -65];
  doplots(x,y,reg);
 case 'indo'
  x=[80 168];y=[-25 25];
  x=[100 140];y=[-15 7];
  doplots(x,y,reg);
 case 'indo2'
  x=[113 135];y=[-13 1];
  doplots(x,y,reg);
 case 'indo3'
  x=[117 128];y=[-12 -5];
  doplots(x,y,reg);
 case 'drake'
  x=[-80 -20];y=[-75 -40];
  doplots(x,y,reg);
 case 'central_america'
  x=[-100 -55];y=[5 32];
  doplots(x,y,reg);
 case 'panama'
  x=[-85 -75];y=[6 11];
  doplots(x,y,reg);
 case 'gibraltar'
  x=[-10 -2];y=[33 38];
  doplots(x,y,reg);
 case 'gibraltar2'
  x=[-9.5 -1.5];y=[33.5 38.5];
  doplots(x,y,reg);
 case 'iceland'
  x=[-35 0];y=[57 72];
  doplots(x,y,reg);
 case 'midatl'
  x=[-35 5];y=[-20 7];
  doplots(x,y,reg);
 case 'queen_elizabeth_islands'
  x=[-130 -55];y=[60 85];
  doplots(x,y,reg);
 case 'ALL'
  makeGeoPlots('indo')
  makeGeoPlots('drake')
  makeGeoPlots('central_america')
  makeGeoPlots('gibraltar')
  makeGeoPlots('iceland')
  makeGeoPlots('midatl')
  makeGeoPlots('queen_elizabeth_islands')
 otherwise
  error('Unknown "reg"')
end

% =======================================================

function [] = doplots(x,y,reg)
Hrng=[-5000 3000];

aspect=(max(y(:))-min(y(:)))/(max(x(:))-min(x(:)));
aspect=.6;
pwidth=4;

figure(1);clf
set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 1 aspect]*pwidth)
[He2,xe2,ye2]=read_sandwell(x,y);
Hrng=drng(He2);
cm_landwater(400,Hrng);
clf;pcolor(xe2,ye2,He2');shading flat;caxis(Hrng);
xlim(x);ylim(y);
labelStraits; drawnow
set(gca,'Position',[.05 .06 .925 .925])
printit(reg,[],'sandwell')
if 0

%ddH=diff(Hrng)/length(colormap)/10;
%HHx=(min(Hrng)-diff(Hrng)/400/2):50:max(Hrng);
%HHx=(0:10*length(colormap))*ddH+min(Hrng);
%clf;pcol(HHx,[0 1],[HHx;HHx]); caxis(Hrng); 
%set(gca,'Position',[0.05 0.05 0.9 0.05],'YTick',[])
%for l=-6000:1000:6000
%  line([1 1]*l,[-1 2],'Color','k');
%end
%printit(reg,[],'colorbar')

[He2,xe2,ye2]=read_etopo1(x,y);
clf;pcolor(xe2,ye2,He2');shading flat;caxis(Hrng);
xlim(x);ylim(y);
labelStraits; drawnow
set(gca,'Position',[.05 .06 .925 .925])
printit(reg,[],'etopo1')

[He2,xe2,ye2]=read_etopo2(x,y);
clf;pcolor(xe2,ye2,He2');shading flat;caxis(Hrng);
xlim(x);ylim(y);
labelStraits; drawnow
set(gca,'Position',[.05 .06 .925 .925])
printit(reg,[],'etopo2')

[He2,xe2,ye2]=read_gebco_08(x,y);
clf;pcolor(xe2,ye2,He2');shading flat;caxis(Hrng);
xlim(x);ylim(y);
labelStraits; drawnow
set(gca,'Position',[.05 .06 .925 .925])
printit(reg,[],'gebco08')

%interpcolormap('bwr','w');
%drawColorbar([0 1],01:.1:1)
%printit(reg,[],'colorbar_bwr')

%interpcolormap('rwb');
%drawColorbar([0 1],01:.1:1)
%printit(reg,[],'colorbar_rwb')

interpcolormap('rjet');
drawColorbar([0 1],01:.1:1)
printit(reg,[],'colorbar_jet')
end
figure(10)
cm_landwater(400,Hrng);
drawColorbar(Hrng,-7000:1000:7000)
printit(reg,[],'colorbar')


for res=[1 .5 .25]
 fname=sprintf('Results/%s_%g.mat',reg,res);
 h=2*res;
 if exist(fname,'file')
  G=load(fname)
 else
  xg=x(1)-h:res:x(2)+h;
  yg=y(1)-h:res:y(2)+h;
  [XG,YG]=ndgrid(xg,yg);
  G=thinwallGenerator(XG,YG,'fv6');
  save(fname,'-struct','G')
 end

%[Hm,Um,Vm,Hf,Uf,Vf,Hmin,Hmax,Have,Umin,Umax,Uave,Vmin,Vmax,Vave,Hz,Uz,Vz,Hfi,Ufi,Vfi]=diagWidthOfChannel(G);
 
 figure(2);clf
 set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 1 aspect]*pwidth)
 cm_landwater(400,Hrng);
 %pcolor(G.X,G.Y,G.H.aave([1:end end],[1:end end]));shading flat;caxis(Hrng);
 pcol(G.X,G.Y,G.H.aave);shading flat;caxis(Hrng);
 hold on;contour(xe2,ye2,He2',[-inf 0],'k');hold off
 xlim(x);ylim(y);drawnow
 set(gca,'Position',[.05 .06 .925 .925])
 printit(reg,res,'aave')

 %plotraster(G.X,G.Y,G.H.aave,G.U.aave,G.V.aave,Hrng)
 %xlim(x);ylim(y);drawnow
 %printit(reg,res,'uvave')

 plotraster(G.X,G.Y,G.H.emin,G.U.emin,G.V.emin,Hrng)
%plotraster(G.X,G.Y,Hmin,Umin,Vmin,Hrng)
 hold on;contour(xe2,ye2,He2',[-inf 0],'k');hold off
 xlim(x);ylim(y);drawnow
 set(gca,'Position',[.05 .06 .925 .925])
 printit(reg,res,'emin')

 %Hpot=min(G.U.emin(1:end-1,:),G.U.emin(2:end,:)); Hpot=min(Hpot,min(G.V.emin(:,1:end-1),G.V.emin(:,2:end)));
 %plotraster(G.X,G.Y,Hpot,G.U.emin,G.V.emin,Hrng)
 %xlim(x);ylim(y);drawnow
 %printit(reg,res,'eminpot')

 figure(3);clf
 set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 1 aspect]*pwidth)
 cm_landwater(400,Hrng);
 plotraster(G.X,G.Y,G.H.aave,G.U.eaave,G.V.eaave,Hrng)
%plotraster(G.X,G.Y,Have,Uave,Vave,Hrng)
 hold on;contour(xe2,ye2,He2',[-inf 0],'k');hold off
 xlim(x);ylim(y);drawnow
 set(gca,'Position',[.05 .06 .925 .925])
 printit(reg,res,'eaave')

 %Hpot=min(G.U.eave(1:end-1,:),G.U.eave(2:end,:)); Hpot=min(Hpot,min(G.V.eave(:,1:end-1),G.V.eave(:,2:end)));
 %plotraster(G.X,G.Y,Hpot,G.U.eave,G.V.eave,Hrng)
 %xlim(x);ylim(y);drawnow
 %printit(reg,res,'eavepot')

 figure(4);clf
 set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 1 aspect]*pwidth)
 cm_landwater(400,Hrng);
 plotraster(G.X,G.Y,G.H.amax,G.U.emax,G.V.emax,Hrng)
%plotraster(G.X,G.Y,Hmax,Umax,Vmax,Hrng)
 hold on;contour(xe2,ye2,He2',[-inf 0],'k');hold off
 xlim(x);ylim(y);drawnow
 set(gca,'Position',[.05 .06 .925 .925])
 printit(reg,res,'emax')

if 1==0
 figure(5);clf
 set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 1 aspect]*pwidth)
%Um=(G.U.aave-G.U.amin)./(G.U.amax-G.U.amin);
%Vm=(G.V.aave-G.V.amin)./(G.V.amax-G.V.amin);
%Hm=(G.H.aave-G.H.emin)./(G.H.amax-G.H.emin);
 plotraster(G.X,G.Y,Hm,Um,Vm,[-1e-13 1]);interpcolormap('bwr','w')
 hold on;contour(xe2,ye2,He2',[-inf 0],'k');hold off
 xlim(x);ylim(y);drawnow
 set(gca,'Position',[.05 .06 .925 .925])
%printit(reg,res,'m')
colorbar;title('m')

 figure(6);clf
 set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 1 aspect]*pwidth)
 plotraster(G.X,G.Y,Hf,Uf,Vf,[-1e-13 1]);interpcolormap('rjet')
 hold on;contour(xe2,ye2,He2',[-inf 0],'w');hold off
 xlim(x);ylim(y);drawnow
 set(gca,'Position',[.05 .06 .925 .925])
%printit(reg,res,'f')
colorbar;title('f')

 figure(7);clf
 set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 1 aspect]*pwidth)
 plotraster(G.X,G.Y,Hz,Uz,Vz);colormap(jet)
 hold on;contour(xe2,ye2,He2',[-inf 0],'w');hold off
 xlim(x);ylim(y);drawnow
 set(gca,'Position',[.05 .06 .925 .925])
%printit(reg,res,'z')
colorbar;title('z')

 figure(8);clf
 set(gcf,'PaperUnits','Inches','PaperPosition',[0 0 1 aspect]*pwidth)
 set(gca,'Position',[.05 .05 .95 .95])
%cm_landwater(400,Hrng); clf;
 plotraster(G.X,G.Y,Hfi,Ufi,Vfi);interpcolormap('rjet')
 set(gca,'Position',[0.05 0.05 .9 .9])
 hold on;contour(xe2,ye2,He2',[-inf 0],'w');hold off
 xlim(x);ylim(y);drawnow
 set(gca,'Position',[.05 .06 .925 .925])
%printit(reg,res,'fi')
colorbar;title('fi')
end % 1==0

%pause

%Hpot=min(G.U.emax(1:end-1,:),G.U.emax(2:end,:)); Hpot=min(Hpot,min(G.V.emax(:,1:end-1),G.V.emax(:,2:end)));
%plotraster(G.X,G.Y,Hpot,G.U.emax,G.V.emax,Hrng)
%xlim(x);ylim(y);drawnow
%printit(reg,res,'emaxpot')

%plotraster(G.X,G.Y,G.H.ave,G.U.eave,G.V.eave,Hrng)
%xlim(x);ylim(y);drawnow
%printit(reg,res,'vizeave')

end

% =======================================================

function [] = printit(reg,res,lbl)
set(gca,'FontSize',6)
if ~isempty(res)
 sres=sprintf('_%g',res);
else
 sres='';
end
fname=sprintf('Images/%s%s_%s',reg,sres,lbl);
print([fname '.png'],'-dpng','-r300')
%print([fname '.eps'],'-depsc2','-r300')
%print([fname '.tif'],'-dtiff','-r300')

% =======================================================

function [] = mycb(Hrng)

gax=axis;
ax=gca; pax=get(ax,'Position');

h=colorbar('Location','South');
op=get(h,'Position');
op(2)=pax(2)+pax(4)*1.02;
dx=op(3); op(1)=op(1)+.2*dx; op(3)=.6*dx;

op=[.05 .05 .95 .95];
set(h,'Position',op);

% =======================================================

function [Hrng] = drng(H)
H=H(:);
[n,d]=hist(H,100);n=n/numel(H);
N=cumsum(n);
hmin=d( find(N>0.01,1,'first') );
hmax=d( find(N<0.99,1,'last') );
hmax=max(hmax,1000);
Hrng=round([hmin hmax]/10)*10;
