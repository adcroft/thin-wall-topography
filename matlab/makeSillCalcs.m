function [] = makeSillCalcs(reg)

switch reg
 case 'denmark'
  x=[-35 -20]; y=[62 70]; x0=[-30 -24]; y0=[64.5 68];
 case 'gibraltar'
  x=[-10 -2]; y=[33 38]; x0=[-7.5 -4.5]; y0=[35.5 36];
 case 'faroe'
  x=[-15 0]; y=[55 67]; x0=[-1 -14]; y0=[65 61];
 case 'romanche'
  x=[-28 -8]; y=[-2 3]; x0=[-15 -27]; y0=[2 -1.8];
 case 'vema'
  x=[-44 -34]; y=[-35 -24]; x0=[-34.7 -39.5]; y0=[-25.5 -34.5];
 case 'lombok'
  x=[113 119]; y=[-12 -6]; x0=[115.47 116.06]; y0=[-9.93 -7.75];
 case 'redsea'
  x=[40 48]; y=[11 16]; x0=[41.612087139845400 47.495432185523548]; y0=[15.955580865603652 11.695899772209568];
 case 'ALL'
  makeSillCalcs('denmark');
  makeSillCalcs('gibraltar');
  makeSillCalcs('faroe');
  makeSillCalcs('romanche');
  makeSillCalcs('vema');
  makeSillCalcs('lombok');
  return
 otherwise
  error('Unknown "reg"')
end

%load ESM2.5G_v2.mat
%dothincalc(reg,G,x0,y0)
[D_etopo2_real,res,D_etopo2_thin,D_etopo2_step] = docalc_etopo2(reg,x,y,x0,y0);
[D_etopo1_real,res,D_etopo1_thin,D_etopo1_ave] = docalc_etopo1(reg,x,y,x0,y0);
[D_sandwell_real,res,D_sandwell_thin,D_sandwell_ave] = docalc_sandwell(reg,x,y,x0,y0);
[D_gebco_real,res,D_gebco_thin,D_gebco_ave] = docalc_gebco08(reg,x,y,x0,y0);
eval(sprintf('save %s_sill_depths.mat',reg));

% =======================================================

function [D_real,res,D_thin,D_ave] = docalc_sandwell(reg,x,y,x0,y0)

[He,xe,ye]=read_sandwell(x+[-1 1],y+[-1 1]);
[D_real,res,D_thin,D_ave] = doallcalcs(reg,He,xe,ye,x,y,x0,y0)

% -------------------------------------------------------

function [D_real,res,D_thin,D_ave] = docalc_etopo2(reg,x,y,x0,y0)

[He,xe,ye]=read_etopo2(x+[-1 1],y+[-1 1]);
[D_real,res,D_thin,D_ave] = doallcalcs(reg,He,xe,ye,x,y,x0,y0)

% -------------------------------------------------------

function [D_real,res,D_thin,D_ave] = docalc_etopo1(reg,x,y,x0,y0)

[He,xe,ye]=read_etopo1(x+[-1 1],y+[-1 1]);
[D_real,res,D_thin,D_ave] = doallcalcs(reg,He,xe,ye,x,y,x0,y0)

% -------------------------------------------------------

function [D_real,res,D_thin,D_ave] = docalc_gebco08(reg,x,y,x0,y0)

[He,xe,ye]=read_gebco_08(x+[-1 1],y+[-1 1]);
[D_real,res,D_thin,D_ave] = doallcalcs(reg,He,xe,ye,x,y,x0,y0)

% -------------------------------------------------------

function [D_real,res,D_thin,D_ave] = doallcalcs(reg,He,xe,ye,x,y,x0,y0)

D_real=donormcalc(reg,He,xe,ye,x0,y0)

res=[];D_thin=[];D_ave=[];

res(end+1)=0.1;
[G,xc,yc]=lgrid(x,y,res(end),He,xe,ye); D_ave(end+1)=donormcalc(reg,G.H.aave,xc,yc,x0,y0); D_thin(end+1,:)=dothincalc(reg,G,x0,y0)

res(end+1)=0.125;
[G,xc,yc]=lgrid(x,y,res(end),He,xe,ye); D_ave(end+1)=donormcalc(reg,G.H.aave,xc,yc,x0,y0); D_thin(end+1,:)=dothincalc(reg,G,x0,y0)

res(end+1)=0.25;
[G,xc,yc]=lgrid(x,y,res(end),He,xe,ye); D_ave(end+1)=donormcalc(reg,G.H.aave,xc,yc,x0,y0); D_thin(end+1,:)=dothincalc(reg,G,x0,y0)

res(end+1)=0.5;
[G,xc,yc]=lgrid(x,y,res(end),He,xe,ye); D_ave(end+1)=donormcalc(reg,G.H.aave,xc,yc,x0,y0); D_thin(end+1,:)=dothincalc(reg,G,x0,y0)

res(end+1)=1;
[G,xc,yc]=lgrid(x,y,res(end),He,xe,ye); D_ave(end+1)=donormcalc(reg,G.H.aave,xc,yc,x0,y0); D_thin(end+1,:)=dothincalc(reg,G,x0,y0)

% -------------------------------------------------------

function [sill_real] = donormcalc(reg,He,xe,ye,x0,y0)

He(isnan(He))=max(He(:)); % NaN's appear at edges due to extrapolation
pcol(xe,ye,He');colorbar;hold on;plot(x0,y0,'w*');hold off;drawnow
[XG,YG] = qgrid(xe,ye);
sill_real=sill_depth(He,[],[],XG,YG,x0,y0);

% -------------------------------------------------------

function [sills] = dothincalc(reg,G,x0,y0)

plotraster(G.X,G.Y,G.H.emin,G.U.emin,G.V.emin,[-5000 2000]);colorbar; hold on;plot(x0,y0,'w*');hold off;drawnow
sill_emin=sill_depth(G.H.emin,G.U.emin,G.V.emin,G.X,G.Y,x0,y0);
sill_ave=sill_depth(G.H.amin,G.U.aave,G.V.aave,G.X,G.Y,x0,y0);
sill_eave=sill_depth(G.H.emin,G.U.eave,G.V.eave,G.X,G.Y,x0,y0);
sills=[sill_emin,sill_ave,sill_eave];

% -------------------------------------------------------

function [XG,YG] = qgrid(xe,ye)
% Approximate corner/grid coordinates
xg=xe-(xe(2)-xe(1))/2;xg(end+1)=xg(end)+(xe(2)-xe(1));
yg=ye-(ye(2)-ye(1))/2;yg(end+1)=yg(end)+(ye(2)-ye(1));
[XG,YG]=ndgrid(xg,yg);

% -------------------------------------------------------

function [G,xc,yc] = lgrid(x,y,res,He,xe,ye)
xg=min(x(:)):res:max(x(:));
yg=min(y(:)):res:max(y(:));
xc=min(x(:)+res/2):res:max(x(:));
yc=min(y(:)+res/2):res:max(y(:));
[XG,YG] = ndgrid(xg,yg);
G=thinwallGenerator(XG,YG,He,xe,ye,'fv6');
