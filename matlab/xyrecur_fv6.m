function [GS]=xyrecur_fv6(Hin,XG,YG,niter,varargin)
% G=xyrecur_fv6(Hin,XG,YG,niter,options)
%
% Applies niter recursions of the "thin wall" algorithm (see Adcroft, OM 2013).
%
% e.g.
%      [G]=xyrecur_fv6(FHmin,FG.X,FG.Y,3,'fv6');

% Copyright 2008-2013 Alistair Adcroft, Princeton University.
%
% This file is part of the thin-wall-topography software suite.
%
% thin-wall-topography is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% thin-wall-topography is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with Foobar. If not, see http://www.gnu.org/licenses/.


txyr=tic;
% Setup methods from varargin
global viznorm vizdebug Hrng Xorig Yorig Horig doColorbar cm lowMem verbose timing
Horig=Hin;
%Horig([1 end],:)=Inf; Horig(:,[1 end])=Inf; % testing impact of boundaries
viznorm=0; vizdebug=0; doColorbar=0;
locdisp(sprintf('xyrecur_fv6: Range in input topog: %f : %f',min(Hin(:)),max(Hin(:))))
Control.doLowerButress=1;
Control.doCorners=1;
Control.doCentral=1;
Control.doExpandCorners=1;
Control.iterate=1;
Control.doConnections=1;
Control.doBoundH=1;
Control.doPotholes=0;
Control.doGenHxHy=1;
Control.showAverage=0;
lowMem=0;
verbose=1;
timing=1;

global iDump jDump widthDump
%iDump=8;jDump=6;widthDump=12;
global xSill ySill
%xSill=[-30 -24]; ySill=[64.5 68];

[tok,meth]=strtok(varargin{1});
if ~strcmp(tok,'fv6')
 error('xyrecur_fv6: 1st arg should be fv6')
end
while length(meth)
 [tok,meth]=strtok(meth);
 switch tok
   case {'buttress'}
    Control.doLowerButress=1-Control.doLowerButress;
  case {'corners'}
    Control.doCorners=1-Control.doCorners;
  case {'central'}
    Control.doCentral=1-Control.doCentral;
  case {'expandcorners'}
    Control.doExpandCorners=1-Control.doExpandCorners;
  case {'pot'}
    Control.doPotholes=1-Control.doPotholes;
  case {'boundh'}
    Control.doBoundH=1-Control.doBoundH;
  case {'connections'}
    Control.doConnections=1-Control.doConnections;
  case {'iter'}
    Control.iterate=1-Control.iterate;
  case {'regen'}
    Control.doGenHxHy=1-Control.doGenHxHy;
  case {'showave'}
    Control.showAverage=1-Control.showAverage;
  case {'lowmem'}
    lowMem=1;
  case {'debug'}
    vizdebug=1-vizdebug;
  case {'viz'}
    viznorm=1-viznorm;
  case {'verbose'}
    verbose=1-verbose;
  case {'timing'}
    timing=1-timing;
  otherwise
    tok, meth
    error('xyrecur_fv6: Unknown method')
  end
end

%debugging [nx,ny]=size(XG);xg=1:nx;yg=1:ny;[XG,YG]=ndgrid(xg,yg);

hres=50;Hrng=[hres*floor(min(Hin(:)/hres)) hres*ceil(max(Hin(:)/hres))];
%gc=gcd(abs(min(Hrng)),max(Hrng)); ncm=factor_products( gc );
%ncm=ncm(:,1);ncm(ncm>300)=[];ncm=max(ncm);
if viznorm | vizdebug
 ncm=256; figure(1); cm=cm_landwater(ncm,Hrng);
%cm=interpcolormap(128,'bgr');
%cm=jet;colormap(cm);
 locdisp(sprintf('xyrecur_fv6: Range in colorbar: %f : %f ,%i',Hrng,ncm))
end

GS.Control=Control; if verbose; Control, end

if vizdebug | viznorm
 Xorig=XG;Yorig=YG; % Keep for plotting later
end

% Repeatedly coarsen grid
G.X=XG; G.Y=YG;
H.emin=Horig; U.emin=[]; V.emin=[]; 
timc=[];
for nit=1:niter,
  tim=tic;
  [G,H,U,V]= coarsen(Control,G,H,U,V,nit);
  timc(end+1)=toc(tim); locdisp([' coarsen: timings' sprintf(' %4.2fs',timc)]);

  if vizdebug & nit<niter
    spn=gcf;
    inp='x';
    while( ~strcmp(inp,'') )
      figure(1);zoom on
%     inp=input('Waiting ... hit enter to continue (z to zoom)','s');
      inp=''
      if strcmp(inp,'z')
        figure(1);ax=axis;zoom off;zoomallfigs(ax,spn)
      end
      figure(1);zoom off
    end
  end
end

GS.X=G.X; GS.Y=G.Y; GS.H=H; GS.U=U; GS.V=V;
txyr=toc(txyr);
if timing
 fprintf('xyrecur_fv6 took %f secs\n',txyr);
end

% -----------------------------------------------------------------------------

function [CG,Hc,Uc,Vc]=coarsen(Control,FG,Hf,Uf,Vf,lvl)
global viznorm vizdebug Hrng Xorig Yorig Horig doColorbar
spn=0; % subplot number

CG.X=FG.X(1:2:end,1:2:end); CG.Y=FG.Y(1:2:end,1:2:end);

[nxf,nyf]=size(Hf.emin);
ic=1:nxf/2;jc=1:nyf/2;
iw=1:2:nxf; ie=2:2:nxf; js=1:2:nyf; jn=2:2:nyf; % Cell indices
Iu=1:2:nxf+1;Jv=1:2:nyf+1; % Edge indices
iu=2:2:nxf; jv=2:2:nyf; % Edge indices in middle of cell
locdisp(sprintf(' coarsen: level %i (%i,%i) -> (%i,%i) ... ',lvl,nxf,nyf,nxf/2,nyf/2));

% Initialize U,V with step topography
if isempty(Uf.emin)
 locdisp('  Umin,Vmin generated from cell centers on entry ...')
 [Uf.emin,Vf.emin]=gen_hxhy(Hf.emin,[],[]);
end
% Allocate memory for ave, max, min, eff
if ~isfield(Hf,'aave')
 Hf.aave=Hf.emin; Hf.amax=Hf.emin; Hf.amin=Hf.emin;
 Uf.lave=Uf.emin; Uf.lmax=Uf.emin; Uf.lmin=Uf.emin;
 Vf.lave=Vf.emin; Vf.lmax=Vf.emin; Vf.lmin=Vf.emin;
 Hf.eave=Hf.emin; Uf.eave=Uf.emin; Vf.eave=Vf.emin;
 Uf.emax=Uf.emin; Vf.emax=Vf.emin;
 Uf.eaave=Uf.emin; Vf.eaave=Vf.emin;
end

if vizdebug | viznorm ;spn=spn+1;
 myfig(spn); myplot(Xorig,Yorig,Horig,[],[],Hrng);plotcolorbar
 overlaygrid(CG.X,CG.Y)
 title(sprintf('Original input data (grid at level %i)',lvl));drawnow
end

if vizdebug ;spn=spn+1;
 myfig(spn); myplot(FG.X,FG.Y,Hf.emin,Uf.emin,Vf.emin,Hrng);plotcolorbar
 overlaygrid(FG.X(1:2:end,1:2:end),FG.Y(1:2:end,1:2:end))
 title(sprintf('Fine grid: Level %i minimum',lvl-1));drawnow
end
checksill('Fine grid',Hf,Uf,Vf,FG)

OriginalDeepestConnections = findDeepestConnections(Uf.emin,Vf.emin);
dump4Cells('Fine grid',Uf.emin,Vf.emin,OriginalDeepestConnections)

if vizdebug ;spn=spn+1;
 myfig(spn); myplot(FG.X,FG.Y,Hf.eave,Uf.eave,Vf.eave,Hrng);plotcolorbar
 overlaygrid(FG.X(1:2:end,1:2:end),FG.Y(1:2:end,1:2:end))
 title(sprintf('Fine grid: Level %i average',lvl-1));drawnow
end

% Push out corners made by ridges
if Control.doCorners
 [Hf,Uf,Vf,msk] = pushCorners(Hf,Uf,Vf); msk2=msk;
 while ~isnan(max(msk2(:))) & Control.iterate
  [Hf,Uf,Vf,msk2] = pushCorners(Hf,Uf,Vf); msk=max(msk,msk2);
  warn(msk2,'pushCorners')
 end
 checksill('pushCorners',Hf,Uf,Vf,FG)
 if vizdebug ;spn=spn+1;
  myfig(spn); myplot(FG.X,FG.Y,Hf.emin,Uf.emin,Vf.emin,Hrng);plotcolorbar
  overlaygrid(FG.X(1:2:end,1:2:end),FG.Y(1:2:end,1:2:end))
  showChangedCells(FG.X,FG.Y,msk);
  title(sprintf('Fine grid after pushing out corner ridges: Level %i minimum',lvl-1));drawnow
 end
 checkDeepestConnections(OriginalDeepestConnections,Uf.emin,Vf.emin,'pushCorners');
end

% Lower tallest buttresses
if Control.doLowerButress
 [Uf.emin,Vf.emin,msk] = lowerTallestButtress(Uf.emin,Vf.emin); msk2=msk;
 while ~isnan(max(msk2(:))) & Control.iterate
  [Uf.emin,Vf.emin,msk2] = lowerTallestButtress(Uf.emin,Vf.emin,'min'); msk=max(msk,msk2);
  warn(msk2,'lowerTallestButtress')
 end
 checksill('lowerTallestButtress',Hf,Uf,Vf,FG)
 if vizdebug ;spn=spn+1;
  myfig(spn); myplot(FG.X,FG.Y,Hf.emin,Uf.emin,Vf.emin,Hrng);plotcolorbar
  overlaygrid(FG.X(1:2:end,1:2:end),FG.Y(1:2:end,1:2:end))
  showChangedCells(FG.X,FG.Y,msk);
  title(sprintf('Fine grid after lowering buttresses: Level %i minimum',lvl-1));drawnow
 end
 checkDeepestConnections(OriginalDeepestConnections,Uf.emin,Vf.emin,'lowerTallestButtress');
 [Uf.eave,Vf.eave,msk] = lowerTallestButtress(Uf.eave,Vf.eave); msk2=msk;
 while ~isnan(max(msk2(:))) & Control.iterate
  [Uf.eave,Vf.eave,msk2] = lowerTallestButtress(Uf.eave,Vf.eave,'eff'); msk=max(msk,msk2);
  warn(msk2,'lowerTallestButtress (effecitve)')
 end
 if vizdebug ;spn=spn+1;
  myfig(spn); myplot(FG.X,FG.Y,Hf.eave,Uf.eave,Vf.eave,Hrng);plotcolorbar
  overlaygrid(FG.X(1:2:end,1:2:end),FG.Y(1:2:end,1:2:end))
  showChangedCells(FG.X,FG.Y,msk);
  title(sprintf('Fine grid after lowering buttresses: Level %i average',lvl-1));drawnow
 end
end

% Handle central ridges
if Control.doCentral
 [Hf,Uf,Vf,msk] = centralRidge(Hf,Uf,Vf); msk2=msk;
 while ~isnan(max(msk2(:))) & Control.iterate
  [Hf,Uf,Vf,msk2] = centralRidge(Hf,Uf,Vf,'min'); msk=max(msk,msk2);
  warn(msk2,'centralRidge')
 end
 checksill('centralRidge',Hf,Uf,Vf,FG)
 if vizdebug ;spn=spn+1;
  myfig(spn); myplot(FG.X,FG.Y,Hf.emin,Uf.emin,Vf.emin,Hrng);plotcolorbar
  overlaygrid(FG.X(1:2:end,1:2:end),FG.Y(1:2:end,1:2:end))
  showChangedCells(FG.X,FG.Y,msk);
  title(sprintf('Fine grid after handling central ridges: Level %i minimum',lvl-1));drawnow
 end
 if Control.doLowerButress
  [Uf.emin,Vf.emin,msk] = lowerTallestButtress(Uf.emin,Vf.emin,'central');
  warn(msk,'centralRidge:lowerTallestButtress')
  checkDeepestConnections(OriginalDeepestConnections,Uf.emin,Vf.emin,'centralRidge');
 end
end

% Expand deepest corner to whole cell
if Control.doExpandCorners
 [Uf,Vf,msk] = expandCorners(Uf,Vf); msk2=msk;
 while ~isnan(max(msk2(:))) & Control.iterate
  [Uf,Vf,msk2] = expandCorners(Uf,Vf); msk=max(msk,msk2);
 %warn(msk2,'expandCorners')
 end
 checksill('expandCorners',Hf,Uf,Vf,FG)
 if vizdebug ;spn=spn+1;
  myfig(spn); myplot(FG.X,FG.Y,Hf.emin,Uf.emin,Vf.emin,Hrng);plotcolorbar
  overlaygrid(FG.X(1:2:end,1:2:end),FG.Y(1:2:end,1:2:end))
  showChangedCells(FG.X,FG.Y,msk);
  title(sprintf('Fine grid after expandCorners: Level %i minimum',lvl-1));drawnow
 end
 checkDeepestConnections(OriginalDeepestConnections,Uf.emin,Vf.emin,'expandCorners');
end

% Apply connections as limits
if Control.doConnections
 [Uf.emin,Vf.emin,msk] = limitNSConnections(OriginalDeepestConnections,Uf.emin,Vf.emin); msk2=msk;
 while ~isnan(max(msk2(:))) & Control.iterate
  [Uf.emin,Vf.emin,msk2] = limitNSConnections(OriginalDeepestConnections,Uf.emin,Vf.emin); msk=max(msk,msk2);
  warn(msk2,'limitNSConnections')
 end
 checksill('limitNSConnections',Hf,Uf,Vf,FG)
 if vizdebug ;spn=spn+1;
  myfig(spn); myplot(FG.X,FG.Y,Hf.emin,Uf.emin,Vf.emin,Hrng);plotcolorbar
  overlaygrid(FG.X(1:2:end,1:2:end),FG.Y(1:2:end,1:2:end))
  showChangedCells(FG.X,FG.Y,msk);
  title(sprintf('Fine grid after applying NS connections: Level %i minimum',lvl-1));drawnow
 end
 checkDeepestConnections(OriginalDeepestConnections,Uf.emin,Vf.emin,'limitNSConnections');

 [Uf.emin,Vf.emin,msk] = limitDConnections(OriginalDeepestConnections,Uf.emin,Vf.emin); msk2=msk;
 while ~isnan(max(msk2(:))) & Control.iterate
  [Uf.emin,Vf.emin,msk2] = limitDConnections(OriginalDeepestConnections,Uf.emin,Vf.emin); msk=max(msk,msk2);
 end
 checksill('limitDConnections',Hf,Uf,Vf,FG)
 if vizdebug ;spn=spn+1;
  myfig(spn); myplot(FG.X,FG.Y,Hf.emin,Uf.emin,Vf.emin,Hrng);plotcolorbar
  overlaygrid(FG.X(1:2:end,1:2:end),FG.Y(1:2:end,1:2:end))
  showChangedCells(FG.X,FG.Y,msk);
  title(sprintf('Fine grid after applying diagonal connections: Level %i minimum',lvl-1));drawnow
 end
 warn(msk2,'limitDConnections')
 checkDeepestConnections(OriginalDeepestConnections,Uf.emin,Vf.emin,'limitDConnections');
end

checkDeepestConnections(OriginalDeepestConnections,Uf.emin,Vf.emin,'Before reduction');

% Four-point average in X-Y
locdisp('  Four-point reduction for cells')
Hc.amin=reduceH('min',Hf.amin);
Hc.amax=reduceH('max',Hf.amax);
Hc.aave=reduceH('ave',Hf.aave);
Hc.emin=reduceH('min',Hf.emin);
Hc.eave=reduceH('ave',Hf.eave);

% Two-point reduction of edge data
locdisp('  Two-point reduction for CUval, CVval')
[Uc.lmin,Vc.lmin]=reduceUV('min',Uf.lmin,Vf.lmin);
[Uc.lmax,Vc.lmax]=reduceUV('max',Uf.lmax,Vf.lmax);
[Uc.lave,Vc.lave]=reduceUV('ave',Uf.lave,Vf.lave);
[Uc.amin,Vc.amin]=reduceUV2('min',Hf.amin);
[Uc.aave,Vc.aave]=reduceUV2('ave',Hf.aave);
[Uc.amax,Vc.amax]=reduceUV2('max',Hf.amax);
[Uc.emin,Vc.emin]=reduceUV('min',Uf.emin,Vf.emin);
[Uc.eave,Vc.eave]=reduceUV('ave',Uf.eave,Vf.eave);
[Uc.emax,Vc.emax]=reduceUV('max',Uf.emax,Vf.emax);
[Uc.eaave,Vc.eaave]=reduceUV2('ave',Hf.aave);
%[Uc.eaave,Vc.eaave]=reduceUV2('ave',Hf.eave);

checksill('Coarse grid reduce',Hf,Uf,Vf,FG)
if vizdebug ;spn=spn+1;
 myfig(spn); myplot(CG.X,CG.Y,Hc.emin,Uc.emin,Vc.emin,Hrng);plotcolorbar
 title(sprintf('Coarse grid after simple reduction: Level %i minimum',lvl));drawnow
end

checkCoarseFineConnections(OriginalDeepestConnections,Uc.emin,Vc.emin,'reduceUV')

if Control.doBoundH
 [Hc.emin,msk]=boundHbyUV(Hc.emin,Uc.emin,Vc.emin);
 if vizdebug | viznorm ;spn=spn+1;
  myfig(spn); myplot(CG.X,CG.Y,Hc.emin,Uc.emin,Vc.emin,Hrng);plotcolorbar
  title(sprintf('Coarse grid after bounding H by Hu,Hv: Level %i minimum',lvl));drawnow
  showChangedCells(CG.X,CG.Y,msk);
  overlaygrid(CG.X(1:2:end,1:2:end),CG.Y(1:2:end,1:2:end))
 end
 % Avereage must be taller than emin
 Uc.eave=max(Uc.eave,Uc.emin); Vc.eave=max(Vc.eave,Vc.emin);
 [Hc.eave,msk]=boundHbyUV(Hc.eave,Uc.eave,Vc.eave); % This step biases the average low
 % Max must be taller than average
 Uc.emax=max(Uc.emax,Uc.eave); Vc.emax=max(Vc.emax,Vc.eave);
 if vizdebug | viznorm ;spn=spn+1;
  myfig(spn); myplot(CG.X,CG.Y,Hc.eave,Uc.eave,Vc.eave,Hrng);plotcolorbar
  title(sprintf('Coarse grid after bounding H by Hu,Hv: Level %i average',lvl));drawnow
  showChangedCells(CG.X,CG.Y,msk);
  overlaygrid(CG.X(1:2:end,1:2:end),CG.Y(1:2:end,1:2:end))
 end
end

% Make CUval,CVval consistent with coarse CHval (minimum)
if Control.doGenHxHy
 locdisp('  Re-generating Hu,Hv at end of algorithm')
 [Uc.emin,Vc.emin,msk]=gen_hxhy(Hc.emin,Uc.emin,Vc.emin);
 checksill('CG gen_hxhy',Hf,Uf,Vf,FG)
 if vizdebug | viznorm ;spn=spn+1;
  myfig(spn); myplot(CG.X,CG.Y,Hc.emin,Uc.emin,Vc.emin,Hrng);plotcolorbar
  title(sprintf('Coarse grid after regeneration of Hu,Hv: Level %i minimum',lvl));drawnow
  showChangedCells(CG.X,CG.Y,msk);
  overlaygrid(CG.X(1:2:end,1:2:end),CG.Y(1:2:end,1:2:end))
 end
 warn(msk,'gen_hxhy')
 % Make CUval,CVval consistent with coarse CHval (effective)
 Hc.eave=max(Hc.eave,Hc.emin);
 Uc.eave=max(Uc.eave,Uc.emin); Vc.eave=max(Vc.eave,Vc.emin);
 [Uc.eave,Vc.eave,msk]=gen_hxhy(Hc.eave,Uc.eave,Vc.eave);
 Uc.eaave=max(Uc.eaave,Uc.emin); Vc.eaave=max(Vc.eaave,Vc.emin);
 if vizdebug | viznorm ;spn=spn+1;
  myfig(spn); myplot(CG.X,CG.Y,Hc.eave,Uc.eave,Vc.eave,Hrng);plotcolorbar
  title(sprintf('Coarse grid after regeneration of Hu,Hv: Level %i average',lvl));drawnow
  showChangedCells(CG.X,CG.Y,msk);
  overlaygrid(CG.X(1:2:end,1:2:end),CG.Y(1:2:end,1:2:end))
 end
end

% Fill pot holes
if Control.doPotholes
 [Hc.emin,msk]=fillPotHoles(Hc.emin,Uc.emin,Vc.emin);
 if vizdebug ;spn=spn+1;
  myfig(spn); myplot(CG.X,CG.Y,Hc.emin,Uc.emin,Vc.emin,Hrng);plotcolorbar
  showChangedCells(CG.X,CG.Y,msk);
  title(sprintf('Coarse grid after filling pot holes: Level %i minimum',lvl));drawnow
 end
 [Hc.eave,msk]=fillPotHoles(Hc.eave,Uc.eave,Vc.eave);
 if vizdebug ;spn=spn+1;
  myfig(spn); myplot(CG.X,CG.Y,Hc.eave,Uc.eave,Vc.eave,Hrng);plotcolorbar
  showChangedCells(CG.X,CG.Y,msk);
  title(sprintf('Coarse grid after filling pot holes: Level %i average',lvl));drawnow
 end
end

% Check that CUval and CVval don't violate basic rules
check_hxhy(Hc.emin,Uc.emin,Vc.emin)
check_hxhy(Hc.eave,Uc.eave,Vc.eave)

if (vizdebug | viznorm) & Control.showAverage ;spn=spn+1;
 myfig(spn); myplot(CG.X,CG.Y,Hc.aave,Uc.lave,Vc.lave,Hrng);plotcolorbar
 title(sprintf('Coarse grid after regeneration of Hu,Hv: Level %i average',lvl));drawnow
 overlaygrid(CG.X(1:2:end,1:2:end),CG.Y(1:2:end,1:2:end))
end

% Make CUval,CVval consistent with coarse CHval (effective)
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
function [Hout] = penmin(Hvec)
% Returns the penultimate minimum value
% Equivalent to second value of sort(Hvec)
 minH=min(Hvec);
 for j=1:length(Hvec)
  if Hvec(j)==minH
   Hout=min(Hvec([1:j-1 j+1:end]));
   return
  end
 end
% -----------------------------------------------------------------------------
function [Hout] = penminV(varargin)
% Returns the penultimate minimum value of multiple arguments
% Equivalent to second value of sort(Hvec)
Nxy=size(varargin{1});
n=length(varargin);
HV=zeros(n,prod(Nxy));
for k=1:n
 HV(k,:)=reshape(varargin{k},[1 prod(Nxy)]);
end
HV=sort(HV);
Hout=reshape(HV(2,:),Nxy);
% -----------------------------------------------------------------------------
function [Hout] = reduceH(cellMeth,Hin)
[nxf,nyf]=size(Hin);
iw=1:2:nxf; ie=2:2:nxf; js=1:2:nyf; jn=2:2:nyf; % Cell indices
switch cellMeth
 case 'ave'
  Hout=((Hin(iw,js)+Hin(ie,jn))+(Hin(iw,jn)+Hin(ie,js)))/4;
 case 'min'
  Hout=mmin(Hin(iw,js),Hin(ie,jn),Hin(iw,jn),Hin(ie,js));
 case 'max'
  Hout=mmax(Hin(iw,js),Hin(ie,jn),Hin(iw,jn),Hin(ie,js));
 otherwise
  locdisp('!!!!!!!!!!!!!!     -none used-')
end
% -----------------------------------------------------------------------------
function [Uout,Vout] = reduceUV(edgeMeth,Uin,Vin)
nxf=size(Vin,1); nyf=size(Uin,2);
iw=1:2:nxf; ie=2:2:nxf; js=1:2:nyf; jn=2:2:nyf; % Cell indices
Iu=1:2:nxf+1;Jv=1:2:nyf+1; % Edge indices
switch edgeMeth
 case 'ave'
  Uout=(Uin(Iu,js)+Uin(Iu,jn))/2; Vout=(Vin(iw,Jv)+Vin(ie,Jv))/2;
 case 'min'
  Uout=min(Uin(Iu,js),Uin(Iu,jn)); Vout=min(Vin(iw,Jv),Vin(ie,Jv));
 case 'max'
  Uout=max(Uin(Iu,js),Uin(Iu,jn)); Vout=max(Vin(iw,Jv),Vin(ie,Jv));
 case 'fave'
  alph=1/2;
  Il=max(1,Iu-1); Ir=min(nxf+1,Iu+1);
  Uout=(Uin(Iu,js)+Uin(Iu,jn))/2*(1-alph) ...
      +(Uin(Il,js)+Uin(Il,jn) ...
      + Uin(Iu,js)+Uin(Ir,jn))/4*alph;
  Jl=max(1,Jv-1); Jr=min(nyf+1,Jv+1);
  Vout=(Vin(iw,Jv)+Vin(ie,Jv))/2*(1-alph) ...
      +(Vin(iw,Jl)+Vin(ie,Jl) ...
      + Vin(iw,Jr)+Vin(ie,Jr))/4*alph;
 otherwise
  locdisp('!!!!!!!!!!!!!!     -none used-')
end
% -----------------------------------------------------------------------------
function [Uout,Vout] = reduceUV2(redMeth,Hin)
% Reduce four fine-grid cell center data to single coarse edge values
nxf=size(Hin,1); nyf=size(Hin,2);
iw=1:2:nxf; ie=2:2:nxf; js=1:2:nyf; jn=2:2:nyf; % Cell indices
ir=min(nxf,1:2:nxf+1); il=max(1,0:2:nxf);
ju=min(nyf,1:2:nyf+1); jd=max(1,0:2:nyf);
switch redMeth
 case 'ave'
  Uout=(Hin(il,jn)+Hin(ir,jn)+Hin(il,js)+Hin(ir,js))/4;
  Vout=(Hin(ie,ju)+Hin(ie,jd)+Hin(iw,ju)+Hin(iw,jd))/4;
 case 'min'
  Uout=min(min(Hin(il,jn),Hin(ir,jn)),min(Hin(il,js),Hin(ir,js)));
  Vout=min(min(Hin(ie,ju),Hin(ie,jd)),min(Hin(iw,ju),Hin(iw,jd)));
 case 'max'
  Uout=max(max(Hin(il,jn),Hin(ir,jn)),max(Hin(il,js),Hin(ir,js)));
  Vout=max(max(Hin(ie,ju),Hin(ie,jd)),max(Hin(iw,ju),Hin(iw,jd)));
 otherwise
  locdisp('!!!!!!!!!!!!!!     -none used-')
end
% -----------------------------------------------------------------------------
function [Hout,Uout,Vout,msk] = pushCorners(Hin,Uin,Vin)
 % A convex corner within a coarse grid cell can be made into a
 % concave corner without changing connectivity across the major
 % parts of the cell. The cross-corner connection for the minor
 % part of the cell may be eliminated.
 [nxf,nyf]=size(Hin.emin); msk=NaN*Hin.emin;
 ic=1:nxf/2;jc=1:nyf/2;
 Hout=Hin;Uout=Uin;Vout=Vin;
 nncnt=0;
 global lowMem
 if lowMem
  for J=jc
   j=2*J-1;
   for I=ic 
    i=2*I-1;
    ncnt=0;
    % SW box
    Hcrnr=min(Uin.emin(i+1,j),Vin.emin(i,j+1)); Hoppo=max(Uin.emin(i+1,j+1),Vin.emin(i+1,j+1));
    if Hcrnr>Hoppo
     Uout.emin(i+1,j)=Hoppo; Vout.emin(i,j+1)=Hoppo;
     Hout.emin(i,j)=Hoppo;
     Uout.emin(i,j)=max(Uout.emin(i,j),Hcrnr);
     Vout.emin(i,j)=max(Vout.emin(i,j),Hcrnr);
     HcrnrAve=(Uin.eave(i+1,j)+Vin.eave(i,j+1))/2;
     Uout.eave(i,j)=max(Uout.eave(i,j),HcrnrAve);
     Vout.eave(i,j)=max(Vout.eave(i,j),HcrnrAve);
     Hout.eave(i,j)=(Hin.eave(i+1,j)+Hin.eave(i,j+1)+Hin.eave(i+1,j+1))/3;
     HcrnrMax=max(Uin.emax(i+1,j),Vin.emax(i,j+1));
     Uout.emax(i,j)=max(Uout.emax(i,j),HcrnrMax);
     Vout.emax(i,j)=max(Vout.emax(i,j),HcrnrMax);
     Hout.emax(i,j)=Hoppo;
     msk(i,j)=1; % debugging
     ncnt=ncnt+1;
    end
    % SE box
    Hcrnr=min(Uin.emin(i+1,j),Vin.emin(i+1,j+1)); Hoppo=max(Uin.emin(i+1,j+1),Vin.emin(i,j+1));
    if Hcrnr>Hoppo
     Uout.emin(i+1,j)=Hoppo; Vout.emin(i+1,j+1)=Hoppo;
     Hout.emin(i+1,j)=Hoppo;
     Uout.emin(i+2,j)=max(Uout.emin(i+2,j),Hcrnr);
     Vout.emin(i+1,j)=max(Vout.emin(i+1,j),Hcrnr);
     HcrnrAve=(Uin.eave(i+1,j)+Vin.eave(i+1,j+1))/2;
     Uout.eave(i+2,j)=max(Uout.eave(i+2,j),HcrnrAve);
     Vout.eave(i+1,j)=max(Vout.eave(i+1,j),HcrnrAve);
     Hout.eave(i+1,j)=(Hin.eave(i,j)+Hin.eave(i,j+1)+Hin.eave(i+1,j+1))/3;
     HcrnrMax=max(Uin.emax(i+1,j),Vin.emax(i+1,j+1));
     Uout.emax(i+2,j)=max(Uout.emax(i+2,j),HcrnrMax);
     Vout.emax(i+1,j)=max(Vout.emax(i+1,j),HcrnrMax);
     Hout.emax(i+1,j)=Hoppo;
     msk(i+1,j)=1; % debugging
     ncnt=ncnt+1;
    end
    % NW box
    Hcrnr=min(Uin.emin(i+1,j+1),Vin.emin(i,j+1)); Hoppo=max(Uin.emin(i+1,j),Vin.emin(i+1,j+1));
    if Hcrnr>Hoppo
     Uout.emin(i+1,j+1)=Hoppo; Vout.emin(i,j+1)=Hoppo;
     Hout.emin(i,j+1)=Hoppo;
     Uout.emin(i,j+1)=max(Uout.emin(i,j+1),Hcrnr);
     Vout.emin(i,j+2)=max(Vout.emin(i,j+2),Hcrnr);
     HcrnrAve=(Uin.eave(i+1,j+1)+Vin.eave(i,j+1))/2;
     Uout.eave(i,j+1)=max(Uout.eave(i,j+1),HcrnrAve);
     Vout.eave(i,j+2)=max(Vout.eave(i,j+2),HcrnrAve);
     Hout.eave(i,j+1)=(Hin.eave(i,j)+Hin.eave(i+1,j)+Hin.eave(i+1,j+1))/3;
     HcrnrMax=max(Uin.emax(i+1,j+1),Vin.emax(i,j+1));
     Uout.emax(i,j+1)=max(Uout.emax(i,j+1),HcrnrMax);
     Vout.emax(i,j+2)=max(Vout.emax(i,j+2),HcrnrMax);
     Hout.emax(i,j+1)=Hoppo;
     msk(i,j+1)=1; % debugging
     ncnt=ncnt+1;
    end
    % NE box
    Hcrnr=min(Uin.emin(i+1,j+1),Vin.emin(i+1,j+1)); Hoppo=max(Uin.emin(i+1,j),Vin.emin(i,j+1));
    if Hcrnr>Hoppo
     Uout.emin(i+1,j+1)=Hoppo; Vout.emin(i+1,j+1)=Hoppo;
     Hout.emin(i+1,j+1)=Hoppo;
     Uout.emin(i+2,j+1)=max(Uout.emin(i+2,j+1),Hcrnr);
     Vout.emin(i+1,j+2)=max(Vout.emin(i+1,j+2),Hcrnr);
     HcrnrAve=(Uin.eave(i+1,j+1)+Vin.eave(i+1,j+1))/2;
     Uout.eave(i+2,j+1)=max(Uout.eave(i+2,j+1),HcrnrAve);
     Vout.eave(i+1,j+2)=max(Vout.eave(i+1,j+2),HcrnrAve);
     Hout.eave(i+1,j+1)=(Hin.eave(i,j)+Hin.eave(i+1,j)+Hin.eave(i,j+1))/3;
     HcrnrMax=max(Uin.emax(i+1,j+1),Vin.emax(i+1,j+1));
     Uout.emax(i+2,j+1)=max(Uout.emax(i+2,j+1),HcrnrMax);
     Vout.emax(i+1,j+2)=max(Vout.emax(i+1,j+2),HcrnrMax);
     Hout.emax(i+1,j+1)=Hoppo;
     msk(i+1,j+1)=1; % debugging
     ncnt=ncnt+1;
    end
    if ncnt>1
      error('This should never happen (while pushing out corners)!')
    end
    nncnt=nncnt+ncnt;
   end % I
  end % J
 else % !lowMem
 moveAveMax=1; % 1 turns on moving, NaN turns off
  Nh=size(Hout.emin);Nu=size(Uout.emin);Nv=size(Vout.emin);
  % SW box
  i=2*ic-1; j=2*jc-1;
  Hcrnr=min(Uin.emin(i+1,j),Vin.emin(i,j+1)); Hoppo=max(Uin.emin(i+1,j+1),Vin.emin(i+1,j+1));
  HcrnrAve=(Uin.eave(i+1,j)+Vin.eave(i,j+1))/2;
  HcrnrMax=max(Uin.emax(i+1,j),Vin.emax(i,j+1));
  HcrnrOave=(Hin.eave(i+1,j)+Hin.eave(i,j+1)+Hin.eave(i+1,j+1))/3;
  [i,j]=find( Hcrnr>Hoppo ); I=2*i-1;J=2*j-1;
  Kh=sub2ind(Nh,I,J);k=sub2ind(size(Hcrnr),i,j);
  Ku=sub2ind(Nu,I,J); Kv=sub2ind(Nv,I,J);
  Uout.emin(sub2ind(Nu,I+1,J))=Hoppo(k);
  Vout.emin(sub2ind(Nv,I,J+1))=Hoppo(k);
  Hout.emin(Kh)=Hoppo(k);
  Uout.emin(Ku)=max(Uout.emin(Ku),Hcrnr(k));
  Vout.emin(Kv)=max(Vout.emin(Kv),Hcrnr(k));
  Uout.eave(Ku)=max(moveAveMax*Uout.eave(Ku),HcrnrAve(k));
  Vout.eave(Kv)=max(moveAveMax*Vout.eave(Kv),HcrnrAve(k));
  Hout.eave(Kh)=HcrnrOave(k);
  Uout.emax(Ku)=max(moveAveMax*Uout.emax(Ku),HcrnrMax(k));
  Vout.emax(Kv)=max(moveAveMax*Vout.emax(Kv),HcrnrMax(k));
  Hout.emax(Kh)=Hoppo(k);
  msk(Kh)=1; % debugging
  % SE box
  i=2*ic-1; j=2*jc-1;
  Hcrnr=min(Uin.emin(i+1,j),Vin.emin(i+1,j+1)); Hoppo=max(Uin.emin(i+1,j+1),Vin.emin(i,j+1));
  HcrnrAve=(Uin.eave(i+1,j)+Vin.eave(i+1,j+1))/2;
  HcrnrMax=max(Uin.emax(i+1,j),Vin.emax(i+1,j+1));
  HcrnrOave=(Hin.eave(i,j)+Hin.eave(i,j+1)+Hin.eave(i+1,j+1))/3;
  [i,j]=find( Hcrnr>Hoppo ); I=2*i-1;J=2*j-1;
  Kh=sub2ind(Nh,I+1,J);k=sub2ind(size(Hcrnr),i,j);
  Ku=sub2ind(Nu,I+2,J); Kv=sub2ind(Nv,I+1,J);
  Uout.emin(sub2ind(Nu,I+1,J))=Hoppo(k);
  Vout.emin(sub2ind(Nv,I+1,J+1))=Hoppo(k);
  Hout.emin(Kh)=Hoppo(k);
  Uout.emin(Ku)=max(Uout.emin(Ku),Hcrnr(k));
  Vout.emin(Kv)=max(Vout.emin(Kv),Hcrnr(k));
  Uout.eave(Ku)=max(moveAveMax*Uout.eave(Ku),HcrnrAve(k));
  Vout.eave(Kv)=max(moveAveMax*Vout.eave(Kv),HcrnrAve(k));
  Hout.eave(Kh)=HcrnrOave(k);
  Uout.emax(Ku)=max(moveAveMax*Uout.emax(Ku),HcrnrMax(k));
  Vout.emax(Kv)=max(moveAveMax*Vout.emax(Kv),HcrnrMax(k));
  Hout.emax(Kh)=Hoppo(k);
  msk(Kh)=1; % debugging
  % NW box
  i=2*ic-1; j=2*jc-1;
  Hcrnr=min(Uin.emin(i+1,j+1),Vin.emin(i,j+1)); Hoppo=max(Uin.emin(i+1,j),Vin.emin(i+1,j+1));
  HcrnrAve=(Uin.eave(i+1,j+1)+Vin.eave(i,j+1))/2;
  HcrnrMax=max(Uin.emax(i+1,j+1),Vin.emax(i,j+1));
  HcrnrOave=(Hin.eave(i,j)+Hin.eave(i+1,j)+Hin.eave(i+1,j+1))/3;
  [i,j]=find( Hcrnr>Hoppo ); I=2*i-1;J=2*j-1;
  Kh=sub2ind(Nh,I,J+1);k=sub2ind(size(Hcrnr),i,j);
  Ku=sub2ind(Nu,I,J+1); Kv=sub2ind(Nv,I,J+2);
  Uout.emin(sub2ind(Nu,I+1,J+1))=Hoppo(k);
  Vout.emin(sub2ind(Nv,I,J+1))=Hoppo(k);
  Uout.emin(Ku)=max(Uout.emin(Ku),Hcrnr(k));
  Vout.emin(Kv)=max(Vout.emin(Kv),Hcrnr(k));
  Hout.emin(Kh)=Hoppo(k);
  Uout.eave(Ku)=max(moveAveMax*Uout.eave(Ku),HcrnrAve(k));
  Vout.eave(Kv)=max(moveAveMax*Vout.eave(Kv),HcrnrAve(k));
  Hout.eave(Kh)=HcrnrOave(k);
  Uout.emax(Ku)=max(moveAveMax*Uout.emax(Ku),HcrnrMax(k));
  Vout.emax(Kv)=max(moveAveMax*Vout.emax(Kv),HcrnrMax(k));
  Hout.emax(Kh)=Hoppo(k);
  msk(Kh)=1; % debugging
  % NE box
  i=2*ic-1; j=2*jc-1;
  Hcrnr=min(Uin.emin(i+1,j+1),Vin.emin(i+1,j+1)); Hoppo=max(Uin.emin(i+1,j),Vin.emin(i,j+1));
  HcrnrAve=(Uin.eave(i+1,j+1)+Vin.eave(i+1,j+1))/2;
  HcrnrMax=max(Uin.emax(i+1,j+1),Vin.emax(i+1,j+1));
  HcrnrOave=(Hin.eave(i,j)+Hin.eave(i+1,j)+Hin.eave(i,j+1))/3;
  [i,j]=find( Hcrnr>Hoppo ); I=2*i-1;J=2*j-1;
  Kh=sub2ind(Nh,I+1,J+1);k=sub2ind(size(Hcrnr),i,j);
  Ku=sub2ind(Nu,I+2,J+1); Kv=sub2ind(Nv,I+1,J+2);
  Uout.emin(sub2ind(Nu,I+1,J+1))=Hoppo(k);
  Vout.emin(sub2ind(Nv,I+1,J+1))=Hoppo(k);
  Uout.emin(Ku)=max(Uout.emin(Ku),Hcrnr(k));
  Vout.emin(Kv)=max(Vout.emin(Kv),Hcrnr(k));
  Hout.emin(Kh)=Hoppo(k);
  Uout.eave(Ku)=max(moveAveMax*Uout.eave(Ku),HcrnrAve(k));
  Vout.eave(Kv)=max(moveAveMax*Vout.eave(Kv),HcrnrAve(k));
  Hout.eave(Kh)=HcrnrOave(k);
  Uout.emax(Ku)=max(moveAveMax*Uout.emax(Ku),HcrnrMax(k));
  Vout.emax(Kv)=max(moveAveMax*Vout.emax(Kv),HcrnrMax(k));
  Hout.emax(Kh)=Hoppo(k);
  msk(Kh)=1; % debugging
 end
 locdisp(sprintf('  cornerPush: %i points (%5.2f%%) "pushed"', ...
   nncnt,100*nncnt/(nxf*nyf) ))
% -----------------------------------------------------------------------------
function [Hout,Uout,Vout,msk] = centralRidge(Hin,Uin,Vin,msg)
% Snaps central ridges into concave edge ridges
 [nxf,nyf]=size(Hin.emin); msk=NaN*Hin.emin;
 ic=1:nxf/2;jc=1:nyf/2;
 Hout=Hin;Uout=Uin;Vout=Vin;
 nncnt=0;
 global lowMem
 moveAveMax=1; % 1 turns on moving, NaN turns off
 if lowMem
  for J=jc
   j=2*J-1;
   for I=ic 
    i=2*I-1;
    ncnt=0;
    % Ridge running east-west, blocking north-south
    Hridge=min(Vin.emin(i,j+1),Vin.emin(i+1,j+1));
    HridgeAve=(Vin.eave(i,j+1)+Vin.eave(i+1,j+1))/2;
    HridgeMax=max(Vin.emax(i,j+1),Vin.emax(i+1,j+1));
    HoppoMax=max(Uin.emin(i+1,j),Uin.emin(i+1,j+1));
    HoppoMin=min(Uin.emin(i+1,j),Uin.emin(i+1,j+1));
    if Hridge>=HoppoMax && Hridge>HoppoMin % This is an east-west ridge
     if Uin.emin(i+1,j+1)>Uin.emin(i+1,j) || ...
        ( Uin.emin(i+1,j)==Uin.emin(i+1,j+1) && ( ...
           Hin.emin(i,j+1)+Hin.emin(i+1,j+1)>Hin.emin(i,j)+Hin.emin(i+1,j) ...
        || Vin.emin(i,j+2)+Vin.emin(i+1,j+2)>Vin.emin(i,j)+Vin.emin(i+1,j) ) )
      Uout.emin(i,j+1)=max(Uout.emin(i,j+1),Hridge);
      Vout.emin(i,j+2)=max(Vout.emin(i,j+2),Hridge);
      Vout.emin(i+1,j+2)=max(Vout.emin(i+1,j+2),Hridge);
      Uout.emin(i+2,j+1)=max(Uout.emin(i+2,j+1),Hridge);
      Uout.eave(i,j+1)=max(moveAveMax*Uout.eave(i,j+1),HridgeAve);
      Vout.eave(i,j+2)=max(moveAveMax*Vout.eave(i,j+2),HridgeAve);
      Vout.eave(i+1,j+2)=max(moveAveMax*Vout.eave(i+1,j+2),HridgeAve);
      Uout.eave(i+2,j+1)=max(moveAveMax*Uout.eave(i+2,j+1),HridgeAve);
      Uout.emax(i,j+1)=max(moveAveMax*Uout.emax(i,j+1),HridgeMax);
      Vout.emax(i,j+2)=max(moveAveMax*Vout.emax(i,j+2),HridgeMax);
      Vout.emax(i+1,j+2)=max(moveAveMax*Vout.emax(i+1,j+2),HridgeMax);
      Uout.emax(i+2,j+1)=max(moveAveMax*Uout.emax(i+2,j+1),HridgeMax);
      Hout.emin(i:i+1,j+1)=HoppoMin;
      Vout.emin(i:i+1,j+1)=HoppoMin; Uout.emin(i+1,j+1)=HoppoMin;
      Hout.eave(i:i+1,j+1)=mean(Hin.eave(i:i+1,j));
      Hout.emax(i:i+1,j+1)=HoppoMin;
      msk(i:i+1,j+1)=1;
     elseif Uin.emin(i+1,j)>Uin.emin(i+1,j+1) || ...
        ( Uin.emin(i+1,j)==Uin.emin(i+1,j+1) && ( ...
           Hin.emin(i,j+1)+Hin.emin(i+1,j+1)<Hin.emin(i,j)+Hin.emin(i+1,j) ...
        || Vin.emin(i,j+2)+Vin.emin(i+1,j+2)<Vin.emin(i,j)+Vin.emin(i+1,j) ) )
      Uout.emin(i,j)=max(Uout.emin(i,j),Hridge);
      Vout.emin(i,j)=max(Vout.emin(i,j),Hridge);
      Vout.emin(i+1,j)=max(Vout.emin(i+1,j),Hridge);
      Uout.emin(i+2,j)=max(Uout.emin(i+2,j),Hridge);
      Uout.eave(i,j)=max(moveAveMax*Uout.eave(i,j),HridgeAve);
      Vout.eave(i,j)=max(moveAveMax*Vout.eave(i,j),HridgeAve);
      Vout.eave(i+1,j)=max(moveAveMax*Vout.eave(i+1,j),HridgeAve);
      Uout.eave(i+2,j)=max(moveAveMax*Uout.eave(i+2,j),HridgeAve);
      Uout.emax(i,j)=max(moveAveMax*Uout.emax(i,j),HridgeMax);
      Vout.emax(i,j)=max(moveAveMax*Vout.emax(i,j),HridgeMax);
      Vout.emax(i+1,j)=max(moveAveMax*Vout.emax(i+1,j),HridgeMax);
      Uout.emax(i+2,j)=max(moveAveMax*Uout.emax(i+2,j),HridgeMax);
      Hout.emin(i:i+1,j)=HoppoMin;
      Vout.emin(i:i+1,j+1)=HoppoMin; Uout.emin(i+1,j)=HoppoMin;
   %%%Hout.eave(i:i+1,j)=mean(Hout.eave(i:i+1,j+1));  %%%%BUG%%%%
      Hout.eave(i:i+1,j)=mean(Hin.eave(i:i+1,j+1));
      Hout.emax(i:i+1,j)=HoppoMin;
      msk(i:i+1,j)=1;
     elseif Uin.emin(i+1,j)==Uin.emin(i+1,j+1)
      Uout.emin(i:i+2,j:j+1)=max(Uout.emin(i:i+2,j:j+1),Hridge);
      Vout.emin(i:i+1,j:j+2)=max(Vout.emin(i:i+1,j:j+2),Hridge);
      Uout.emin(i+1,j:j+1)=Hridge;
      Vout.emin(i:i+1,j+1)=Hridge;
      Hout.emin(i:i+1,j:j+1)=Hridge*NaN;
      msk(i:i+1,j:j+1)=1;
      disp('WARNING: east-west ridge with equal side ridges!')
     else
      error('Should never ever get here (e-w)!')
     end
     ncnt=ncnt+1;
    end % east-west ridge
    % Ridge running north-south, blocking east-west
    Hridge=min(Uin.emin(i+1,j),Uin.emin(i+1,j+1));
    HridgeAve=(Uin.eave(i+1,j)+Uin.eave(i+1,j+1))/2;
    HridgeMax=max(Uin.emax(i+1,j),Uin.emax(i+1,j+1));
    HoppoMax=max(Vin.emin(i,j+1),Vin.emin(i+1,j+1));
    HoppoMin=min(Vin.emin(i,j+1),Vin.emin(i+1,j+1));
    if Hridge>=HoppoMax && Hridge>HoppoMin % This is a north-south ridge
     if Vin.emin(i+1,j+1)>Vin.emin(i,j+1) || ...
        ( Vin.emin(i,j+1)==Vin.emin(i+1,j+1) && ( ...
           Hin.emin(i+1,j)+Hin.emin(i+1,j+1)>Hin.emin(i,j)+Hin.emin(i,j+1) ...
        || Uin.emin(i+2,j)+Uin.emin(i+2,j+1)>Uin.emin(i,j)+Uin.emin(i,j+1) ) )
      Vout.emin(i+1,j)=max(Vout.emin(i+1,j),Hridge);
      Uout.emin(i+2,j)=max(Uout.emin(i+2,j),Hridge);
      Uout.emin(i+2,j+1)=max(Uout.emin(i+2,j+1),Hridge);
      Vout.emin(i+1,j+2)=max(Vout.emin(i+1,j+2),Hridge);
      Vout.eave(i+1,j)=max(moveAveMax*Vout.eave(i+1,j),HridgeAve);
      Uout.eave(i+2,j)=max(moveAveMax*Uout.eave(i+2,j),HridgeAve);
      Uout.eave(i+2,j+1)=max(moveAveMax*Uout.eave(i+2,j+1),HridgeAve);
      Vout.eave(i+1,j+2)=max(moveAveMax*Vout.eave(i+1,j+2),HridgeAve);
      Vout.emax(i+1,j)=max(moveAveMax*Vout.emax(i+1,j),HridgeMax);
      Uout.emax(i+2,j)=max(moveAveMax*Uout.emax(i+2,j),HridgeMax);
      Uout.emax(i+2,j+1)=max(moveAveMax*Uout.emax(i+2,j+1),HridgeMax);
      Vout.emax(i+1,j+2)=max(moveAveMax*Vout.emax(i+1,j+2),HridgeMax);
      Hout.emin(i+1,j:j+1)=HoppoMin;
      Uout.emin(i+1,j:j+1)=HoppoMin; Vout.emin(i+1,j+1)=HoppoMin;
      Hout.eave(i+1,j:j+1)=mean(Hin.eave(i,j:j+1));
      Hout.emax(i+1,j:j+1)=HoppoMin;  %%%%BUG should be max?
      msk(i+1,j:j+1)=1;
     elseif Vin.emin(i,j+1)>Vin.emin(i+1,j+1) || ...
        ( Vin.emin(i,j+1)==Vin.emin(i+1,j+1) && ( ...
           Hin.emin(i+1,j)+Hin.emin(i+1,j+1)<Hin.emin(i,j)+Hin.emin(i,j+1) ...
        || Uin.emin(i+2,j)+Uin.emin(i+2,j+1)<Uin.emin(i,j)+Uin.emin(i,j+1) ) )
      Vout.emin(i,j)=max(Vout.emin(i,j),Hridge);
      Uout.emin(i,j)=max(Uout.emin(i,j),Hridge);
      Uout.emin(i,j+1)=max(Uout.emin(i,j+1),Hridge);
      Vout.emin(i,j+2)=max(Vout.emin(i,j+2),Hridge);
      Vout.eave(i,j)=max(moveAveMax*Vout.eave(i,j),HridgeAve);
      Uout.eave(i,j)=max(moveAveMax*Uout.eave(i,j),HridgeAve);
      Uout.eave(i,j+1)=max(moveAveMax*Uout.eave(i,j+1),HridgeAve);
      Vout.eave(i,j+2)=max(moveAveMax*Vout.eave(i,j+2),HridgeAve);
      Vout.emax(i,j)=max(moveAveMax*Vout.emax(i,j),HridgeMax);
      Uout.emax(i,j)=max(moveAveMax*Uout.emax(i,j),HridgeMax);
      Uout.emax(i,j+1)=max(moveAveMax*Uout.emax(i,j+1),HridgeMax);
      Vout.emax(i,j+2)=max(moveAveMax*Vout.emax(i,j+2),HridgeMax);
      Hout.emin(i,j:j+1)=HoppoMin;
      Uout.emin(i+1,j:j+1)=HoppoMin; Vout.emin(i,j+1)=HoppoMin;
      Hout.eave(i,j:j+1)=mean(Hin.eave(i+1,j:j+1));
      Hout.emax(i,j:j+1)=HoppoMin;
      msk(i,j:j+1)=1;
     elseif Vin.emin(i,j+1)==Vin.emin(i+1,j+1)
      Vout.emin(i:i+1,j:j+2)=max(Vout.emin(i:i+1,j:j+2),Hridge);
      Uout.emin(i:i+2,j:j+1)=max(Uout.emin(i:i+2,j:j+1),Hridge);
      Vout.emin(i:i+1,j+1)=Hridge;
      Uout.emin(i+1,j:j+1)=Hridge;
      Hout.emin(i:i+1,j:j+1)=Hridge*NaN;
      msk(i:i+1,j:j+1)=1;
      disp('WARNING: north-south ridge with equal side ridges!')
     else
      error('Should never ever get here (n-s)!')
     end
     ncnt=ncnt+1;
    end % north-south ridge
    if ncnt>1
      error('This should never happen (while handling central ridges)!')
    end
    nncnt=nncnt+ncnt;
   end % I
  end % J
 else % !lowMem
  Nu=size(Uout.emin);Nv=size(Vout.emin);Nh=size(Hout.emin);
  i=2*ic-1; j=2*jc-1;
    % Ridge running east-west, blocking north-south
    Hridge=min(Vin.emin(i,j+1),Vin.emin(i+1,j+1));
    HridgeAve=(Vin.eave(i,j+1)+Vin.eave(i+1,j+1))/2;
    HridgeMax=max(Vin.emax(i,j+1),Vin.emax(i+1,j+1));
    HoppoMax=max(Uin.emin(i+1,j),Uin.emin(i+1,j+1));
    HoppoMin=min(Uin.emin(i+1,j),Uin.emin(i+1,j+1));
  mskEW=Hridge>=HoppoMax & Hridge>HoppoMin;
  mskC=0*mskEW;
  [i,j]=find( mskEW & ( Uin.emin(i+1,j+1)>Uin.emin(i+1,j) ...
            | ( Uin.emin(i+1,j)==Uin.emin(i+1,j+1) ...
              & ( Hin.emin(i,j+1)+Hin.emin(i+1,j+1)>Hin.emin(i,j)+Hin.emin(i+1,j) ...
                | Vin.emin(i,j+2)+Vin.emin(i+1,j+2)>Vin.emin(i,j)+Vin.emin(i+1,j) ...
            ) ) ) );
  I=2*i-1;J=2*j-1; k=sub2ind(size(Hridge),i,j);
      Uout.emin(sub2ind(Nu,I,J+1))=max(Uout.emin(sub2ind(Nu,I,J+1)),Hridge(k));
      Vout.emin(sub2ind(Nv,I,J+2))=max(Vout.emin(sub2ind(Nv,I,J+2)),Hridge(k));
      Vout.emin(sub2ind(Nv,I+1,J+2))=max(Vout.emin(sub2ind(Nv,I+1,J+2)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+2,J+1))=max(Uout.emin(sub2ind(Nu,I+2,J+1)),Hridge(k));
      Uout.eave(sub2ind(Nu,I,J+1))=max(moveAveMax*Uout.eave(sub2ind(Nu,I,J+1)),HridgeAve(k));
      Vout.eave(sub2ind(Nv,I,J+2))=max(moveAveMax*Vout.eave(sub2ind(Nv,I,J+2)),HridgeAve(k));
      Vout.eave(sub2ind(Nv,I+1,J+2))=max(moveAveMax*Vout.eave(sub2ind(Nv,I+1,J+2)),HridgeAve(k));
      Uout.eave(sub2ind(Nu,I+2,J+1))=max(moveAveMax*Uout.eave(sub2ind(Nu,I+2,J+1)),HridgeAve(k));
      Uout.emax(sub2ind(Nu,I,J+1))=max(moveAveMax*Uout.emax(sub2ind(Nu,I,J+1)),HridgeMax(k));
      Vout.emax(sub2ind(Nv,I,J+2))=max(moveAveMax*Vout.emax(sub2ind(Nv,I,J+2)),HridgeMax(k));
      Vout.emax(sub2ind(Nv,I+1,J+2))=max(moveAveMax*Vout.emax(sub2ind(Nv,I+1,J+2)),HridgeMax(k));
      Uout.emax(sub2ind(Nu,I+2,J+1))=max(moveAveMax*Uout.emax(sub2ind(Nu,I+2,J+1)),HridgeMax(k));
      Hout.emin(sub2ind(Nh,I,J+1))=HoppoMin(k);
      Hout.emin(sub2ind(Nh,I+1,J+1))=HoppoMin(k);
      Vout.emin(sub2ind(Nv,I,J+1))=HoppoMin(k);
      Vout.emin(sub2ind(Nv,I+1,J+1))=HoppoMin(k);
      Uout.emin(sub2ind(Nu,I+1,J+1))=HoppoMin(k);
      Hout.eave(sub2ind(Nh,I,J+1))=(Hin.eave(sub2ind(Nh,I,J))+Hin.eave(sub2ind(Nh,I+1,J)))/2;
      Hout.eave(sub2ind(Nh,I+1,J+1))=(Hin.eave(sub2ind(Nh,I,J))+Hin.eave(sub2ind(Nh,I+1,J)))/2;
      Hout.emax(sub2ind(Nh,I,J+1))=HoppoMin(k);
      Hout.emax(sub2ind(Nh,I+1,J+1))=HoppoMin(k);
      msk(sub2ind(Nh,I,J+1))=1;
      msk(sub2ind(Nh,I+1,J+1))=1;
      mskC(k)=1;
  i=2*ic-1; j=2*jc-1;
  [i,j]=find( mskEW & ( Uin.emin(i+1,j)>Uin.emin(i+1,j+1) ...
            | ( Uin.emin(i+1,j)==Uin.emin(i+1,j+1) ...
              & ( Hin.emin(i,j+1)+Hin.emin(i+1,j+1)<Hin.emin(i,j)+Hin.emin(i+1,j) ...
                | Vin.emin(i,j+2)+Vin.emin(i+1,j+2)<Vin.emin(i,j)+Vin.emin(i+1,j) ...
            ) ) ) );
  I=2*i-1;J=2*j-1; k=sub2ind(size(Hridge),i,j);
      Uout.emin(sub2ind(Nu,I,J))=max(Uout.emin(sub2ind(Nu,I,J)),Hridge(k));
      Vout.emin(sub2ind(Nv,I,J))=max(Vout.emin(sub2ind(Nv,I,J)),Hridge(k));
      Vout.emin(sub2ind(Nv,I+1,J))=max(Vout.emin(sub2ind(Nv,I+1,J)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+2,J))=max(Uout.emin(sub2ind(Nu,I+2,J)),Hridge(k));
      Uout.eave(sub2ind(Nu,I,J))=max(moveAveMax*Uout.eave(sub2ind(Nu,I,J)),HridgeAve(k));
      Vout.eave(sub2ind(Nv,I,J))=max(moveAveMax*Vout.eave(sub2ind(Nv,I,J)),HridgeAve(k));
      Vout.eave(sub2ind(Nv,I+1,J))=max(moveAveMax*Vout.eave(sub2ind(Nv,I+1,J)),HridgeAve(k));
      Uout.eave(sub2ind(Nu,I+2,J))=max(moveAveMax*Uout.eave(sub2ind(Nu,I+2,J)),HridgeAve(k));
      Uout.emax(sub2ind(Nu,I,J))=max(moveAveMax*Uout.emax(sub2ind(Nu,I,J)),HridgeMax(k));
      Vout.emax(sub2ind(Nv,I,J))=max(moveAveMax*Vout.emax(sub2ind(Nv,I,J)),HridgeMax(k));
      Vout.emax(sub2ind(Nv,I+1,J))=max(moveAveMax*Vout.emax(sub2ind(Nv,I+1,J)),HridgeMax(k));
      Uout.emax(sub2ind(Nu,I+2,J))=max(moveAveMax*Uout.emax(sub2ind(Nu,I+2,J)),HridgeMax(k));
      Hout.emin(sub2ind(Nh,I,J))=HoppoMin(k);
      Hout.emin(sub2ind(Nh,I+1,J))=HoppoMin(k);
      Vout.emin(sub2ind(Nv,I,J+1))=HoppoMin(k);
      Vout.emin(sub2ind(Nv,I+1,J+1))=HoppoMin(k);
      Uout.emin(sub2ind(Nu,I+1,J))=HoppoMin(k);
      Hout.eave(sub2ind(Nh,I,J))=(Hin.eave(sub2ind(Nh,I,J+1))+Hin.eave(sub2ind(Nh,I+1,J+1)))/2;
      Hout.eave(sub2ind(Nh,I+1,J))=(Hin.eave(sub2ind(Nh,I,J+1))+Hin.eave(sub2ind(Nh,I+1,J+1)))/2;
      Hout.emax(sub2ind(Nh,I,J))=HoppoMin(k);
      Hout.emax(sub2ind(Nh,I+1,J))=HoppoMin(k);
      msk(sub2ind(Nh,I,J))=1;
      msk(sub2ind(Nh,I+1,J))=1;
      mskC(k)=1;
  i=2*ic-1; j=2*jc-1;
  [i,j]=find( mskEW & mskC==0 & ( Uin.emin(i+1,j)==Uin.emin(i+1,j+1) ) );
  I=2*i-1;J=2*j-1; k=sub2ind(size(Hridge),i,j);
      Uout.emin(sub2ind(Nu,I,J))=max(Uout.emin(sub2ind(Nu,I,J)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+2,J))=max(Uout.emin(sub2ind(Nu,I+2,J)),Hridge(k));
      Uout.emin(sub2ind(Nu,I,J+1))=max(Uout.emin(sub2ind(Nu,I,J+1)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+2,J+1))=max(Uout.emin(sub2ind(Nu,I+2,J+1)),Hridge(k));
      Vout.emin(sub2ind(Nv,I,J))=max(Vout.emin(sub2ind(Nv,I,J)),Hridge(k));
      Vout.emin(sub2ind(Nv,I,J+2))=max(Vout.emin(sub2ind(Nv,I,J+2)),Hridge(k));
      Vout.emin(sub2ind(Nv,I+1,J))=max(Vout.emin(sub2ind(Nv,I+1,J)),Hridge(k));
      Vout.emin(sub2ind(Nv,I+1,J+2))=max(Vout.emin(sub2ind(Nv,I+1,J+2)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+1,J))=Hridge(k);
      Uout.emin(sub2ind(Nu,I+1,J+1))=Hridge(k);
      Vout.emin(sub2ind(Nv,I,J+1))=Hridge(k);
      Vout.emin(sub2ind(Nv,I+1,J+1))=Hridge(k);
      Hout.emin(sub2ind(Nh,I,J))=Hridge(k)*NaN;
      Hout.emin(sub2ind(Nh,I+1,J))=Hridge(k)*NaN;
      Hout.emin(sub2ind(Nh,I,J+1))=Hridge(k)*NaN;
      Hout.emin(sub2ind(Nh,I+1,J+1))=Hridge(k)*NaN;
      msk(sub2ind(Nh,I,J))=1;
      msk(sub2ind(Nh,I+1,J))=1;
      msk(sub2ind(Nh,I,J+1))=1;
      msk(sub2ind(Nh,I+1,J+1))=1;

  i=2*ic-1; j=2*jc-1;
    % Ridge running north-south, blocking east-west
    Hridge=min(Uin.emin(i+1,j),Uin.emin(i+1,j+1));
    HridgeAve=(Uin.eave(i+1,j)+Uin.eave(i+1,j+1))/2;
    HridgeMax=max(Uin.emax(i+1,j),Uin.emax(i+1,j+1));
    HoppoMax=max(Vin.emin(i,j+1),Vin.emin(i+1,j+1));
    HoppoMin=min(Vin.emin(i,j+1),Vin.emin(i+1,j+1));
  mskEW=Hridge>=HoppoMax & Hridge>HoppoMin;
  mskC=0*mskEW;
  [i,j]=find( mskEW & ( Vin.emin(i+1,j+1)>Vin.emin(i,j+1) ...
            | ( Vin.emin(i,j+1)==Vin.emin(i+1,j+1) ...
              & ( Hin.emin(i+1,j)+Hin.emin(i+1,j+1)>Hin.emin(i,j)+Hin.emin(i,j+1) ...
                | Uin.emin(i+2,j)+Uin.emin(i+2,j+1)>Uin.emin(i,j)+Uin.emin(i,j+1) ...
            ) ) ) );
  I=2*i-1;J=2*j-1; k=sub2ind(size(Hridge),i,j);
      Vout.emin(sub2ind(Nv,I+1,J))=max(Vout.emin(sub2ind(Nv,I+1,J)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+2,J))=max(Uout.emin(sub2ind(Nu,I+2,J)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+2,J+1))=max(Uout.emin(sub2ind(Nu,I+2,J+1)),Hridge(k));
      Vout.emin(sub2ind(Nv,I+1,J+2))=max(Vout.emin(sub2ind(Nv,I+1,J+2)),Hridge(k));
      Vout.eave(sub2ind(Nv,I+1,J))=max(moveAveMax*Vout.eave(sub2ind(Nv,I+1,J)),HridgeAve(k));
      Uout.eave(sub2ind(Nu,I+2,J))=max(moveAveMax*Uout.eave(sub2ind(Nu,I+2,J)),HridgeAve(k));
      Uout.eave(sub2ind(Nu,I+2,J+1))=max(moveAveMax*Uout.eave(sub2ind(Nu,I+2,J+1)),HridgeAve(k));
      Vout.eave(sub2ind(Nv,I+1,J+2))=max(moveAveMax*Vout.eave(sub2ind(Nv,I+1,J+2)),HridgeAve(k));
      Vout.emax(sub2ind(Nv,I+1,J))=max(moveAveMax*Vout.emax(sub2ind(Nv,I+1,J)),HridgeMax(k));
      Uout.emax(sub2ind(Nu,I+2,J))=max(moveAveMax*Uout.emax(sub2ind(Nu,I+2,J)),HridgeMax(k));
      Uout.emax(sub2ind(Nu,I+2,J+1))=max(moveAveMax*Uout.emax(sub2ind(Nu,I+2,J+1)),HridgeMax(k));
      Vout.emax(sub2ind(Nv,I+1,J+2))=max(moveAveMax*Vout.emax(sub2ind(Nv,I+1,J+2)),HridgeMax(k));
      Hout.emin(sub2ind(Nh,I+1,J))=HoppoMin(k);
      Hout.emin(sub2ind(Nh,I+1,J+1))=HoppoMin(k);
      Uout.emin(sub2ind(Nu,I+1,J))=HoppoMin(k);
      Uout.emin(sub2ind(Nu,I+1,J+1))=HoppoMin(k);
      Vout.emin(sub2ind(Nv,I+1,J+1))=HoppoMin(k);
      Hout.eave(sub2ind(Nh,I+1,J))=(Hin.eave(sub2ind(Nh,I,J))+Hin.eave(sub2ind(Nh,I,J+1)))/2;
      Hout.eave(sub2ind(Nh,I+1,J+1))=(Hin.eave(sub2ind(Nh,I,J))+Hin.eave(sub2ind(Nh,I,J+1)))/2;
      Hout.emax(sub2ind(Nh,I+1,J))=HoppoMin(k);
      Hout.emax(sub2ind(Nh,I+1,J+1))=HoppoMin(k);
      msk(sub2ind(Nh,I+1,J))=1;
      msk(sub2ind(Nh,I+1,J+1))=1;
      mskC(k)=1;
  i=2*ic-1; j=2*jc-1;
  [i,j]=find( mskEW & ( Vin.emin(i,j+1)>Vin.emin(i+1,j+1) ...
            | ( Vin.emin(i,j+1)==Vin.emin(i+1,j+1) ...
              & ( Hin.emin(i+1,j)+Hin.emin(i+1,j+1)<Hin.emin(i,j)+Hin.emin(i,j+1) ...
                | Uin.emin(i+2,j)+Uin.emin(i+2,j+1)<Uin.emin(i,j)+Uin.emin(i,j+1) ...
            ) ) ) );
  I=2*i-1;J=2*j-1; k=sub2ind(size(Hridge),i,j);
      Vout.emin(sub2ind(Nv,I,J))=max(Vout.emin(sub2ind(Nv,I,J)),Hridge(k));
      Uout.emin(sub2ind(Nu,I,J))=max(Uout.emin(sub2ind(Nu,I,J)),Hridge(k));
      Uout.emin(sub2ind(Nu,I,J+1))=max(Uout.emin(sub2ind(Nu,I,J+1)),Hridge(k));
      Vout.emin(sub2ind(Nv,I,J+2))=max(Vout.emin(sub2ind(Nv,I,J+2)),Hridge(k));
      Vout.eave(sub2ind(Nv,I,J))=max(moveAveMax*Vout.eave(sub2ind(Nv,I,J)),HridgeAve(k));
      Uout.eave(sub2ind(Nu,I,J))=max(moveAveMax*Uout.eave(sub2ind(Nu,I,J)),HridgeAve(k));
      Uout.eave(sub2ind(Nu,I,J+1))=max(moveAveMax*Uout.eave(sub2ind(Nu,I,J+1)),HridgeAve(k));
      Vout.eave(sub2ind(Nv,I,J+2))=max(moveAveMax*Vout.eave(sub2ind(Nv,I,J+2)),HridgeAve(k));
      Vout.emax(sub2ind(Nv,I,J))=max(moveAveMax*Vout.emax(sub2ind(Nv,I,J)),HridgeMax(k));
      Uout.emax(sub2ind(Nu,I,J))=max(moveAveMax*Uout.emax(sub2ind(Nu,I,J)),HridgeMax(k));
      Uout.emax(sub2ind(Nu,I,J+1))=max(moveAveMax*Uout.emax(sub2ind(Nu,I,J+1)),HridgeMax(k));
      Vout.emax(sub2ind(Nv,I,J+2))=max(moveAveMax*Vout.emax(sub2ind(Nv,I,J+2)),HridgeMax(k));
      Hout.emin(sub2ind(Nh,I,J))=HoppoMin(k);
      Hout.emin(sub2ind(Nh,I,J+1))=HoppoMin(k);
      Uout.emin(sub2ind(Nu,I+1,J))=HoppoMin(k);
      Uout.emin(sub2ind(Nu,I+1,J+1))=HoppoMin(k);
      Vout.emin(sub2ind(Nv,I,J+1))=HoppoMin(k);
      Hout.eave(sub2ind(Nh,I,J))=(Hin.eave(sub2ind(Nh,I+1,J))+Hin.eave(sub2ind(Nh,I+1,J+1)))/2;
      Hout.eave(sub2ind(Nh,I,J+1))=(Hin.eave(sub2ind(Nh,I+1,J))+Hin.eave(sub2ind(Nh,I+1,J+1)))/2;
      Hout.emax(sub2ind(Nh,I,J))=HoppoMin(k);
      Hout.emax(sub2ind(Nh,I,J+1))=HoppoMin(k);
      msk(sub2ind(Nh,I,J))=1;
      msk(sub2ind(Nh,I,J+1))=1;
      mskC(k)=1;
  i=2*ic-1; j=2*jc-1;
  [i,j]=find( mskEW & mskC==0 & ( Vin.emin(i,j+1)==Vin.emin(i+1,j+1) ) );
  I=2*i-1;J=2*j-1; k=sub2ind(size(Hridge),i,j);
      Uout.emin(sub2ind(Nu,I,J))=max(Uout.emin(sub2ind(Nu,I,J)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+2,J))=max(Uout.emin(sub2ind(Nu,I+2,J)),Hridge(k));
      Uout.emin(sub2ind(Nu,I,J+1))=max(Uout.emin(sub2ind(Nu,I,J+1)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+2,J+1))=max(Uout.emin(sub2ind(Nu,I+2,J+1)),Hridge(k));
      Vout.emin(sub2ind(Nv,I,J))=max(Vout.emin(sub2ind(Nv,I,J)),Hridge(k));
      Vout.emin(sub2ind(Nv,I,J+2))=max(Vout.emin(sub2ind(Nv,I,J+2)),Hridge(k));
      Vout.emin(sub2ind(Nv,I+1,J))=max(Vout.emin(sub2ind(Nv,I+1,J)),Hridge(k));
      Vout.emin(sub2ind(Nv,I+1,J+2))=max(Vout.emin(sub2ind(Nv,I+1,J+2)),Hridge(k));
      Uout.emin(sub2ind(Nu,I+1,J))=Hridge(k);
      Uout.emin(sub2ind(Nu,I+1,J+1))=Hridge(k);
      Vout.emin(sub2ind(Nv,I,J+1))=Hridge(k);
      Vout.emin(sub2ind(Nv,I+1,J+1))=Hridge(k);
      Hout.emin(sub2ind(Nh,I,J))=Hridge(k)*NaN;
      Hout.emin(sub2ind(Nh,I+1,J))=Hridge(k)*NaN;
      Hout.emin(sub2ind(Nh,I,J+1))=Hridge(k)*NaN;
      Hout.emin(sub2ind(Nh,I+1,J+1))=Hridge(k)*NaN;
      msk(sub2ind(Nh,I,J))=1;
      msk(sub2ind(Nh,I+1,J))=1;
      msk(sub2ind(Nh,I,J+1))=1;
      msk(sub2ind(Nh,I+1,J+1))=1;

 end
 if ~exist('msg','var')
  locdisp(sprintf('  centralRidge: %i points (%5.2f%%) "pushed"', ...
   nncnt,100*nncnt/(nxf*nyf) ))
 else
  if nncnt>0
   locdisp(sprintf('  centralRidge: %i points (%5.2f%%) "pushed" (%s)', ...
   nncnt,100*nncnt/(nxf*nyf) ,msg))
  end
 end
% -----------------------------------------------------------------------------
function [Uout,Vout,msk] = lowerTallestButtress(Uin,Vin,msg)
% Lowers any buttress to the highest neighbours since this does not
% affect connectivity
 nxf=size(Vin,1); nyf=size(Uin,2); msk=NaN*zeros(nxf,nyf);
 ic=1:nxf/2;jc=1:nyf/2;
 Uout=Uin;Vout=Vin;
 global lowMem
 if lowMem
  nncnt=0;
  for J=jc
   j=2*J-1;
   for I=ic 
    i=2*I-1;
    ncnt=0;
    % W
    Ho=max([Uin(i+1,j) Vin(i+1,j+1) Uin(i+1,j+1)]);
    if Vout(i,j+1)>Ho
     Vout(i,j+1)=Ho; ncnt=ncnt+1; msk(i,j:j+1)=1;
    end
    % S
    Ho=max([Vin(i,j+1) Uin(i+1,j+1) Vin(i+1,j+1)]);
    if Uout(i+1,j)>Ho
     Uout(i+1,j)=Ho; ncnt=ncnt+1; msk(i:i+1,j)=1;
    end
    % E
    Ho=max([Uin(i+1,j) Vin(i,j+1) Uin(i+1,j+1)]);
    if Vout(i+1,j+1)>Ho
     Vout(i+1,j+1)=Ho; ncnt=ncnt+1; msk(i+1,j:j+1)=1;
    end
    % N
    Ho=max([Vin(i,j+1) Uin(i+1,j) Vin(i+1,j+1)]);
    if Uout(i+1,j+1)>Ho
     Uout(i+1,j+1)=Ho; ncnt=ncnt+1; msk(i:i+1,j+1)=1;
    end
    if ncnt>1
     error('lowerTallestButtress: This should be impossible!')
    end
    nncnt=nncnt+ncnt;
   end % I
  end % J
 else % !lowMem
  % W
  i=2*ic-1; j=2*jc-1;
  Ho=max(max(Uin(i+1,j),Vin(i+1,j+1)),Uin(i+1,j+1));
  [i,j]=find(Vout(i,j+1)>Ho); I=2*i-1;J=2*j-1;
  Vout(sub2ind(size(Vout),I,J+1))=Ho(sub2ind(size(Ho),i,j));
  msk(sub2ind(size(msk),I,J))=1; msk(sub2ind(size(msk),I,J+1))=1;
  % E
  i=2*ic-1; j=2*jc-1;
  Ho=max(max(Uin(i+1,j),Vin(i,j+1)),Uin(i+1,j+1));
  [i,j]=find(Vout(i+1,j+1)>Ho); I=2*i-1;J=2*j-1;
  Vout(sub2ind(size(Vout),I+1,J+1))=Ho(sub2ind(size(Ho),i,j));
  msk(sub2ind(size(msk),I+1,J))=1; msk(sub2ind(size(msk),I+1,J+1))=1;
  % S
  i=2*ic-1; j=2*jc-1;
  Ho=max(max(Vin(i,j+1),Uin(i+1,j+1)),Vin(i+1,j+1));
  [i,j]=find(Uout(i+1,j)>Ho); I=2*i-1;J=2*j-1;
  Uout(sub2ind(size(Uout),I+1,J))=Ho(sub2ind(size(Ho),i,j));
  msk(sub2ind(size(msk),I,J))=1; msk(sub2ind(size(msk),I+1,J))=1;
  % N
  i=2*ic-1; j=2*jc-1;
  Ho=max(max(Vin(i,j+1),Uin(i+1,j)),Vin(i+1,j+1));
  [i,j]=find(Uout(i+1,j+1)>Ho); I=2*i-1;J=2*j-1;
  Uout(sub2ind(size(Uout),I+1,J+1))=Ho(sub2ind(size(Ho),i,j));
  msk(sub2ind(size(msk),I,J+1))=1; msk(sub2ind(size(msk),I+1,J+1))=1;
  nncnt=length(find(~isnan(msk(:))));
 end
 if ~exist('msg','var')
  locdisp(sprintf('  lowerTallestButtress: %i points (%5.2f%%) "lowered"', ...
   nncnt,100*nncnt/(nxf*nyf) ))
 else
  if nncnt>0
   locdisp(sprintf('  lowerTallestButtress: %i points (%5.2f%%) "lowered" (%s)', ...
   nncnt,100*nncnt/(nxf*nyf) ,msg))
  end
 end
% -----------------------------------------------------------------------------
function [Uout,Vout,msk] = expandCorners(Uin,Vin)
 [nxf,nyf]=size(Vin.emin); nyf=nyf-1; msk=NaN*zeros(nxf,nyf);
 ic=1:nxf/2;jc=1:nyf/2;
 Uout=Uin;Vout=Vin;
 global lowMem
 if lowMem
  nncnt=0;
  for J=jc
   j=2*J-1;
   for I=ic 
    i=2*I-1;
    ncnt=0;
    % Outer pairs max
    HswOMx=max(Uin.emin(i,j),Vin.emin(i,j));
    HseOMx=max(Uin.emin(i+2,j),Vin.emin(i+1,j));
    HnwOMx=max(Uin.emin(i,j+1),Vin.emin(i,j+2));
    HneOMx=max(Uin.emin(i+2,j+1),Vin.emin(i+1,j+2));
    % Inner pairs min
    HswIMn=min(Uin.emin(i+1,j),Vin.emin(i,j+1));
    HseIMn=min(Uin.emin(i+1,j),Vin.emin(i+1,j+1));
    HnwIMn=min(Uin.emin(i+1,j+1),Vin.emin(i,j+1));
    HneIMn=min(Uin.emin(i+1,j+1),Vin.emin(i+1,j+1));
    % Inner pairs max
    HswIMx=max(Uin.emin(i+1,j),Vin.emin(i,j+1));
    HseIMx=max(Uin.emin(i+1,j),Vin.emin(i+1,j+1));
    HnwIMx=max(Uin.emin(i+1,j+1),Vin.emin(i,j+1));
    HneIMx=max(Uin.emin(i+1,j+1),Vin.emin(i+1,j+1));
    % Deepest outer corner
    crnrMins=min([HswOMx,HseOMx,HnwOMx,HneOMx]);
    crnrPenMins=penmin([HswOMx,HseOMx,HnwOMx,HneOMx]);
    if crnrPenMins>crnrMins % Only for unique deepest corners
     if HswOMx==crnrMins && HswOMx<HswIMn
      Uout.emin(i+1,j:j+1)=HswOMx; Vout.emin(i:i+1,j+1)=HswOMx; % Set interior +
      newVal=min(HseIMx,HnwIMx); % Possible diagonal ridge forming corner
      Uout.emin(i,j+1)=max(Uout.emin(i,j+1),newVal);
      Vout.emin(i,j+2)=max(Vout.emin(i,j+2),newVal);
      Vout.emin(i+1,j+2)=max(Vout.emin(i+1,j+2),newVal);
      Uout.emin(i+2,j+1)=max(Uout.emin(i+2,j+1),newVal);
      Uout.emin(i+2,j)=max(Uout.emin(i+2,j),newVal);
      Vout.emin(i+1,j)=max(Vout.emin(i+1,j),newVal);
      ncnt=ncnt+1; msk(i,j)=1;
     end % SW
     if HneOMx==crnrMins && HneOMx<HneIMn
      Uout.emin(i+1,j:j+1)=HneOMx; Vout.emin(i:i+1,j+1)=HneOMx; % Set interior +
      newVal=min(HseIMx,HnwIMx); % Possible diagonal ridge forming corner
      Uout.emin(i+2,j)=max(Uout.emin(i+2,j),newVal);
      Vout.emin(i+1,j)=max(Vout.emin(i+1,j),newVal);
      Vout.emin(i,j)=max(Vout.emin(i,j),newVal);
      Uout.emin(i,j)=max(Uout.emin(i,j),newVal);
      Uout.emin(i,j+1)=max(Uout.emin(i,j+1),newVal);
      Vout.emin(i,j+2)=max(Vout.emin(i,j+2),newVal);
      ncnt=ncnt+1; msk(i+1,j+1)=1;
     end % NE
     if HnwOMx==crnrMins && HnwOMx<HnwIMn
      Uout.emin(i+1,j:j+1)=HnwOMx; Vout.emin(i:i+1,j+1)=HnwOMx; % Set interior +
      newVal=min(HneIMx,HswIMx); % Possible diagonal ridge forming corner
      Vout.emin(i+1,j+2)=max(Vout.emin(i+1,j+2),newVal);
      Uout.emin(i+2,j+1)=max(Uout.emin(i+2,j+1),newVal);
      Uout.emin(i+2,j)=max(Uout.emin(i+2,j),newVal);
      Vout.emin(i+1,j)=max(Vout.emin(i+1,j),newVal);
      Vout.emin(i,j)=max(Vout.emin(i,j),newVal);
      Uout.emin(i,j)=max(Uout.emin(i,j),newVal);
      ncnt=ncnt+1; msk(i,j+1)=1;
     end % NW
     if HseOMx==crnrMins && HseOMx<HseIMn
      Uout.emin(i+1,j:j+1)=HseOMx; Vout.emin(i:i+1,j+1)=HseOMx; % Set interior +
      newVal=min(HneIMx,HswIMx); % Possible diagonal ridge forming corner
      Vout.emin(i,j)=max(Vout.emin(i,j),newVal);
      Uout.emin(i,j)=max(Uout.emin(i,j),newVal);
      Uout.emin(i,j+1)=max(Uout.emin(i,j+1),newVal);
      Vout.emin(i,j+2)=max(Vout.emin(i,j+2),newVal);
      Vout.emin(i+1,j+2)=max(Vout.emin(i+1,j+2),newVal);
      Uout.emin(i+2,j+1)=max(Uout.emin(i+2,j+1),newVal);
      ncnt=ncnt+1; msk(i+1,j)=1;
     end % SE
    end
    if ncnt>1
      locdisp('expandCorners: This should NEVER happen (while expanding big corners)!')
    end
    nncnt=nncnt+ncnt;
   end % I
  end % J
 else % !lowMem
  Nu=size(Uout.emin);Nv=size(Vout.emin);Nh=min(Nu,Nv);
  i=2*ic-1; j=2*jc-1;
    % Outer pairs max
    HswOMx=max(Uin.emin(i,j),Vin.emin(i,j));
    HseOMx=max(Uin.emin(i+2,j),Vin.emin(i+1,j));
    HnwOMx=max(Uin.emin(i,j+1),Vin.emin(i,j+2));
    HneOMx=max(Uin.emin(i+2,j+1),Vin.emin(i+1,j+2));
    % Inner pairs min
    HswIMn=min(Uin.emin(i+1,j),Vin.emin(i,j+1));
    HseIMn=min(Uin.emin(i+1,j),Vin.emin(i+1,j+1));
    HnwIMn=min(Uin.emin(i+1,j+1),Vin.emin(i,j+1));
    HneIMn=min(Uin.emin(i+1,j+1),Vin.emin(i+1,j+1));
    % Inner pairs max
    HswIMx=max(Uin.emin(i+1,j),Vin.emin(i,j+1));
    HseIMx=max(Uin.emin(i+1,j),Vin.emin(i+1,j+1));
    HnwIMx=max(Uin.emin(i+1,j+1),Vin.emin(i,j+1));
    HneIMx=max(Uin.emin(i+1,j+1),Vin.emin(i+1,j+1));
    % Deepest outer corner
    crnrMins=min(min(HswOMx,HseOMx),min(HnwOMx,HneOMx));
    crnrPenMins=penminV(HswOMx,HseOMx,HnwOMx,HneOMx);
  % SW
  [i,j]=find( crnrPenMins>crnrMins & HswOMx==crnrMins & HswOMx<HswIMn ); I=2*i-1;J=2*j-1;
  k=sub2ind(size(crnrMins),i,j);
      Uout.emin(sub2ind(Nu,I+1,J))=HswOMx(k); % Set interior +
      Uout.emin(sub2ind(Nu,I+1,J+1))=HswOMx(k); % Set interior +
      Vout.emin(sub2ind(Nv,I,J+1))=HswOMx(k); % Set interior +
      Vout.emin(sub2ind(Nv,I+1,J+1))=HswOMx(k); % Set interior +
      newVal=min(HseIMx(k),HnwIMx(k)); % Possible diagonal ridge forming corner
      Uout.emin(sub2ind(Nu,I,J+1))=max(Uout.emin(sub2ind(Nu,I,J+1)),newVal);
      Vout.emin(sub2ind(Nv,I,J+2))=max(Vout.emin(sub2ind(Nv,I,J+2)),newVal);
      Vout.emin(sub2ind(Nv,I+1,J+2))=max(Vout.emin(sub2ind(Nv,I+1,J+2)),newVal);
      Uout.emin(sub2ind(Nu,I+2,J+1))=max(Uout.emin(sub2ind(Nu,I+2,J+1)),newVal);
      Uout.emin(sub2ind(Nu,I+2,J))=max(Uout.emin(sub2ind(Nu,I+2,J)),newVal);
      Vout.emin(sub2ind(Nv,I+1,J))=max(Vout.emin(sub2ind(Nv,I+1,J)),newVal);
      newVal=(HseIMx(k)+HnwIMx(k))/2; % Possible diagonal ridge forming corner
      Uout.eave(sub2ind(Nu,I,J+1))=max(Uout.eave(sub2ind(Nu,I,J+1)),newVal);
      Vout.eave(sub2ind(Nv,I,J+2))=max(Vout.eave(sub2ind(Nv,I,J+2)),newVal);
      Vout.eave(sub2ind(Nv,I+1,J+2))=max(Vout.eave(sub2ind(Nv,I+1,J+2)),newVal);
      Uout.eave(sub2ind(Nu,I+2,J+1))=max(Uout.eave(sub2ind(Nu,I+2,J+1)),newVal);
      Uout.eave(sub2ind(Nu,I+2,J))=max(Uout.eave(sub2ind(Nu,I+2,J)),newVal);
      Vout.eave(sub2ind(Nv,I+1,J))=max(Vout.eave(sub2ind(Nv,I+1,J)),newVal);
      newVal=max(HseIMx(k),HnwIMx(k)); % Possible diagonal ridge forming corner
      Uout.emax(sub2ind(Nu,I,J+1))=max(Uout.emax(sub2ind(Nu,I,J+1)),newVal);
      Vout.emax(sub2ind(Nv,I,J+2))=max(Vout.emax(sub2ind(Nv,I,J+2)),newVal);
      Vout.emax(sub2ind(Nv,I+1,J+2))=max(Vout.emax(sub2ind(Nv,I+1,J+2)),newVal);
      Uout.emax(sub2ind(Nu,I+2,J+1))=max(Uout.emax(sub2ind(Nu,I+2,J+1)),newVal);
      Uout.emax(sub2ind(Nu,I+2,J))=max(Uout.emax(sub2ind(Nu,I+2,J)),newVal);
      Vout.emax(sub2ind(Nv,I+1,J))=max(Vout.emax(sub2ind(Nv,I+1,J)),newVal);
      msk(sub2ind(Nh,I,J))=1;
  % NE
  [i,j]=find( crnrPenMins>crnrMins & HneOMx==crnrMins & HneOMx<HneIMn ); I=2*i-1;J=2*j-1;
  k=sub2ind(size(crnrMins),i,j);
      Uout.emin(sub2ind(Nu,I+1,J))=HneOMx(k); % Set interior +
      Uout.emin(sub2ind(Nu,I+1,J+1))=HneOMx(k); % Set interior +
      Vout.emin(sub2ind(Nv,I,J+1))=HneOMx(k); % Set interior +
      Vout.emin(sub2ind(Nv,I+1,J+1))=HneOMx(k); % Set interior +
      newVal=min(HseIMx(k),HnwIMx(k)); % Possible diagonal ridge forming corner
      Uout.emin(sub2ind(Nu,I+2,J))=max(Uout.emin(sub2ind(Nu,I+2,J)),newVal);
      Vout.emin(sub2ind(Nv,I+1,J))=max(Vout.emin(sub2ind(Nv,I+1,J)),newVal);
      Vout.emin(sub2ind(Nv,I,J))=max(Vout.emin(sub2ind(Nv,I,J)),newVal);
      Uout.emin(sub2ind(Nu,I,J))=max(Uout.emin(sub2ind(Nu,I,J)),newVal);
      Uout.emin(sub2ind(Nu,I,J+1))=max(Uout.emin(sub2ind(Nu,I,J+1)),newVal);
      Vout.emin(sub2ind(Nv,I,J+2))=max(Vout.emin(sub2ind(Nv,I,J+2)),newVal);
      newVal=(HseIMx(k)+HnwIMx(k))/2; % Possible diagonal ridge forming corner
      Uout.eave(sub2ind(Nu,I+2,J))=max(Uout.eave(sub2ind(Nu,I+2,J)),newVal);
      Vout.eave(sub2ind(Nv,I+1,J))=max(Vout.eave(sub2ind(Nv,I+1,J)),newVal);
      Vout.eave(sub2ind(Nv,I,J))=max(Vout.eave(sub2ind(Nv,I,J)),newVal);
      Uout.eave(sub2ind(Nu,I,J))=max(Uout.eave(sub2ind(Nu,I,J)),newVal);
      Uout.eave(sub2ind(Nu,I,J+1))=max(Uout.eave(sub2ind(Nu,I,J+1)),newVal);
      Vout.eave(sub2ind(Nv,I,J+2))=max(Vout.eave(sub2ind(Nv,I,J+2)),newVal);
      newVal=max(HseIMx(k),HnwIMx(k)); % Possible diagonal ridge forming corner
      Uout.emax(sub2ind(Nu,I+2,J))=max(Uout.emax(sub2ind(Nu,I+2,J)),newVal);
      Vout.emax(sub2ind(Nv,I+1,J))=max(Vout.emax(sub2ind(Nv,I+1,J)),newVal);
      Vout.emax(sub2ind(Nv,I,J))=max(Vout.emax(sub2ind(Nv,I,J)),newVal);
      Uout.emax(sub2ind(Nu,I,J))=max(Uout.emax(sub2ind(Nu,I,J)),newVal);
      Uout.emax(sub2ind(Nu,I,J+1))=max(Uout.emax(sub2ind(Nu,I,J+1)),newVal);
      Vout.emax(sub2ind(Nv,I,J+2))=max(Vout.emax(sub2ind(Nv,I,J+2)),newVal);
      msk(sub2ind(Nh,I+1,J+1))=1;
  % NW
  [i,j]=find( crnrPenMins>crnrMins & HnwOMx==crnrMins & HnwOMx<HnwIMn ); I=2*i-1;J=2*j-1;
  k=sub2ind(size(crnrMins),i,j);
      Uout.emin(sub2ind(Nu,I+1,J))=HnwOMx(k); % Set interior +
      Uout.emin(sub2ind(Nu,I+1,J+1))=HnwOMx(k); % Set interior +
      Vout.emin(sub2ind(Nv,I,J+1))=HnwOMx(k); % Set interior +
      Vout.emin(sub2ind(Nv,I+1,J+1))=HnwOMx(k); % Set interior +
      newVal=min(HneIMx(k),HswIMx(k)); % Possible diagonal ridge forming corner
      Vout.emin(sub2ind(Nv,I+1,J+2))=max(Vout.emin(sub2ind(Nv,I+1,J+2)),newVal);
      Uout.emin(sub2ind(Nu,I+2,J+1))=max(Uout.emin(sub2ind(Nu,I+2,J+1)),newVal);
      Uout.emin(sub2ind(Nu,I+2,J))=max(Uout.emin(sub2ind(Nu,I+2,J)),newVal);
      Vout.emin(sub2ind(Nv,I+1,J))=max(Vout.emin(sub2ind(Nv,I+1,J)),newVal);
      Vout.emin(sub2ind(Nv,I,J))=max(Vout.emin(sub2ind(Nv,I,J)),newVal);
      Uout.emin(sub2ind(Nu,I,J))=max(Uout.emin(sub2ind(Nu,I,J)),newVal);
      newVal=(HneIMx(k)+HswIMx(k))/2; % Possible diagonal ridge forming corner
      Vout.eave(sub2ind(Nv,I+1,J+2))=max(Vout.eave(sub2ind(Nv,I+1,J+2)),newVal);
      Uout.eave(sub2ind(Nu,I+2,J+1))=max(Uout.eave(sub2ind(Nu,I+2,J+1)),newVal);
      Uout.eave(sub2ind(Nu,I+2,J))=max(Uout.eave(sub2ind(Nu,I+2,J)),newVal);
      Vout.eave(sub2ind(Nv,I+1,J))=max(Vout.eave(sub2ind(Nv,I+1,J)),newVal);
      Vout.eave(sub2ind(Nv,I,J))=max(Vout.eave(sub2ind(Nv,I,J)),newVal);
      Uout.eave(sub2ind(Nu,I,J))=max(Uout.eave(sub2ind(Nu,I,J)),newVal);
      newVal=max(HneIMx(k),HswIMx(k)); % Possible diagonal ridge forming corner
      Vout.emax(sub2ind(Nv,I+1,J+2))=max(Vout.emax(sub2ind(Nv,I+1,J+2)),newVal);
      Uout.emax(sub2ind(Nu,I+2,J+1))=max(Uout.emax(sub2ind(Nu,I+2,J+1)),newVal);
      Uout.emax(sub2ind(Nu,I+2,J))=max(Uout.emax(sub2ind(Nu,I+2,J)),newVal);
      Vout.emax(sub2ind(Nv,I+1,J))=max(Vout.emax(sub2ind(Nv,I+1,J)),newVal);
      Vout.emax(sub2ind(Nv,I,J))=max(Vout.emax(sub2ind(Nv,I,J)),newVal);
      Uout.emax(sub2ind(Nu,I,J))=max(Uout.emax(sub2ind(Nu,I,J)),newVal);
      msk(sub2ind(Nh,I,J+1))=1;
  % SE
  [i,j]=find( crnrPenMins>crnrMins & HseOMx==crnrMins & HseOMx<HseIMn ); I=2*i-1;J=2*j-1;
  k=sub2ind(size(crnrMins),i,j);
      Uout.emin(sub2ind(Nu,I+1,J))=HseOMx(k); % Set interior +
      Uout.emin(sub2ind(Nu,I+1,J+1))=HseOMx(k); % Set interior +
      Vout.emin(sub2ind(Nv,I,J+1))=HseOMx(k); % Set interior +
      Vout.emin(sub2ind(Nv,I+1,J+1))=HseOMx(k); % Set interior +
      newVal=min(HneIMx(k),HswIMx(k)); % Possible diagonal ridge forming corner
      Vout.emin(sub2ind(Nv,I,J))=max(Vout.emin(sub2ind(Nv,I,J)),newVal);
      Uout.emin(sub2ind(Nu,I,J))=max(Uout.emin(sub2ind(Nu,I,J)),newVal);
      Uout.emin(sub2ind(Nu,I,J+1))=max(Uout.emin(sub2ind(Nu,I,J+1)),newVal);
      Vout.emin(sub2ind(Nv,I,J+2))=max(Vout.emin(sub2ind(Nv,I,J+2)),newVal);
      Vout.emin(sub2ind(Nv,I+1,J+2))=max(Vout.emin(sub2ind(Nv,I+1,J+2)),newVal);
      Uout.emin(sub2ind(Nu,I+2,J+1))=max(Uout.emin(sub2ind(Nu,I+2,J+1)),newVal);
      newVal=(HneIMx(k)+HswIMx(k))/2; % Possible diagonal ridge forming corner
      Vout.eave(sub2ind(Nv,I,J))=max(Vout.eave(sub2ind(Nv,I,J)),newVal);
      Uout.eave(sub2ind(Nu,I,J))=max(Uout.eave(sub2ind(Nu,I,J)),newVal);
      Uout.eave(sub2ind(Nu,I,J+1))=max(Uout.eave(sub2ind(Nu,I,J+1)),newVal);
      Vout.eave(sub2ind(Nv,I,J+2))=max(Vout.eave(sub2ind(Nv,I,J+2)),newVal);
      Vout.eave(sub2ind(Nv,I+1,J+2))=max(Vout.eave(sub2ind(Nv,I+1,J+2)),newVal);
      Uout.eave(sub2ind(Nu,I+2,J+1))=max(Uout.eave(sub2ind(Nu,I+2,J+1)),newVal);
      newVal=max(HneIMx(k),HswIMx(k)); % Possible diagonal ridge forming corner
      Vout.emax(sub2ind(Nv,I,J))=max(Vout.emax(sub2ind(Nv,I,J)),newVal);
      Uout.emax(sub2ind(Nu,I,J))=max(Uout.emax(sub2ind(Nu,I,J)),newVal);
      Uout.emax(sub2ind(Nu,I,J+1))=max(Uout.emax(sub2ind(Nu,I,J+1)),newVal);
      Vout.emax(sub2ind(Nv,I,J+2))=max(Vout.emax(sub2ind(Nv,I,J+2)),newVal);
      Vout.emax(sub2ind(Nv,I+1,J+2))=max(Vout.emax(sub2ind(Nv,I+1,J+2)),newVal);
      Uout.emax(sub2ind(Nu,I+2,J+1))=max(Uout.emax(sub2ind(Nu,I+2,J+1)),newVal);
      msk(sub2ind(Nh,I+1,J))=1;

  nncnt=length(find(~isnan(msk(:))));
 end
 if ~exist('msg','var')
  locdisp(sprintf('  expandCorners: %i points (%5.2f%%) "filled"', ...
   nncnt,100*nncnt/(4*nxf*nyf) ))
 else
  if nncnt>0
   locdisp(sprintf('  expandCorners: %i points (%5.2f%%) "filled" (%s)', ...
    nncnt,100*nncnt/(4*nxf*nyf),msg))
  end
 end
% -----------------------------------------------------------------------------
function [Hout,msk] = fillPotHoles(Hin,Uin,Vin)
 [nx,ny]=size(Hin); msk=NaN*Hin;
 Hout=Hin;
 ncnt=0;
 global lowMem
 if lowMem
  for j=1:ny
   for i=1:nx
    He=min([Uin(i,j),Uin(i+1,j),Vin(i,j),Vin(i,j+1)]);
    if Hin(i,j)<He
     Hout(i,j)=He; ncnt=ncnt+1; msk(i,j)=1;
    end
   end
  end
 else % !lowMem
  i=1:nx; j=1:ny;
  He=min(min(Uin(i,j),Uin(i+1,j)),min(Vin(i,j),Vin(i,j+1)));
  Hout=max(Hin,He);
  msk( find(Hout>Hin) )=1; nncnt=length(find(~isnan(msk(:))));
 end
 locdisp(sprintf('  fillPotHoles: %i points (%5.2f%%) "filled"', ...
   ncnt,100*ncnt/(nx*ny) ))
% -----------------------------------------------------------------------------
function [Hout,msk] = boundHbyUV(Hin,Uin,Vin)
 [nx,ny]=size(Hin); msk=NaN*Hin;
 Hout=Hin;
 ncnt=0;
 global lowMem
 if lowMem
  for j=1:ny
   for i=1:nx
    He=min([Uin(i,j),Uin(i+1,j),Vin(i,j),Vin(i,j+1)]);
    if Hin(i,j)>He
     Hout(i,j)=He; ncnt=ncnt+1; msk(i,j)=1;
    end
   end
  end
 else % !lowMem
  i=1:nx; j=1:ny;
  He=min(min(Uin(i,j),Uin(i+1,j)),min(Vin(i,j),Vin(i,j+1)));
  Hout=min(Hin,He);
  msk( find(Hout>Hin) )=1; nncnt=length(find(~isnan(msk(:))));
 end
 locdisp(sprintf('  boundHbyUV: %i points (%5.2f%%) "filled"', ...
   ncnt,100*ncnt/(nx*ny) ))
% -----------------------------------------------------------------------------
function [Uout,Vout,msk] = limitNSConnections(Deepest,Uin,Vin)
 [nxf,nyf]=size(Vin); nyf=nyf-1; msk=NaN*zeros(nxf,nyf);
 ic=1:nxf/2;jc=1:nyf/2;
 Uout=Uin;Vout=Vin;
 global lowMem
 if 1%lowMem
  nncnt=0;
  for J=jc
   j=2*J-1;
   for I=ic 
    i=2*I-1;
    ncnt=0;
    % Coarse width ridges
    Hew_s=min(Vin(i,j  ),Vin(i+1,j));
    Hew_n=min(Vin(i,j+2),Vin(i+1,j+2));
    Hns_w=min(Uin(i  ,j),Uin(i  ,j+1));
    Hns_e=min(Uin(i+2,j),Uin(i+2,j+1));
 %  Uout(i+1,j:j+1)=min(Uout(i+1,j:j+1),min([Hew_s Hew_n Hns_w Hns_e])); Vout(i:i+1,j+1)=min(Vout(i:i+1,j+1),min([Hew_s Hew_n Hns_w Hns_e])); % Throw away DEBUGging
    if Deepest.NS(I,J)>max(Hew_n,Hew_s) % Limit NS connectivity
     if Hew_n>Hew_s
      Vout(i:i+1,j+2)=max(Vout(i:i+1,j+2),Deepest.NS(I,J)); % Block N
      msk(i,j+1)=1;
     elseif Hew_n<Hew_s
      Vout(i:i+1,j)=max(Vout(i:i+1,j),Deepest.NS(I,J)); % Block S
      msk(i+1,j+1)=1;
     else
      Vout(i:i+1,j+2)=max(Vout(i:i+1,j+2),Deepest.NS(I,J)); % Block N
      Vout(i:i+1,j)=max(Vout(i:i+1,j),Deepest.NS(I,J)); % Block S
      msk(i:i+1,j+1)=1;
     end
     ncnt=ncnt+1;
    end
    if Deepest.EW(I,J)>max(Hns_e,Hns_w) % Limit EW connectivity
     if Hns_e>Hns_w
      Uout(i+2,j:j+1)=max(Uout(i+2,j:j+1),Deepest.EW(I,J)); % Block E
      msk(i,j)=1;
     elseif Hns_e<Hns_w
      Uout(i,j:j+1)=max(Uout(i,j:j+1),Deepest.EW(I,J)); % Block W
      msk(i+1,j)=1;
     else
      Uout(i+2,j:j+1)=max(Uout(i+2,j:j+1),Deepest.EW(I,J)); % Block E
      Uout(i,j:j+1)=max(Uout(i,j:j+1),Deepest.EW(I,J)); % Block W
      msk(i:i+1,j)=1;
     end
     ncnt=ncnt+1;
    end
    if ncnt>2
      locdisp('limitNSConnections: This should NEVER happen!')
    end
    nncnt=nncnt+ncnt;
   end % I
  end % J
 else % !lowMem
  Nu=size(Uin);Nv=size(Vin);Nh=min(Nu,Nv);Ch=Nh/2;
  i=2*ic-1; j=2*jc-1; I=ic; J=jc;
    % Coarse width ridges
    Hew_s=min(Vin(i,j  ),Vin(i+1,j));
    Hew_n=min(Vin(i,j+2),Vin(i+1,j+2));
    Hns_w=min(Uin(i  ,j),Uin(i  ,j+1));
    Hns_e=min(Uin(i+2,j),Uin(i+2,j+1));
    mskD=Deepest.NS>max(Hew_n,Hew_s); % Limit NS connectivity
    [I,J]=find(mskD & Hew_n>Hew_s); i=2*I-1; j=2*J-1; K=sub2ind(Ch,I,J);
      Vout(sub2ind(Nv,i,j+2))=max(Vout(sub2ind(Nv,i,j+2)),Deepest.NS(K)); % Block N
      Vout(sub2ind(Nv,i+1,j+2))=max(Vout(sub2ind(Nv,i+1,j+2)),Deepest.NS(K)); % Block N
      msk(sub2ind(Nh,i,j+1))=1;
    [I,J]=find(mskD & Hew_n<Hew_s); i=2*I-1; j=2*J-1; K=sub2ind(Ch,I,J);
      Vout(sub2ind(Nv,i,j))=max(Vout(sub2ind(Nv,i,j)),Deepest.NS(K)); % Block S
      Vout(sub2ind(Nv,i+1,j))=max(Vout(sub2ind(Nv,i+1,j)),Deepest.NS(K)); % Block S
      msk(sub2ind(Nh,i+1,j+1))=1;
    [I,J]=find(mskD & Hew_n==Hew_s); i=2*I-1; j=2*J-1; K=sub2ind(Ch,I,J);
      Vout(sub2ind(Nv,i,j+2))=max(Vout(sub2ind(Nv,i,j+2)),Deepest.NS(K)); % Block N
      Vout(sub2ind(Nv,i+1,j+2))=max(Vout(sub2ind(Nv,i+1,j+2)),Deepest.NS(K)); % Block N
      Vout(sub2ind(Nv,i,j))=max(Vout(sub2ind(Nv,i,j)),Deepest.NS(K)); % Block S
      Vout(sub2ind(Nv,i+1,j))=max(Vout(sub2ind(Nv,i+1,j)),Deepest.NS(K)); % Block S
      msk(sub2ind(Nh,i,j+1))=1;
      msk(sub2ind(Nh,i+1,j+1))=1;
    mskD=Deepest.EW>max(Hns_e,Hns_w); % Limit EW connectivity
    [I,J]=find(mskD & Hns_e>Hns_w); i=2*I-1; j=2*J-1; K=sub2ind(Ch,I,J);
      Uout(sub2ind(Nv,i+2,j))=max(Vout(sub2ind(Nv,i+2,j)),Deepest.EW(K)); % Block E
      Uout(sub2ind(Nv,i+2,j+1))=max(Vout(sub2ind(Nv,i+2,j+1)),Deepest.EW(K)); % Block E
      msk(sub2ind(Nh,i,j))=1;
    [I,J]=find(mskD & Hns_e<Hns_w); i=2*I-1; j=2*J-1; K=sub2ind(Ch,I,J);
      Uout(sub2ind(Nv,i,j))=max(Vout(sub2ind(Nv,i,j)),Deepest.EW(K)); % Block W
      Uout(sub2ind(Nv,i,j+1))=max(Vout(sub2ind(Nv,i,j+1)),Deepest.EW(K)); % Block W
      msk(sub2ind(Nh,i+1,j))=1;
    [I,J]=find(mskD & Hns_e==Hns_w); i=2*I-1; j=2*J-1; K=sub2ind(Ch,I,J);
      Uout(sub2ind(Nv,i+2,j))=max(Vout(sub2ind(Nv,i+2,j)),Deepest.EW(K)); % Block E
      Uout(sub2ind(Nv,i+2,j+1))=max(Vout(sub2ind(Nv,i+2,j+1)),Deepest.EW(K)); % Block E
      Uout(sub2ind(Nv,i,j))=max(Vout(sub2ind(Nv,i,j)),Deepest.EW(K)); % Block W
      Uout(sub2ind(Nv,i,j+1))=max(Vout(sub2ind(Nv,i,j+1)),Deepest.EW(K)); % Block W
      msk(sub2ind(Nh,i,j))=1;
      msk(sub2ind(Nh,i+1,j))=1;

   nncnt=sum(msk(:));
 end
 if ~exist('msg','var')
  locdisp(sprintf('  limitNSConnections: %i points (%5.2f%%) "filled"', ...
   nncnt,100*nncnt/(4*nxf*nyf) ))
 else
  if nncnt>0
   locdisp(sprintf('  limitNSConnections: %i points (%5.2f%%) "filled" (%s)', ...
    nncnt,100*nncnt/(4*nxf*nyf),msg))
  end
 end
% -----------------------------------------------------------------------------
function [Uout,Vout,msk] = limitDConnections(Deepest,Uin,Vin)
 [nxf,nyf]=size(Vin); nyf=nyf-1; msk=NaN*zeros(nxf,nyf);
 ic=1:nxf/2;jc=1:nyf/2;
 Uout=Uin;Vout=Vin;
 global lowMem
 if 1 % lowMem
  nncnt=0;
  for J=jc
   j=2*J-1;
   for I=ic 
    i=2*I-1;
    ncnt=0;
    % Coarse width ridges
    Hew_s=min(Vin(i,j  ),Vin(i+1,j));
    Hew_n=min(Vin(i,j+2),Vin(i+1,j+2));
    Hns_w=min(Uin(i  ,j),Uin(i  ,j+1));
    Hns_e=min(Uin(i+2,j),Uin(i+2,j+1));
    % Corner connections
    Hsw=max(Hew_s,Hns_w);
    Hnw=max(Hew_n,Hns_w);
    Hne=max(Hew_n,Hns_e);
    Hse=max(Hew_s,Hns_e);
   %Uout(i+1,j:j+1)=min(Uout(i+1,j:j+1),min([Hew_s Hew_n Hns_w Hns_e]))*NaN; Vout(i:i+1,j+1)=min(Vout(i:i+1,j+1),min([Hew_s Hew_n Hns_w Hns_e]))*NaN; % Throw away DEBUGging
    if Deepest.SW(I,J)>max(Hew_s,Hns_w) % Limit SW connectivity
     if Hns_w>Hew_s
      Uout(i,j:j+1)=max(Uout(i,j:j+1),Deepest.SW(I,J)); % Block W
     elseif Hns_w<Hew_s
      Vout(i:i+1,j)=max(Vout(i:i+1,j),Deepest.SW(I,J)); % Block S
     else
      Uout(i,j:j+1)=max(Uout(i,j:j+1),Deepest.SW(I,J)); % Block W
      Vout(i:i+1,j)=max(Vout(i:i+1,j),Deepest.SW(I,J)); % Block S
     end
     ncnt=ncnt+1;
     msk(i,j)=1;
    end
    if Deepest.NW(I,J)>max(Hew_n,Hns_w) % Limit NW connectivity
     if Hns_w>Hew_n
      Uout(i,j:j+1)=max(Uout(i,j:j+1),Deepest.NW(I,J)); % Block W
     elseif Hns_w<Hew_n
      Vout(i:i+1,j+2)=max(Vout(i:i+1,j+2),Deepest.NW(I,J)); % Block N
     else
      Uout(i,j:j+1)=max(Uout(i,j:j+1),Deepest.NW(I,J)); % Block W
      Vout(i:i+1,j+2)=max(Vout(i:i+1,j+2),Deepest.NW(I,J)); % Block N
     end
     ncnt=ncnt+1;
     msk(i,j+1)=1;
    end
    if Deepest.NE(I,J)>max(Hew_n,Hns_e) % Limit NE connectivity
     if Hns_e>Hew_n
      Uout(i+2,j:j+1)=max(Uout(i+2,j:j+1),Deepest.NE(I,J)); % Block E
     elseif Hns_e<Hew_n
      Vout(i:i+1,j+2)=max(Vout(i:i+1,j+2),Deepest.NE(I,J)); % Block N
     else
      Uout(i+2,j:j+1)=max(Uout(i+2,j:j+1),Deepest.NE(I,J)); % Block E
      Vout(i:i+1,j+2)=max(Vout(i:i+1,j+2),Deepest.NE(I,J)); % Block N
     end
     ncnt=ncnt+1;
     msk(i+1,j+1)=1;
    end
    if Deepest.SE(I,J)>max(Hew_s,Hns_e) % Limit SE connectivity
     if Hns_e>Hew_s
      Uout(i+2,j:j+1)=max(Uout(i+2,j:j+1),Deepest.SE(I,J)); % Block E
     elseif Hns_e<Hew_s
      Vout(i:i+1,j)=max(Vout(i:i+1,j),Deepest.SE(I,J)); % Block S
     else
      Uout(i+2,j:j+1)=max(Uout(i+2,j:j+1),Deepest.SE(I,J)); % Block E
      Vout(i:i+1,j)=max(Vout(i:i+1,j),Deepest.SE(I,J)); % Block S
     end
     ncnt=ncnt+1;
     msk(i+1,j)=1;
    end
    if ncnt>4
      locdisp('limitDConnections: This should NEVER happen!')
    end
    nncnt=nncnt+ncnt;
   end % I
  end % J
 else % !lowMem
  i=2*ic-1; j=2*jc-1;
  error('Not implemented yet')
 end
 if ~exist('msg','var')
  locdisp(sprintf('  limitDConnections: %i points (%5.2f%%) "filled"', ...
   nncnt,100*nncnt/(4*nxf*nyf) ))
 else
  if nncnt>0
   locdisp(sprintf('  limitDConnections: %i points (%5.2f%%) "filled" (%s)', ...
    nncnt,100*nncnt/(4*nxf*nyf),msg))
  end
 end
% -----------------------------------------------------------------------------
function [] = dump4Cells(label,Uin,Vin,Deepest)
global iDump jDump widthDump
if ~isempty(iDump) && size(Vin,1)==2*widthDump
 i=iDump*2-1;j=jDump*2-1;
 fprintf('%s: i=%i,j=%i\n',label,i,j)
 fprintf('\t      %5i       %5i\n',Vin(i:i+1,j+2))
 fprintf('\t%5i       %5i       %5i\n',Uin(i:i+2,j+1))
 fprintf('\t      %5i      %5i\n',Vin(i:i+1,j+1))
 fprintf('\t%5i       %5i       %5i\n',Uin(i:i+2,j))
 fprintf('\t      %5i       %5i\n',Vin(i:i+1,j))
 if ~isempty(Deepest)
  fprintf('%s: i=%i,j=%i\n',label,iDump,jDump)
  fprintf('NS=%i\t',Deepest.NS(iDump,jDump))
  fprintf('EW=%i\t',Deepest.EW(iDump,jDump))
  fprintf('SW=%i\t',Deepest.SW(iDump,jDump))
  fprintf('NW=%i\t',Deepest.NW(iDump,jDump))
  fprintf('NE=%i\t',Deepest.NE(iDump,jDump))
  fprintf('SE=%i\n',Deepest.SE(iDump,jDump))
 end
end
% -----------------------------------------------------------------------------
function [] = dump1Cell(label,Uin,Vin)
global iDump jDump widthDump
if ~isempty(iDump) && size(Vin,1)==widthDump
 i=iDump;j=jDump;
 fprintf('%s: i=%i,j=%i\n',label,i,j)
 fprintf('\t            %5i\n',Vin(i,j+1))
 fprintf('\t%5i                   %5i\n',Uin(i:i+1,j))
 fprintf('\t            %5i\n',Vin(i,j))
end
% -----------------------------------------------------------------------------
function [deepest] = findDeepestConnections(Uin,Vin)
deepest.NS=findDeepestNSconnection(Uin,Vin);
deepest.EW=findDeepestEWconnection(Uin,Vin);
deepest.SW=findDeepestSWconnection(Uin,Vin);
deepest.NW=findDeepestNWconnection(Uin,Vin);
deepest.NE=findDeepestNEconnection(Uin,Vin);
deepest.SE=findDeepestSEconnection(Uin,Vin);
% -----------------------------------------------------------------------------
function [] = checkDeepestConnections(OriginalDeepest,Uin,Vin,label)
Deepest=findDeepestConnections(Uin,Vin);
dump4Cells(label,Uin,Vin,Deepest)
if min(Deepest.NS(:)-OriginalDeepest.NS(:))<0 | min(Deepest.EW(:)-OriginalDeepest.EW(:))<0 ...
 | min(Deepest.NE(:)-OriginalDeepest.NE(:))<0 | min(Deepest.SW(:)-OriginalDeepest.SW(:))<0 ...
 | min(Deepest.NW(:)-OriginalDeepest.NW(:))<0 | min(Deepest.SE(:)-OriginalDeepest.SE(:))<0
 fprintf('checkDeepestConnections: Violation at %s\n',label)
 stats(Deepest.NS-OriginalDeepest.NS,'NS')
 stats(Deepest.EW-OriginalDeepest.EW,'EW')
 stats(Deepest.NE-OriginalDeepest.NE,'NE')
 stats(Deepest.SW-OriginalDeepest.SW,'SW')
 stats(Deepest.NW-OriginalDeepest.NW,'NW')
 stats(Deepest.SE-OriginalDeepest.SE,'SE')
 keyboard
 error('checkDeepestConnections: violation')
end
% -----------------------------------------------------------------------------
function [] = checkCoarseFineConnections(FineDeepest,Uin,Vin,label)
dump1Cell(label,Uin,Vin)
NS=max(Vin(:,1:end-1),Vin(:,2:end));
EW=max(Uin(1:end-1,:),Uin(2:end,:));
SW=max(Uin(1:end-1,:),Vin(:,1:end-1));
NW=max(Uin(1:end-1,:),Vin(:,2:end));
NE=max(Uin(2:end,:),Vin(:,2:end));
SE=max(Uin(2:end,:),Vin(:,1:end-1));
global iDump jDump widthDump
if ~isempty(iDump) && size(Vin,1)==widthDump
 fprintf('%s: i=%i,j=%i\n',label,iDump,jDump)
 fprintf('NS=%i\t',NS(iDump,jDump))
 fprintf('EW=%i\t',EW(iDump,jDump))
 fprintf('SW=%i\t',SW(iDump,jDump))
 fprintf('NW=%i\t',NW(iDump,jDump))
 fprintf('NE=%i\t',NE(iDump,jDump))
 fprintf('SE=%i\n',SE(iDump,jDump))
end
if min(NS(:)-FineDeepest.NS(:))<0 ...
 | min(EW(:)-FineDeepest.EW(:))<0 ...
 | min(SW(:)-FineDeepest.SW(:))<0 ...
 | min(NW(:)-FineDeepest.NW(:))<0 ...
 | min(NE(:)-FineDeepest.NE(:))<0 ...
 | min(SE(:)-FineDeepest.SE(:))<0
 fprintf('checkCoarseFineConnections: Violation at %s\n',label)
 stats(NS-FineDeepest.NS,'NS')
 stats(EW-FineDeepest.EW,'EW')
 stats(SW-FineDeepest.SW,'SW')
 stats(NW-FineDeepest.NW,'NW')
 stats(NE-FineDeepest.NE,'NE')
 stats(SE-FineDeepest.SE,'SE')
 keyboard
end
% -----------------------------------------------------------------------------
function [deepest] = findDeepestNSconnection(Uin,Vin)
% Finds the deepest pathway across the four fine grid cells
% of eight possible paths in a north-south direction.
% Deepest is defined as highest along any individual path
% but deepest (lowset) of any choice.
 nx=size(Vin,1);ny=size(Uin,2);
 deepest=zeros(nx/2,ny/2);
 global lowMem
 if lowMem
  for J=1:ny/2
   j=2*J-1;
   for I=1:nx/2
    i=2*I-1;
    % Four pathways originating south-west
    %  Two paths across W ridge
     H_n_ne=max(Uin(i+1,j+1),Vin(i+1,j+2)); % Path from NW cell to NE exit
     H_nw=min(H_n_ne,Vin(i,j+2)); % Path from NW cell to either N exit
     H_w=max(H_nw,Vin(i,j+1)); % Path through W ridge to either N exit
    %  Two paths across S ridge
     H_n_nw=max(Vin(i,j+2),Uin(i+1,j+1)); % Path from NE cell to NW exit
     H_ne=min(H_n_nw,Vin(i+1,j+2)); % Path from NE cell to either N exit
     H_s=max(H_ne,Vin(i+1,j+1)); % Path through E ridge to either N exit
    H_sw=max(Vin(i,j), min(H_s,H_w)); % Path across SW entry

    % Four pathways originating south-east
    %  Two paths across E ridge
     H_ne=min(Vin(i+1,j+2),H_n_nw); % Path from NE cell to either N exit
     H_e=max(H_ne,Vin(i+1,j+1)); % Path through E ridge to either N exit
    %  Two paths across S ridge
     H_nw=min(Vin(i,j+2),H_n_ne); % Path from NW cell to either N exit
     H_w=max(H_nw,Vin(i,j+1)); % Path through W to either N exit
     H_s=max(H_w,Uin(i+1,j)); % Path through S to either N exit
    H_se=max(Vin(i+1,j), min(H_s,H_e)); % Path across SE entry
    deepest(I,J)=min(H_sw,H_se);
   end
  end
 else % !lowMem
  ic=1:nx/2; jc=1:ny/2;
  i=2*ic-1; j=2*jc-1;
    % Four pathways originating south-west
    %  Two paths across W ridge
     H_n_ne=max(Uin(i+1,j+1),Vin(i+1,j+2)); % Path from NW cell to NE exit
     H_nw=min(H_n_ne,Vin(i,j+2)); % Path from NW cell to either N exit
     H_w=max(H_nw,Vin(i,j+1)); % Path through W ridge to either N exit
    %  Two paths across S ridge
     H_n_nw=max(Vin(i,j+2),Uin(i+1,j+1)); % Path from NE cell to NW exit
     H_ne=min(H_n_nw,Vin(i+1,j+2)); % Path from NE cell to either N exit
     H_s=max(H_ne,Vin(i+1,j+1)); % Path through E ridge to either N exit
    H_sw=max(Vin(i,j), min(H_s,H_w)); % Path across SW entry

    % Four pathways originating south-east
    %  Two paths across E ridge
     H_ne=min(Vin(i+1,j+2),H_n_nw); % Path from NE cell to either N exit
     H_e=max(H_ne,Vin(i+1,j+1)); % Path through E ridge to either N exit
    %  Two paths across S ridge
     H_nw=min(Vin(i,j+2),H_n_ne); % Path from NW cell to either N exit
     H_w=max(H_nw,Vin(i,j+1)); % Path through W to either N exit
     H_s=max(H_w,Uin(i+1,j)); % Path through S to either N exit
    H_se=max(Vin(i+1,j), min(H_s,H_e)); % Path across SE entry
    deepest=min(H_sw,H_se);
 end
% -----------------------------------------------------------------------------
function [deepest] = findDeepestEWconnection(Uin,Vin)
deepest=findDeepestNSconnection(Vin',Uin')';
% -----------------------------------------------------------------------------
function [deepest] = findDeepestSWconnection(Uin,Vin)
% Finds the deepest pathway across the four fine grid cells
% of six possible paths in a south-west direction.
% Deepest is defined as highest along any individual path
% but deepest (lowset) of any choice.
 nx=size(Vin,1);ny=size(Uin,2);
 deepest=zeros(nx/2,ny/2);
 global lowMem
 if lowMem
  for J=1:ny/2
   j=2*J-1;
   for I=1:nx/2
    i=2*I-1;
    % Common paths: across W ridge to NW exit
    Hblue=max(Uin(i,j+1),Vin(i,j+1)); % Path from SW cell to NW exit
    Hred=max(Uin(i,j+1),max(Uin(i+1,j+1),Vin(i+1,j+1))); % Path from SE cell to NW exit
    % Three pathways originating south-west
    Hzig=max(Uin(i+1,j),Hred); % Longest path from SW cell to NW ext
    H_w=min(Uin(i,j), min(Hzig,Hblue)); % Path from SW cell to either W exit
    H_sw=max(H_w,Vin(i,j)); % Path across SW edge to either W exit

    % Three pathways originating south-east
    Hq=max(Uin(i,j),Vin(i,j+1)); % Path from NW cell to SW exit
    Hqq=min(Hq, Uin(i,j+1)); % Path from NW cell to W exit
    Hredder=max(Hqq,max(Uin(i+1,j+1),Vin(i+1,j+1))); % Modified red path
    H_s=max(min(Hblue,Uin(i,j)),Uin(i+1,j)); % Path through S ridge to either W exit
    H_e=min(H_s,Hredder); % Path from SE cell to either W exit
    H_se=max(H_e,Vin(i+1,j)); % Path across SE edge to either W exit
    deepest(I,J)=min(H_sw,H_se);
   end
  end
 else % !lowMem
  ic=1:nx/2; jc=1:ny/2;
  i=2*ic-1; j=2*jc-1;
    % Common paths: across W ridge to NW exit
    Hblue=max(Uin(i,j+1),Vin(i,j+1)); % Path from SW cell to NW exit
    Hred=max(Uin(i,j+1),max(Uin(i+1,j+1),Vin(i+1,j+1))); % Path from SE cell to NW exit
    % Three pathways originating south-west
    Hzig=max(Uin(i+1,j),Hred); % Longest path from SW cell to NW ext
    H_w=min(Uin(i,j), min(Hzig,Hblue)); % Path from SW cell to either W exit
    H_sw=max(H_w,Vin(i,j)); % Path across SW edge to either W exit
    % Three pathways originating south-east
    Hq=max(Uin(i,j),Vin(i,j+1)); % Path from NW cell to SW exit
    Hqq=min(Hq, Uin(i,j+1)); % Path from NW cell to W exit
    Hredder=max(Hqq,max(Uin(i+1,j+1),Vin(i+1,j+1))); % Modified red path
    H_s=max(min(Hblue,Uin(i,j)),Uin(i+1,j)); % Path through S ridge to either W exit
    H_e=min(H_s,Hredder); % Path from SE cell to either W exit
    H_se=max(H_e,Vin(i+1,j)); % Path across SE edge to either W exit
    deepest=min(H_sw,H_se);
 end
% -----------------------------------------------------------------------------
function [deepest] = findDeepestNEconnection(Uin,Vin)
deepest=findDeepestSWconnection(Uin(end:-1:1,end:-1:1),Vin(end:-1:1,end:-1:1));
deepest=deepest(end:-1:1,end:-1:1);
% -----------------------------------------------------------------------------
function [deepest] = findDeepestNWconnection(Uin,Vin)
deepest=findDeepestSWconnection(Uin(:,end:-1:1),Vin(:,end:-1:1));
deepest=deepest(:,end:-1:1);
% -----------------------------------------------------------------------------
function [deepest] = findDeepestSEconnection(Uin,Vin)
deepest=findDeepestSWconnection(Uin(end:-1:1,:),Vin(end:-1:1,:));
deepest=deepest(end:-1:1,:);
% -----------------------------------------------------------------------------
function [m] = mmin(varargin)
m=varargin{1};
for n=2:nargin
 m=min(m,varargin{n});
end
% -----------------------------------------------------------------------------
function [m] = mmax(varargin)
m=varargin{1};
for n=2:nargin
 m=max(m,varargin{n});
end
% -----------------------------------------------------------------------------
function [] = warn(msk,msg)
n=length(find(~isnan(msk(:))));
if n>0
 msg=(sprintf('** WARNING **: %i events found for %s',n,msg));
 locdisp(msg)
%keyboard
%error(msg)
end
% -----------------------------------------------------------------------------
function [] = myfig(n)
sp1=3;sp2=2;
%subplot(sp1,sp2,n)
figure(n)
clf;%set(gca,'Position',[0 0 1 1]+0.05*[1 1 -2 -2])
% -----------------------------------------------------------------------------
function [] = plotcolorbar()
global doColorbar
if doColorbar
 colorbar('h')
end
% -----------------------------------------------------------------------------
function [] = showChangedCells(X,Y,M);
hold on
x=(X(1:end-1,1:end-1)+X(2:end,2:end))/2;
y=(Y(1:end-1,1:end-1)+Y(2:end,2:end))/2;
msize=4;
if min(size(M))<260
 msize=5;
end
plot( x(:).*M(:), y(:).*M(:), 'r.', 'MarkerSize', msize );
hold off
% -----------------------------------------------------------------------------
function [] = myplot(XG,YG,H,U,V,Hrng)
global cm
colormap(cm)
%plotthinwalls(XG,YG,H,U,V,Hrng);
if isempty(U)
 pcolor(XG,YG,H([1:end end],[1:end end]));shading flat;caxis(Hrng);
else
 plotraster(XG,YG,H,U,V,Hrng);caxis(Hrng)
end
drawnow
% -----------------------------------------------------------------------------
function [] = locdisp(msg)
global verbose
if verbose; disp(msg); end
% -----------------------------------------------------------------------------
function [] = checksill(label,Hf,Uf,Vf,G)
global xSill ySill
if ~isempty(xSill)
 figure(9)
 s=sill_depth(Hf.amin,Uf.emin,Vf.emin,G.X,G.Y,xSill,ySill);
 fprintf('%i sill depth after %s\n',s,label);
end
% -----------------------------------------------------------------------------
