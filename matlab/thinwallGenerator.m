function [GS] = thinwallGenerator(xg,yg,varargin)
% G = thinwallGenerator(xg,yg)
% G = thinwallGenerator(xg,yg,method)
% G = thinwallGenerator(xg,yg,Hin,xin,yin)
% G = thinwallGenerator(xg,yg,Hin,xin,yin,method)
%
% xg,yg  - coordinates of grid-cell edges (_EDGES_ not centers)
%          either vectors or matrices but NOT 3+ dimensions
% G      - Structure containing thinwall data on grid (xg,yg)
%
% Optional:
% Hin    - High-resolution data
% xin    - Longitude of Hin
% yin    - Latitude of Hin
% method - Method/options for 'fv6'
%
% Interpolates high-resolution gridded topography to a lower-resolution
% grid using a thin-wall connectivity-preserving algoirthm.
%
% Example:
%
% >> xg=-30:4:12; yg=0:4:12;
% >> G = thinwallGenerator(xg,yg)
% >> [XG,YG] = ndgrid(xg,yg);
% >> G = thinwallGenerator(XG,YG,'fv6')
%
% or
%
% >> [Hin,xin,yin] = read_etopo1(-180,180,-90,90);
% >> [Hin,xin,yin] = read_etopo1([-180 180],[-90 90]);
% >> [H2,Hx2,Hy2] = thinwallGenerator(xg,yg,Hin,xin,yin)

% Copyright 2008-2013 Alistair Adcroft, Princeton University.
%
% This file is part of the thin-wall-topography software suite.
%
% thin-wall-topography is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% thin-wall-topography is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with Foobar. If not, see http://www.gnu.org/licenses/.

debugmode=0;
meth='fv6';

lon_lo=min(xg(:)); lon_hi=max(xg(:)); lat_lo=min(yg(:)); lat_hi=max(yg(:));
disp( sprintf('thinwallGenerator: range of coordinates (%g,%g,%g,%g)', ...
         [lon_lo lon_hi lat_lo lat_hi]))

% Get the ETOPO5 data on a conveniently sized grid
if nargin == 2 | nargin == 3
 %disp('thinwallGenerator: extracting ETOPO5 data ...');
 %[Hin,xin,yin]=read_etopo5([lon_lo lon_hi],[lat_lo lat_hi]);
 %disp('thinwallGenerator: extracting ETOPO2 data ...');
 %[Hin,xin,yin]=read_etopo2([lon_lo lon_hi],[lat_lo lat_hi]);
  locdisp('thinwallGenerator: extracting Sandwell data ...');
  [Hin,xin,yin]=read_sandwell([lon_lo lon_hi],[lat_lo lat_hi]);
 %locdisp('thinwallGenerator: extracting Gebco data ...');
 %[Hin,xin,yin]=read_gebco_08([lon_lo lon_hi],[lat_lo lat_hi]);
 xin=xin(:);yin=yin(:)';
 if nargin == 3
  meth=varargin{1};
 end
elseif nargin == 5 | nargin == 6
 locdisp('thinwallGenerator: using hi-res data supplied by argument');
 Hin=varargin{1};
 xin=varargin{2};
 yin=varargin{3};
 if min(xg(:)) < min(xin(:)) | max(xg(:)) > max(xin(:))
  error('Range of xg is larger than xin!')
 end
 if min(yg(:)) < min(yin(:)) | max(yg(:)) > max(yin(:))
  error('Range of yg is larger than yin!')
 end
 i0=find( xin(1:end-1) <= lon_lo & lon_lo < xin(2:end) );
 if isempty(i0)
  i0=find( xin(1:end-1) <= lon_lo+360 & lon_lo+360 < xin(2:end) );
 end
 if isempty(i0)
  i0=find( xin(1:end-1) <= lon_lo-360 & lon_lo-360 < xin(2:end) );
 end
 i1=find( xin(1:end-1) < lon_hi & lon_hi <= xin(2:end) )+1;
 if isempty(i1)
  i1=find( xin(1:end-1) < lon_hi+360 & lon_hi+360 <= xin(2:end) )+1;
 end
 if isempty(i1)
  i1=find( xin(1:end-1) < lon_hi-360 & lon_hi-360 <= xin(2:end) )+1;
 end
 j0=find( yin(1:end-1) <= lat_lo & lat_lo < yin(2:end) );
 j1=find( yin(1:end-1) < lat_hi & lat_hi <= yin(2:end) )+1;
 if i0<i1
  Hin=Hin(i0:i1,j0:j1); xin=xin(i0:i1); yin=yin(j0:j1);
 else
  Hin=Hin([i0:end 1:i1],j0:j1); xin=xin([i0:end 1:i1]); yin=yin(j0:j1);
 end
 if nargin == 6
  meth=varargin{4};
 end
else
 error('thinwallGenerator: too many arguments!')
end

locdisp(sprintf('thinwallGenerator: hi-res dimensions are %i x %i',size(Hin)));
if debugmode ~= 0
 save ETOPO5.mat Hin xin yin
end

% Interpolate to fine grid
locdisp('thinwallGenerator: creating fine grid ...');
dxfmin=min( xin(2:end,1)-xin(1:end-1,1) );
dyfmin=min( yin(1,2:end)-yin(1,1:end-1) );
%dxfmin=dxfmin/2; dyfmin=dyfmin/2; % Over do the fine grid?
[xgf,ygf,xcf,ycf,nfinepow] = generate_finegrid(xg,yg,dxfmin,dyfmin);
locdisp(sprintf('thinwallGenerator: fine grid dimensions are %i x %i',size(xcf,1),size(ycf,2)));
if debugmode ~= 0
 save HFGRID xgf ygf xcf ycf nfinepow
end

locdisp('thinwallGenerator: interpolating hi-res to fine grid ...');
 Hfine=interp2(yin,xin,Hin,ycf,xcf,'nearest');
if debugmode ~= 0
 save HFINE.mat Hfine xcf ycf
end
clear Hin xin yin xcf ycf

% Recursively coarsen grid
locdisp('thinwallGenerator: recursively coarsening from fine grid to target ...');
if strcmp(meth,'fv')
 [H2,Hx2,Hy2,x2,y2]=xyrecur_fv(Hfine,[],[],xgf,ygf,nfinepow,varargin{5:end});
elseif strcmp(meth,'fv2')
 [H2,Hx2,Hy2,x2,y2]=xyrecur_fv2(Hfine,[],[],xgf,ygf,nfinepow,varargin{5:end});
elseif strcmp(strtok(meth),'fv3')
 [H2,Hx2,Hy2,x2,y2]=xyrecur_fv3(Hfine,[],[],xgf,ygf,nfinepow,meth);
elseif strcmp(strtok(meth),'fv4')
 [GS]=xyrecur_fv4(Hfine,xgf,ygf,nfinepow,meth);
elseif strcmp(strtok(meth),'fv5')
 [GS]=xyrecur_fv5(Hfine,xgf,ygf,nfinepow,meth);
elseif strcmp(strtok(meth),'fv6')
 [GS]=xyrecur_fv6(Hfine,xgf,ygf,nfinepow,meth);
else
 error('Unknown method')
end
clear Hfine xcfine ycfine

% =======================================
function [] = locdisp(str)
disp(str)
