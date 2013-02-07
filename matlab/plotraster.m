function [] = plotraster(varargin)
% plotraster(H,Hx,Hy)
% plotraster(H,Hx,Hy,[hmin hmax])
% plotraster(x,y,H,Hx,Hy)
% plotraster(x,y,H,Hx,Hy,[hmin hmax])

% Copyright 2008-2013 Alistair Adcroft, Princeton University.
%
% This file is part of the thin-wall-topography software suite.
%
% thin-wall-topography is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% thin-wall-topography is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with Foobar. If not, see http://www.gnu.org/licenses/.

Htop=Inf;
Hbot=-Inf;
Eps=0.05;
Eps=0.10;

if nargin==3 | nargin==4
 H=varargin{1};
 Hx=varargin{2};
 Hy=varargin{3};
 [nx,ny]=size(H);
 x=1:nx+1;
 y=1:ny+1;
 if nargin==4; Hbot=min(varargin{4}); Htop=max(varargin{4}); end
elseif nargin==5 | nargin==6
 x=varargin{1};
 y=varargin{2};
 H=varargin{3};
 Hx=varargin{4};
 Hy=varargin{5};
 [nx,ny]=size(H);
 if nargin==6; Hbot=min(varargin{6}); Htop=max(varargin{6}); end
else
 error('Wrong number of arguments!')
end

if isempty(Hx)
 [Hx,Hy]=gen_hxhy(H,[],[]);
end

coordsAre2D=prod(prod(abs([size(x)-size(H); size(y)-size(H)])))==1;
if coordsAre2D
 xf=x;
elseif length(x)==size(H,1)
 xf=x;xf(end+1)=2*x(end)-x(end-1);
elseif length(x)==size(H,1)+1
 xf=x(:)';
else
 error('Oops X?')
end
if coordsAre2D
 yf=y;
elseif length(y)==size(H,2)
 yf=y;yf(end+1)=2*y(end)-y(end-1);
elseif length(y)==size(H,2)+1
 yf=y(:)';
else
 error('Oops Y?')
end

jh=isnan(H);H=min(Htop,max(Hbot,H));
j=isnan(Hx);Hx=min(Htop,max(Hbot,Hx));Hx(j)=NaN;
j=isnan(Hy);Hy=min(Htop,max(Hbot,Hy));Hy(j)=NaN;
H(jh)=NaN;
if Htop == 0 | Hbot == 0
 H=sq(H);
 Hx=sq(Hx);
 Hy=sq(Hy);
end

% Create new coordinates
if sum(abs(size(xf)-size(yf)))~=0
 error('Need 2D coords for now')
end
[nxg,nyg]=size(xf);
i=[1-Eps:nxg-Eps; 1+Eps:nxg+Eps]; i=i(:); i(1)=1; i(end)=nxg;
j=[1-Eps:nyg-Eps; 1+Eps:nyg+Eps]; j=j(:); j(1)=1; j(end)=nyg;
[I,J]=ndgrid(i,j);
XG=interp2(1:nyg,1:nxg,xf,J,I);
YG=interp2(1:nyg,1:nxg,yf,J,I);

Hras=NaN*zeros(2*nxg-1,2*nyg-1);
Hras(2:2:end-1,2:2:end-1)=H;
Hras(1:2:end,2:2:end-1)=Hx;
Hras(2:2:end-1,1:2:end)=Hy;
%Hras(3:2:end-2,3:2:end-2)=1*(...
%               (Hx(2:end-1,1:end-1)+Hx(2:end-1,2:end)) ...
%              +(Hy(1:end-1,2:end-1)+Hy(2:end,2:end-1)) )/4;
Hras(3:2:end-2,3:2:end-2)=max(...
               max(Hx(2:end-1,1:end-1),Hx(2:end-1,2:end)) ...
              ,max(Hy(1:end-1,2:end-1),Hy(2:end,2:end-1)) );
pcolor(XG,YG,Hras([1:end end],[1:end end]));shading flat

if isinf(Hbot)
 hmin=min([min(H(:)) min(Hx(:)) min(Hy(:))]);
else
 hmin=Hbot;
end
if isinf(Htop)
 hmax=max([max(H(:)) max(Hx(:)) max(Hy(:))]);
else
 hmax=Htop;
end
caxis([hmin hmax]);
if isinf(Hbot) & isinf(Htop)
 cm_landwater(length(colormap),[hmin hmax]);
end
drawnow
