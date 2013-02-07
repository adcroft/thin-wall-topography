function [hh] = pcol(varargin);
% pcol  Corrected pseudocolor plot (aka pcolor).
%
% Similar to pcolor except that pcol() doesn't drop the last column
% or row of data. It uses the shading flat method by default.
%
% pcol(data)
% pcol(x,y,data)
%
% where data(j,i) has netcdf-convention dimensions (y,x).
%
% If no coordiantes are provided, cell-centered registration is assumed.
% If the x,y coordinates are provided the coordinates must both be either
% vectors or matrices.
% If the dimensions of the coordinates match the data the cell-centered
% registration is assumed.
% If the dimensions of the coordinates is one larger than the data then
% grid registration is assumed.
%
% e.g. All the following plot 3x5 cells with centers at integer coordinates:
% xc=.5:4.5,yc=.5:2.5 
% pcol(xc,yc,rand(3,5))
% xg=0:5,yg=0:3
% pcol(xg,yg,rand(3,5))
% [XC,YC]=meshgrid(xc,yc);
% pcol(XC,YC,rand(3,5))
% [XG,YG]=meshgrid(xg,yg);
% pcol(XG,YG,rand(3,5))
%
% See also  PCOLOR, IMAGESC
%
% Written by Alistair Adcroft (Princeton Univ.), 2008

if nargin == 1
%hh=imagesc(data);
 data=varargin{1};
 hh=pcolor(.5:size(data,2)+.5,.5:size(data,1)+.5,data([1:end 1],[1:end 1]));
else
 x=varargin{1};
 y=varargin{2};
 data=varargin{3};
 if max(size(x))==length(x(:)) % 1-dimensional coordinates
  x=x(:)'; y=y(:)';
  if length(x)==size(data,2) % Cell-centered coordinates
   x=[3*x(1)-x(2) x(1:end-1)+x(2:end) 3*x(end)-x(end-1)]/2;
   y=[3*y(1)-y(2) y(1:end-1)+y(2:end) 3*y(end)-y(end-1)]/2;
  elseif length(x)==size(data,2)+1 % Grid-centered coordinates
  else
   error('Coordinate dimensions do not match the data dimensions')
  end
  hh=pcolor(x,y,data([1:end 1],[1:end 1]));
 else % 1-dimensional coordinates
  if sum(abs(size(x)-size(data)))==0 % Cell-centered coordinates
  %disp('Note from pcol: Cell-centered curvilinear coordinates are not fully implemented')
   x=[3*x(:,1)-x(:,2) x(:,1:end-1)+x(:,2:end) 3*x(:,end)-x(:,end-1)]/2;
   x=[3*x(1,:)-x(2,:);x(1:end-1,:)+x(2:end,:);3*x(end,:)-x(end-1,:)]/2;
   y=[3*y(1,:)-y(2,:);y(1:end-1,:)+y(2:end,:);3*y(end,:)-y(end-1,:)]/2;
   y=[3*y(:,1)-y(:,2) y(:,1:end-1)+y(:,2:end) 3*y(:,end)-y(:,end-1)]/2;
  elseif sum(abs(size(x)-size(data)))==2 % Cell-centered coordinates
  else
   error('Coordinate dimensions do not match the data dimensions')
  end
  hh=pcolor(x,y,data([1:end 1],[1:end 1]));
 end
end
shading flat;
