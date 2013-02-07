function [xgf,ygf,xcf,ycf,nfinepow] = generate_finegrid(xg,yg,dxemin,dyemin)
% [xgf,ygf,xcf,ycf,nfinepow] = generate_finegrid(xg,yg,dxmin,dymin)
%
% Repeatedly divides grid xg,yg until every element is smaller than dxmin,dymin

if min(size(xg)) == 1 & min(size(yg)) == 1 % 1-D coordinates
 [xg,yg]=ndgrid(xg(:),yg(:));
elseif sum(abs(size(xg)-size(yg)))==0 % 2-D coordinates
else
 error('generate_finegrid: Not sure about the input coordinates!')
end

dxfmin=mdel(xg);
dyfmin=mdel(yg');
locdisp( sprintf(' generate_finegrid: initial grid range (%g,%g,%g,%g)', ...
  min(xg(:)),max(xg(:)),min(yg(:)), max(yg(:))) )
locdisp( sprintf(' generate_finegrid: initial grid resolution (%g,%g)', ...
  dxfmin*60,dyfmin*60 ) )
locdisp( sprintf(' generate_finegrid: target grid resolution (%g,%g)', ...
  dxemin*60,dyemin*60 ) )

fudgeFac = 0.9999;
if dxfmin > fudgeFac*dxemin | dyfmin > fudgeFac*dyemin
 [xgf,ygf]=divide_grid(xg,yg);
 nfinepow=1;
 for l=2:9
  locdisp( sprintf(' generate_finegrid: trial grid range (%g,%g,%g,%g)', ...
   min(xgf(:)),max(xgf(:)),min(ygf(:)), max(ygf(:))) )
  dxfmin=mdel(xgf);
  dyfmin=mdel(ygf');
  locdisp( sprintf(' generate_finegrid: trial grid resolution (%g,%g)', ...
   dxfmin*60,dyfmin*60 ) )
 %locdisp( sprintf(' generate_finegrid: refining grid (%i,%i) (%g,%g)', ...
 %     size(xgf,1),size(ygf,2),dxfmin,dyfmin) )
  if fudgeFac*dxfmin <= dxemin & fudgeFac*dyfmin <= dyemin
   break
  end
  [xgf,ygf]=divide_grid(xgf,ygf);
  nfinepow=nfinepow+1;
 end
end

if size(xgf,2) == 1 & size(ygf,1) == 1
 xcf=bar(xgf(1:end-1),xgf(2:end));
 ycf=bar(ygf(1:end-1),ygf(2:end));
else
 xcf=bar(bar(xgf(1:end-1,1:end-1,:),xgf(2:end,1:end-1,:)), ...
         bar(xgf(1:end-1,2:end,:),xgf(2:end,2:end,:)));
 ycf=bar(bar(ygf(1:end-1,1:end-1,:),ygf(2:end,1:end-1,:)), ...
         bar(ygf(1:end-1,2:end,:),ygf(2:end,2:end,:)));
end
locdisp( sprintf(' generate_finegrid: final grid (%i,%i) (%g,%g)', ...
      size(xgf,1),size(ygf,2), dxfmin*60, dyfmin*60 ))
locdisp(sprintf(' generate_finegrid: there are %i levels of refinement',nfinepow))


% =============================================================================

function [dxmin] = mdel(xg)

dx=abs( xg([2:end],:)-xg([1:end-1],:) );
dx=dx(:)';
dxmin=min(dx( find(dx ~= 0) ));
dx=abs( xg(:,[2:end])-xg(:,[1:end-1]) );
dx=dx(:)';
dxmin=min([dxmin dx( find(dx ~= 0) )]);

% =============================================================================

function [b] = bar(x1,x2)

xx=max(x1,x2);
xn=min(x1,x2);
d=xx-xn;
b=(xx+xn)/2;
%c=xn; c( find(c<0) )=c( find(c<0) )+360;
%c=(xx+c)/2;
%b( find(d>300) )=c( find(d>300) );
%b( find(b>180) )=b( find(b>180) )-360;

% =============================================================================

function [xgf,ygf] = divide_grid(xg,yg)

nxc=size(xg,1)-1;
nyc=size(yg,2)-1;
nxf=2*nxc;
nyf=2*nyc;

if size(xg,2) == 1
error
 xgf(1:2:nxf+1)=xg;
 xgf(2:2:nxf)=bar(xg(1:nxc),xg(2:nxc+1));
 xgf=xgf';
else
 xgf(1:2:nxf+1,1:2:nyf+1)=xg(:,:);
 xgf(2:2:nxf,1:2:nyf+1)=bar(xg(1:nxc,:),xg(2:nxc+1,:));
 xgf(1:2:nxf+1,2:2:nyf)=bar(xg(:,1:nyc),xg(:,2:nyc+1));
 xgf(2:2:nxf,2:2:nyf)=bar(bar(xg(1:nxc,1:nyc),xg(1:nxc,2:nyc+1)) ...
                         ,bar(xg(2:nxc+1,1:nyc),xg(2:nxc+1,2:nyc+1)));
end

if size(yg,1) == 1
 ygf(1:2:nyf+1)=yg;
 ygf(2:2:nyf)=bar(yg(1:nyc),yg(2:nyc+1));
else
 ygf(1:2:nxf+1,1:2:nyf+1)=yg(:,:);
 ygf(2:2:nxf,1:2:nyf+1)=bar(yg(1:nxc,:),yg(2:nxc+1,:));
 ygf(1:2:nxf+1,2:2:nyf)=bar(yg(:,1:nyc),yg(:,2:nyc+1));
 ygf(2:2:nxf,2:2:nyf)=bar(bar(yg(1:nxc,1:nyc),yg(1:nxc,2:nyc+1)) ...
                         ,bar(yg(2:nxc+1,1:nyc),yg(2:nxc+1,2:nyc+1)));
end

%xgf( find(xgf>180) )=xgf( find(xgf>180) )-360;

% ===================

function [] = locdisp(msg)
disp(msg)
