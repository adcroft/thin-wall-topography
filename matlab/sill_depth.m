function [sill,D] = sill_depth(H,Hx,Hy,x,y,xx,yy)
% sill = sill_depth(H,Hx,Hy,x,y,[x0 x1],[y0 y1])
%
% [He2,xe2,ye2]=read_etopo2([-35 0],[57 72]);
% xg=xe2-(xe2(2)-xe2(1))/2;xg(end+1)=xg(end)+(xe2(2)-xe2(1));
% yg=ye2-(ye2(2)-ye2(1))/2;yg(end+1)=yg(end)+(ye2(2)-ye2(1));
% [XG,YG]=ndgrid(xg,yg);
% x0=-33.5,y0=62.5,x1=-18.5,y1=69.5 % Denmark Strait
% s=sill_depth(He2,[],[],XG,YG,x0,y0,x1,y1)
% s=sill_depth(He2,[],[],XG,YG,[-33.5 -18.5],[62.5 69.5])
% load Results/iceland_0.125.mat
% [s8,D8]=sill_depth(G.H.eff,G.U.eff,G.V.eff,G.X,G.Y,x0,y0,x1,y1)
% s8a=sill_depth(G.H.ave,[],[],G.X,G.Y,x0,y0,x1,y1)
% load Results/iceland_0.25.mat
% [s4,D4]=sill_depth(G.H.eff,G.U.eff,G.V.eff,G.X,G.Y,x0,y0,x1,y1)
% s4a=sill_depth(G.H.ave,[],[],G.X,G.Y,x0,y0,x1,y1)
% load Results/iceland_0.5.mat
% [s2,D2]=sill_depth(G.H.eff,G.U.eff,G.V.eff,G.X,G.Y,x0,y0,x1,y1)
% s2a=sill_depth(G.H.ave,[],[],G.X,G.Y,x0,y0,x1,y1)
% load Results/iceland_1.mat
% [s1,D1]=sill_depth(G.H.eff,G.U.eff,G.V.eff,G.X,G.Y,x0,y0,x1,y1)
% s1a=sill_depth(G.H.ave,[],[],G.X,G.Y,x0,y0,x1,y1)
%
% x0=-14.5,y0=60.5,x1=-1.5,y1=62.5 % Faroe Bank Channel
% x0=-14.5,y0=63.5,x1=-8.5,y1=64.5 % Iceland-Faroe Bank
% x0=-75.5,y0=-60.5,x1=-21.5,y1=-57.5 % Drake passage
% x0=-4.5,y0=35.5,x1=-7.5,y1=35.5 % Gibraltar Strait

x0=xx(1); x1=xx(2);
y0=yy(1); y1=yy(2);

cax=caxis;

if isempty(Hx)
 [Hx,Hy]=gen_hxhy(H,[],[]);
 %pcolor(x,y,H([1:end end],[1:end end]));shading flat;colorbar
%pcol(x,y,H);shading flat;colorbar
else
%plotthinwalls(x,y,H,Hx,Hy);shading flat;colorbar
end
%cm_landwater(200,caxis);
hold on;plot([x0 x1],[y0 y1],'r.');hold off;drawnow

Hx([1 end],:)=Inf;
Hy(:,[1 end])=Inf;

[i0,j0]=find_ij(x,y,x0,y0);
[i1,j1]=find_ij(x,y,x1,y1);
D=NaN*H;Dx=NaN*Hx;Dy=NaN*Hy;
D(i0,j0)=H(i0,j0);

ii=min(i0,i1):max(i0,i1); jj=min(j0,j1):max(j0,j1);
[nx,ny]=size(H);
%i=1:nx;j=1:ny;
%im1=[1 1:(nx-1)];jm1=[1 1:(ny-1)];
%ip1=[2:nx nx];jp1=[2:ny ny];

%D = crawler(H,Hx,Hy,i0,j0);
%  pcolor(x,y,D([1:end end],[1:end end]));shading flat;colorbar
%  hold on;plot([x0 x1],[y0 y1],'r.');hold off
%sill=D(i1,j1);
%return

nits=max(abs(j0-j1),abs(i0-i1))*1.5; n=0; dd=1;
while n<nits | dd~=0
 n=n+1; oD=D;
 M=1+0*D;
 Mx=max(M([1:end end],:),M([1 1:end],:));
 Dx=min(D([1 1:end],:),D([1:end end],:));
 Dx=max(Dx,Hx).*Mx;
 D=min(D,min(Dx([1:end-1],:),Dx(2:end,:)));
 My=max(M(:,[1:end end]),M(:,[1 1:end]));
 Dy=min(D(:,[1 1:end]),D(:,[1:end end]));
 Dy=max(Dy,Hy).*My;
 D=min(D,min(Dy(:,[1:end-1]),Dy(:,2:end)));
 if mod(n,20)==0
 %pcolor(x,y,D([1:end end],[1:end end]));shading flat;colorbar
  pcol(x,y,D);shading flat;colorbar
  %cm_landwater(200,caxis);
  hold on;plot([x0 x1],[y0 y1],'r.');hold off
 %plot(1:61,  D(i0-30:i0+30,j0),'rx',...
 %    .5:61.5,Hx(i0-30:i0+31,j0),'k.',...
 %    .5:61.5,Dx(i0-30:i0+31,j0),'r.')
  drawnow
 end
 oD(isnan(oD))=0;
 dd=D(ii,jj)-oD(ii,jj); dd(isnan(dd))=0; dd=max(abs(dd(:)));
end
sill=D(i1,j1);
%D(D==sill)=NaN;
pcol(x,y,D);shading flat;colorbar
caxis([-1 1]*.5+sill);
hold on;plot([x0 x1],[y0 y1],'r.');hold off;drawnow

function [i0,j0] = find_ij(x,y,x0,y0)
xb=barj(bari(x));
yb=barj(bari(y));
r2=(xb-x0).^2+(yb-y0).^2;
minr2=min(r2(:));
[i0,j0]=find(r2==minr2,1,'first');
if length(i0)>1
 i0
 j0
 x(i0,j0)
 y(i0,j0)
 r2(i0,j0)
 error('Too many hits!')
end

function [bi] = bari(a)
bi=(a(1:end-1,:)+a(2:end,:))/2;
function [bj] = barj(a)
bj=(a(:,1:end-1,:)+a(:,2:end))/2;

function [D] = crawler(H,Hx,Hy,i0,j0)
[ni,nj]=size(H);
D=NaN*H;
D(i0,j0)=H(i0,j0);
[nx,ny]=size(H);
i=1:nx;j=1:ny;
im1=[1 1:(nx-1)];jm1=[1 1:(ny-1)];
ip1=[2:nx nx];jp1=[2:ny ny];

n=1;c=0;
while n
 c=c+1;n=0;

 % Look left
 [I,J]=find( ~isnan(D(im1,:)) & ( (~isnan(D) & D>max(D(im1,:),Hx(i,:))) | isnan(D) ) );
 D(sub2ind([ni nj],I,J))=max( D(sub2ind([ni nj],I-1,J)), Hx(sub2ind([ni+1 nj],I,J)) );
 n=n+length(I);
 % Look down
 [I,J]=find( ~isnan(D(:,jm1)) & ( (~isnan(D) & D>max(D(:,jm1),Hy(:,j)) | isnan(D)) ) );
 D(sub2ind([ni nj],I,J))=max( D(sub2ind([ni nj],I,J-1)), Hy(sub2ind([ni nj+1],I,J)) );
 n=n+length(I);
 % Look right
 [I,J]=find( ~isnan(D(ip1,:)) & ( (~isnan(D) & D>max(D(ip1,:),Hx(i+1,:)) | isnan(D)) ) );
 D(sub2ind([ni nj],I,J))=max( D(sub2ind([ni nj],I+1,J)), Hx(sub2ind([ni+1 nj],I+1,J)) );
 n=n+length(I);
 % Look up
 [I,J]=find( ~isnan(D(:,jp1)) & ( (~isnan(D) & D>max(D(:,jp1),Hy(:,j+1)) | isnan(D)) ) );
 D(sub2ind([ni nj],I,J))=max( D(sub2ind([ni nj],I,J+1)), Hy(sub2ind([ni nj+1],I,J+1)) );
 n=n+length(I);

 if mod(c,20)==0
  pcolor(D');shading flat;caxis([-3000 500]);colorbar;drawnow
 end
end

function [D] = pcrawler(H,Hx,Hy,i0,j0)
global D Dx Dy ni nj rec maxrec
[ni,nj]=size(H);
D=NaN*H;
D(i0,j0)=H(i0,j0);
Dx=Hx;Dy=Hy;

rec=0; maxrec=0;
crawl_right(i0+1,j0);
maxrec

function [] = crawl_right(i,j)
global D Dx Dy ni nj rec maxrec
rec=rec+1, maxrec=max(maxrec,rec);
if i>ni || D(i,j)<=max(D(i-1,j),Dx(i,j))
 rec=rec-1
 return
end
io=i;jo=j;
while i<=ni
 if isnan(D(i,j))
  D(i,j)=max(D(i-1,j),Dx(i,j));
  i=i+1;
 elseif D(i,j)>max(D(i-1,j),Dx(i,j))
  D(i,j)=max(D(i-1,j),Dx(i,j));
  crawl_right(i+1,j);
  crawl_up(i,j+1);
  crawl_down(i,j-1);
 else
  break
 end
end % while i<=ni
if mod(rec,10)==0
   pcolor(D([1:end end],[1:end end])');shading flat;colorbar;drawnow
end
crawl_up(io,jo+1);
crawl_down(io,jo-1);
rec=rec-1

function [] = crawl_up(i,j)
global D Dx Dy ni nj rec maxrec
rec=rec+1, maxrec=max(maxrec,rec);
if j>nj || D(i,j)<=max(D(i,j-1),Dy(i,j))
 rec=rec-1
 return
end
io=i;jo=j;
while j<=nj
 if isnan(D(i,j)) || D(i,j)>max(D(i,j-1),Dy(i,j))
  D(i,j)=max(D(i,j-1),Dy(i,j));
  j=j+1;
 elseif D(i,j)>max(D(i,j-1),Dy(i,j))
  D(i,j)=max(D(i,j-1),Dy(i,j));
  crawl_up(i,j+1);
  crawl_left(i-1,j);
  crawl_right(i+1,j);
 else
  break
 end
end % while j<=nj
if mod(rec,10)==0
   pcolor(D([1:end end],[1:end end])');shading flat;colorbar;drawnow
end
crawl_left(io-1,jo);
crawl_right(io+1,jo);
rec=rec-1

function [] = crawl_left(i,j)
global D Dx Dy ni nj rec maxrec
rec=rec+1, maxrec=max(maxrec,rec);
if i<1 || D(i,j)<=max(D(i+1,j),Dx(i+1,j))
 rec=rec-1
 return
end
io=i;jo=j;
while i>=1
 if isnan(D(i,j)) || D(i,j)>max(D(i+1,j),Dx(i+1,j))
  D(i,j)=max(D(i+1,j),Dx(i+1,j));
  i=i-1;
 else
  break
 end
end % while i>=1
if mod(rec,10)==0
   pcolor(D([1:end end],[1:end end])');shading flat;colorbar;drawnow
end
crawl_up(io,jo+1);
crawl_down(io,jo-1);
rec=rec-1

function [] = crawl_down(i,j)
global D Dx Dy ni nj rec maxrec
rec=rec+1, maxrec=max(maxrec,rec);
if j<1 || D(i,j)<=max(D(i,j+1),Dy(i,j+1))
 rec=rec-1
 return
end
io=i;jo=j;
while j>=1
 if isnan(D(i,j)) || D(i,j)>max(D(i,j+1),Dy(i,j+1))
  D(i,j)=max(D(i,j+1),Dy(i,j+1));
  j=j-1;
 else
  break
 end
end % while j>=1
if mod(rec,10)==0
   pcolor(D([1:end end],[1:end end])');shading flat;colorbar;drawnow
end
crawl_left(io-1,jo);
crawl_right(io+1,jo);
rec=rec-1
