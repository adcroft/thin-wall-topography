function [h,x,y,ii,jj] = read_etopo5(xi,yi);
% [h,x,y] = read_etopo5(xi,yi);
%
% Read the 5 minute ETOPO5 dataset to a regular grid that encompasses the
% end values of (XI,YI)

x0=xi(1); x1=xi(end); y0=yi(1); y1=yi(end);

DatasetDir='/home/aja/projects/thinwalls/matlab/datasets/';
FileName='etopo5_new.cdf';
nc=netcdf([DatasetDir FileName],'read');

Xe=nc{'X'}(:);
if x0>=0
 if x1>x0
  i1=find( Xe>=x0 & Xe<=x1 );xshift1=0;
  i2=[];xshift2=0;
 else
  i1=find( Xe>=x0 ); xshift1=-360;
  i2=find( Xe<=x1 ); xshift2=0;
 end
else
 if x1<0
  i1=find( Xe>=x0+360 & Xe<=x1+360 );xshift1=-360;
  i2=[];xshift2=-360;
 else
  i1=find( Xe>=x0+360 );xshift1=-360;
  i2=find( Xe<=x1 );xshift2=0;
 end
end
%[min(i1) max(i1) min(i2) max(i2)]
x=[Xe(i1)'+xshift1 Xe(i2)'+xshift2]';
Ye=nc{'Y'}(:);
jj=find( Ye>=y0 & Ye<=y1 );
y=Ye(jj)';
ii=[i1' i2']';
[min(i1) max(i1) min(i2) max(i2) min(jj) max(jj)]
h=nc{'bath'}(jj,i1)'; % (i,j) indexing
%h=nc{'bath'}(jj,i1); % (j,i) indexing
if ~isempty(i2)
 h(length(i1)+(1:length(i2)),:)=nc{'bath'}(jj,i2)'; % (i,j) indexing
 %h(:,length(i1)+(1:length(i2)))=nc{'bath'}(jj,i2); % (j,i) indexing
end
