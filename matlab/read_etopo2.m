function [h,x,y,ii,jj] = read_etopo2(xi,yi);
% [h,x,y] = read_etopo2(xi,yi);
%
% Read the 2 minute ETOPO2 dataset to a regular grid that encompasses the
% end values of (XI,YI)

x0=xi(1); x1=xi(end); y0=yi(1); y1=yi(end);
while x1<=x0 % translate relative bounds to +ve style
 x1=x1+360;
end

DatasetDir='/home/aja/projects/thinwalls/matlab/datasets';
nc=netcdf([DatasetDir '/ETOPO2v2.cdf'],'read');

Xe0=nc{'XF'}(1); Xe1=nc{'XF'}(end); Xe=nc{'X'}(:);
while x0>Xe1 % shift Xe to be in phase with xi coordinates
 Xe0=Xe0+360; Xe1=Xe1+360; Xe=Xe+360;
end
while x0<Xe0 % shift Xe to be in phase with xi coordinates
 Xe0=Xe0-360; Xe1=Xe1-360; Xe=Xe-360;
end

if x1<=Xe1 % xi fits with in range Xe
 i1=find( Xe([2:end end])>=x0 & Xe([1 1:end-1])<=x1 );
 i2=[];xshift2=0;
else % xi wraps off of the right
 i1=find( Xe([2:end end])>=x0 );
 i2=find( Xe([1 1:end-1])<=x1-360 ); xshift2=360;
end
%[min(i1) max(i1) min(i2) max(i2)]
x=[Xe(i1)' Xe(i2)'+xshift2]';
Ye=nc{'Y'}(:);
jj=find( Ye([2:end end])>=y0 & Ye([1 1:end-1])<=y1 );
y=Ye(jj)';
ii=[i1' i2']';
%[min(i1) max(i1) min(i2) max(i2) min(jj) max(jj)]
h=nc{'H'}(jj,i1)'; % (i,j) indexing
%h=nc{'bath'}(jj,i1); % (j,i) indexing
if ~isempty(i2)
 h(length(i1)+(1:length(i2)),:)=nc{'H'}(jj,i2)'; % (i,j) indexing
 %h(:,length(i1)+(1:length(i2)))=nc{'H'}(jj,i2); % (j,i) indexing
end

close(nc)
