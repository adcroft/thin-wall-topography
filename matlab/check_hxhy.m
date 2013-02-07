function []=check_hxhy(He,Hx,Hy)
% check_hxhy(H,Hu,Hv)
%
% Checks that edge values, Hu and Hv, are
% everywhere higher than the neighboring
% center values, H.
%
% e.g.
%   check_hxhy(H,Hu,Hv)

[nx,ny]=size(He);
%im=[nx 1:nx];jm=[ny 1:ny]; ip=[1:nx 1];jp=[1:ny 1]; % Periodic
im=[1 1:nx];jm=[1 1:ny]; ip=[1:nx nx];jp=[1:ny ny]; % Non-periodic

ii = find(Hx<max( He(im,:,:) , He(ip,:,:) ));
jj = find(Hy<max( He(:,jm,:) , He(:,jp,:) ));
if sum([ length(ii) length(jj) ])~=0
 error('H, Hu and Hv are inconsistent')
end
