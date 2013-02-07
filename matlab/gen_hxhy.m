function [Hx,Hy,msk]=gen_hxhy(He,Hu,Hv)
% [Hu,Hv]=gen_hxhy(He,[],[]);
% [Hu,Hv]=gen_hxhy(He,Hu,Hv);
%
% Returns the depths at U and V points.
%
% e.g.
%      [Hx,Hy]=gen_hxhy(He,Hu,Hv);

[nx,ny]=size(He);
%im=[nx 1:nx];jm=[ny 1:ny]; ip=[1:nx 1];jp=[1:ny 1]; % Periodic
im=[1 1:nx];jm=[1 1:ny]; ip=[1:nx nx];jp=[1:ny ny]; % Non-periodic

%Hx=max( He([end 1:end],:,:) , He([1:end 1],:,:) );
Hx=max( He(im,:,:) , He(ip,:,:) );
if ~isempty(Hu)
 Hx=max(Hx,Hu);
 n=length( find(Hx~=Hu) );
 Mu=NaN*Hu;Mu( find(Hx~=Hu) )=1;
else
 n=prod(size(Hx));
 Mu=1+0*Hx;
end

%Hy=max( He(:,[end 1:end],:) , He(:,[1:end 1],:) );
Hy=max( He(:,jm,:) , He(:,jp,:) );
if ~isempty(Hv)
 Hy=max(Hy,Hv);
 n=length( find(Hy~=Hv) );
 Mv=NaN*Hv;Mv( find(Hy~=Hv) )=1;
else
 n=n+prod(size(Hy));
 Mv=1+0*Hy;
end
nxy=prod(size(Hx))+prod(size(Hy));
%disp(sprintf('  gen_hxhy: %i points changed (%4.1f%%)',n,n/nxy*100))

msk=max(Mu(1:end-1,:),Mu(2:end,:));
msk=max(msk,Mv(:,1:end-1));
msk=max(msk,Mv(:,2:end));
