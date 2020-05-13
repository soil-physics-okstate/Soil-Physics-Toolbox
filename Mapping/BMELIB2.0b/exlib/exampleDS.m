% exampleDS                 - Downscaling example

randn('state',2);

%
% mean and covariance of the local field X(s)
%
Zmean=2;
c01=4.0e-1;
ar1=300;
c02=0.2e-1;
ar2=20;
covmodel={'exponentialC','exponentialC'};
covparam={[c01 ar1],[c02 ar2]};

%
% fine grid of points for the local field X(s)
%
Ls1=300;
Ls2=300;
ns1=12;
ns2=12;
[S1 S2]=meshgrid([0:Ls1/ns1:Ls1],[0:Ls2/ns2:Ls2]);
s1=S1(:);
s2=S2(:);
[n1 n2]=size(S1);

%
% Sequential simulation of the X(s) at the fine grid
%
nhmax=15;
dmax=100;
[z]=simuseq([s1 s2],covmodel,covparam,nhmax,dmax);
z=z-mean(z)+Zmean;z(z<0)=-z(z<0);z=z-mean(z)+Zmean;
figure;
subplot(2,2,2);
pcolor(reshape(s1,n1,n2),reshape(s2,n1,n2),reshape(z,n1,n2));
shading flat;
cax=caxis;
colormap('gray');
colorbar;
xlabel('s_1 (Km)');
ylabel('s_2 (Km)');
title('(b)');

%
% coarse grid for Z(s), and calculation of Z(s) by aggragation
% of the values of X(s) at the fine grid
%
ns1u=3;
ns2u=3;
s1gu=0:Ls1/ns1u:Ls1+Ls1/ns1u;
s2gu=0:Ls2/ns2u:Ls2+Ls2/ns2u;
for i1=1:length(s1gu)-1,
  for i2=1:length(s2gu)-1,
    S1u(i1,i2)=s1gu(i1);
    S2u(i1,i2)=s2gu(i2);
    Zu(i1,i2)=mean( z( s1gu(i1) <= s1 & s1 < s1gu(i1+1) & s2gu(i2) <= s2 & s2 < s2gu(i2+1) ) );
  end;
end;
subplot(2,2,1);
pcolor(S1u,S2u,Zu);
shading flat;
caxis(cax);
colormap('gray');
colorbar;
xlabel('s_1 (Km)');
ylabel('s_2 (Km)');
title('(a)');

%
% upscaling covariance of X(s) to covariance of Z(s) by summation
%
[Vs1 Vs2]=meshgrid([0:Ls1/ns1:Ls1/ns1u-Ls1/ns1],[0:Ls2/ns2:Ls2/ns2u-Ls2/ns2]);
r=[0:20:300];
cr=coord2K(0,r(:),covmodel,covparam);
for ir=1:length(r)
  cruMat=coord2K([Vs1(:) Vs2(:)],[Vs1(:)+r(ir) Vs2(:)],covmodel,covparam);
  cru(ir)=mean(mean(cruMat));
end;
subplot(2,1,2);
hold on;
hcr=plot(r,cr,'k-');
hcru=plot(r,cru,'r.-');
legend([hcru hcr],'County covariance','Local (zip code) covariance');
xlabel('Spatial lag, r (Km)');
ylabel('Covariance');
title('(c)');

print -deps exampleDS 