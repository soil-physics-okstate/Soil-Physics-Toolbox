function MODELSLIBtest(testtype);

% MODELSLIBtest             - Test for the modelslib directory (Jan 1, 2001)
%
% SYNTAX :
%
% MODELSLIBtest(testtype);
% 
% INPUT :
%
% testtype    scalar     selects what type of covariance to test
%                          1 - to test spatial covariance
%                          2 - for space/time separable covariance
%                          3 - for space/time non separable covariance
%                        the default value is 1

if nargin<1, 
  testtype=1;
  disp('Note: Syntax is MODELSLIBtest(testtype)');
  disp('Using default value testtype=1');
end;  


ccMat=[2 0.5;0.5 1];
nv=length(ccMat);
nd=2;
cNugget=0.2;
cc=1;
as=1;at=2;

switch testtype
case 1,  
  disp('Testing spatial covariance');
  covmodel={'nuggetC','gaussianC'};
  covparam={ {cNugget*ccMat,[]}, ...   
             {cc*ccMat, [as at]} };
case 2,  
  disp('Testing space/time separarable covariance');
  covmodel={'nuggetC/nuggetC','gaussianC/exponentialC'};
  covparam={ {cNugget*ccMat,[]}, ...   
             {cc*ccMat, [as at]} };
case 3,  
  disp('Testing space/time non separabpe covariance');
  covmodel={'nuggetCST','gaussianCST'};
  covparam={ {cNugget*ccMat,[]}, ...   
             {cc*ccMat, [as as/at]} };
otherwise, 
  disp('MODELSLIBtest: testtype=4 not allowed, check help for MODELSLIBtest.');
  return;
end;

if testtype==1
  c1=zeros(1,nd);
  n1=1;
  sg=(0:0.02:1.5)';
  nsg=length(sg);
  c2=[sg zeros(nsg,nd-1)];
  figure;
  for i=1:nv
    for j=1:nv
      Kst{i,j}=coord2K({c1,i},{c2,j*ones(nsg,1)},covmodel,covparam);
      subplot(nv,nv,j+(i-1)*nv);
      plot(sg,Kst{i,j});
      if j==1, ylabel('C(r)'); end;
      if i==nv, xlabel('r'); end;
    end;
  end;
  disp('test complete');
  return;
end;

c1=zeros(1,nd+1);
n1=1;
sg=(0:0.02:1.5)';
nsg=length(sg);
tg=(0:0.04:3)';
ntg=length(tg);
c2=[kron(ones(ntg,1),sg) zeros(nsg*ntg,nd-1) kron(tg,ones(nsg,1))];
n2=nsg*ntg;
for i=1:nv
  for j=1:i           
    Kst{i,j}=coord2K({c1,i},{c2,j*ones(n2,1)},covmodel,covparam);
    KstMat{i,j}=reshape(Kst{i,j},nsg,ntg);
  end;
end;
sgMat=reshape(c2(:,1),nsg,ntg);
tgMat=reshape(c2(:,end),nsg,ntg);

figure;
for i=1:nv
  for j=1:i 
    subplot(nv,nv,j+(i-1)*nv);
    plot(sgMat(:,1),KstMat{i,j}(:,1));
    title(sprintf('(%d,%d)',i,j));
    if j==1, ylabel('C(r,t=0)'); end;
    if i==nv, xlabel('r'); end;
  end;
end;

figure;
for i=1:nv
  for j=1:i 
    subplot(nv,nv,j+(i-1)*nv);
    plot(tgMat(1,:),KstMat{i,j}(1,:));
    title(sprintf('(%d,%d)',i,j));
    if j==1, ylabel('C(r=0,t)'); end;
    if i==nv, xlabel('t'); end;
  end;
end;

figure;
for i=1:nv
  for j=1:i 
    subplot(nv,nv,j+(i-1)*nv);
    pcolor(sgMat,tgMat,KstMat{i,j});
    colorbar;
    shading interp;
    title(sprintf('(%d,%d)',i,j));
    if j==1, ylabel('t'); end;
    if i==nv, xlabel('r'); end;
  end;
end;

disp('test complete');
