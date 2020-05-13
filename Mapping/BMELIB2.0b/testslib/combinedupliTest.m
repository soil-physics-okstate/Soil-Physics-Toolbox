% combineduplitest            - Test the combinedupli function

rand('state',2);
randn('state',2);
nh=100;
nsU=30;
nsG=20;
d=2;
ns=nsU+nsG;
ch=rand(nh,d);
zh=rand(nh,1);
cs=rand(ns,d);
zlow=-.1-rand(nsU,1);
zup=.1+rand(nsU,1);
[softpdftype,nlU,limiU,probdensU]=probaUniform(zlow,zup);
zm=rand(nsG,1)-0.5;
zv=0.1+2*rand(nsG,1);
[softpdftype,nlG,limiG,probdensG]=probaGaussian(zm,zv);
[softpdftype,nl,limi,probdens]=probacat(softpdftype,nlU,limiU,probdensU,...
    softpdftype,nlG,limiG,probdensG);
nhPerm=randperm(nh);
n=1;
for i=2:5
  ir=nhPerm(n:n+i);
  for j=2:i
    ch(ir(j),:)=ch(ir(1),:);
  end
  n=n+i+1;
end
nsPerm=randperm(ns);
n=1;
for i=2:5
  ir=nsPerm(n:n+i);
  for j=2:i
    cs(ir(j),:)=cs(ir(1),:);
  end
  n=n+i+1;
end
for i=2:5
  ir=nsPerm(n:n+i);
  for j=2:i
    cs(ir(j),:)=ch(i,:);
  end
  n=n+i+1;
end

[chC,zhC] = combinedupli(ch,zh);   
[chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,'ave'); 
[chC,zhC] = combinedupli(ch,zh,'min'); 
[chC,zhC] = combinedupli(ch,zh,'max'); 
[chC,zhC] = combinedupli(ch,zh,0.5); 
[chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,'hist'); 
[chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli([],[],'ave',cs,softpdftype,nl,limi,probdens); 
[chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli([],[],'hist',cs,softpdftype,nl,limi,probdens); 
[chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,'ave',cs,softpdftype,nl,limi,probdens); 
[chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,0.5,cs,softpdftype,nl,limi,probdens); 
[chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,'hist',cs,softpdftype,nl,limi,probdens); 
if exist('tinv')==2
  [chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,'studentT'); 
  [chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,'studentT',cs,softpdftype,nl,limi,probdens); 
end
disp('test complete');
