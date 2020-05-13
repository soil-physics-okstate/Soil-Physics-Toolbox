function pairsindextest(testAccuracy,testSpeed)
% pairsindextest                 - Test the pairsindex function
%
% SYNTAX :
%
% pairsindextest(testAccuracy,testSpeed)
%
% INPUT :
%
% testAccuracy  scalar    1 to test accuracy of pairsindex, 0 o.w.
%                         default is 1
% testSpeed     scalar    1 to test speed of pairsindex, 0 o.w.
%                         default is 1

if nargin<1, testAccuracy=1; end;
if nargin<2, testSpeed=1; end
rand('state',0);

dist={'coord2dist'};

%%%%%%% Test accuracy by comparing to pairsplot

if testAccuracy
  n=100;
  disp(sprintf('n=%d',n));
  c1=rand(n,2);
  c2=1.1+rand(n+1,2);
  rLag=[2.1 2.2];
  rLagTol=[.001 .001];
  nr=length(rLag);
  
  for ir=1:nr,
    figure;
    [index{ir}]=pairsplot(c1,c2,[rLag(ir)-rLagTol(ir) rLag(ir)+rLagTol(ir)]);
    drawnow;
  end;
  
  [idxpairskron]=pairsindex(c1,c2,rLag,rLagTol,dist,'kron');
  for ir=1:nr
    if isequal(sortrows(idxpairskron{ir}),sortrows(index{ir})),
      disp(sprintf('test %d complete for kron',ir));
    else
      disp(sprintf('Error in test %d for kron',ir));
    end
  end
  
  [idxpairskronloop]=pairsindex(c1,c2,rLag,rLagTol,dist,'kronloop');
  for ir=1:nr
    if isequal(sortrows(idxpairskronloop{ir}),sortrows(index{ir})),
      disp(sprintf('test %d complete for kronloop',ir));
    else
      disp(sprintf('Error in test %d for kronloop',ir));
    end
  end
  
  [idxpairssuperblock]=pairsindex(c1,c2,rLag,rLagTol,dist,'superblock');
  for ir=1:nr
    if isequal(sortrows(idxpairssuperblock{ir}),sortrows(index{ir})),
      disp(sprintf('test %d complete for superblock',ir));
    else
      disp(sprintf('Error in test %d for superblock',ir));
    end
  end
  
end


%%%%%%% Test speed by comparing to pairsplot

if testSpeed  
  figure;
  hold on;
  clear n;
  %n{1}=[100:100:700 800:400:2800];
  %n{2}=[100:100:700 800:400:4800];
  %n{3}=[100:100:700 800:400:3200];
  n{1}=[100:100:700 800:400:2000];
  n{2}=[100:100:700 800:400:2000];
  n{3}=[100:100:700 800:400:2000];
  title('Computational work for different methods in pairsindex');
  ylabel('CPU time (seconds)');
  xlabel('number of points');
  hold on;
  lt={'o-b','*-r','+-m'};
  for m=1:3
    rand('state',1);
    switch m
    case 1, method='kron';
    case 2, method='kronloop';
    case 3, method='superblock';
    end
    for i=1:length(n{m})
      disp(sprintf('n=%d',n{m}(i)));
      c1=rand(n{m}(i),2);
      c2=rand(n{m}(i)+10,2);
      rLag=[0.5];
      rLagTol=[0.001];
      t0 = clock;
      [idxpairs]=pairsindex(c1,c2,rLag,rLagTol,dist,method);
      t{m}(i)=etime(clock,t0);
      disp(sprintf('method=%s, number of pairs=%d, time =%f',method,size(idxpairs{end},1),t{m}(i)));
      plot(n{m}(1:i),t{m}(1:i),lt{m});
      drawnow;
    end
  end
  for m=1:3
      h(m)=plot(n{m},t{m},lt{m});
  end
  legend(h,{'kron','kronloop','superblock'});
  save test1 n t
end

