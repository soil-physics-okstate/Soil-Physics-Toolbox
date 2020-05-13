% example11_2               - Example VI.5.5 : Part 2
%
% This is the second part of example 11 (please refer to
% the example11_1.m file for more explanations).
%
% A random sample set of the cumulative probability distribution
% function obtained using i) kriging with exact 50 simulated hard
% data, ii) BME with 50 interval data, and iii) indicator kriging
% with interval data are plotted and compared

clear;

load example11_1;
monotonicity=sum(min(diff(cdfInd'))>=0);

dz=z(2)-z(1);
rand('seed',3);
pos=randperm(400);

set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[1 1 15 19]);
set(gcf,'Position',[1 1 15 19]);
for i=1:12,
  subplot(4,3,i);
  plot(classe,cdfInd(pos(i),:),'b-','LineWidth',1.5);hold on;
  plot(z,cumsum(pdfKrig(:,pos(i)).*dz),'b--');
  plot(z,cumsum(pdfBME(:,pos(i)).*dz),'b-');
  axis([-3 3 0 1]);
  set(gca,'XTick',-3:1:3);
  set(gca,'YTick',0:0.25:1);
  set(gca,'FontSize',8');
  if mod(i-1,3)==0
    ylabel('cdf');
  end
  if i>9
    xlabel('X');
  end;
  box on;
  hold off;
end;

print -deps example11_2
save example11_2
