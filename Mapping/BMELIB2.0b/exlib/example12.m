% example12                 - Example VI.6.1 : Merging Pdfs
%
% This shows the combination of a BME posterior pdf and
% an exogeneous a priori pdf, when :
% a) the exogeneous pdf is uniform
% b) the exogeneous pdf is normal
%
% For the normal case, the resulting pdf is normal again.
% This is proved by superimposing a Gaussian pdf with
% mean and variance obtained from theorical developments

clear
z=(-5:0.01:7)';
mG=0;vG=1;
pdfK=gausspdf(z,[mG vG]);

%%%%%% Case a)

pdfE=uniformpdf(z,[-1 6]);
pdfC=pdfprod(z,pdfK,pdfE);

figure;
set(gcf,'PaperUnits','centimeters');
set(gcf,'Units','centimeters');
set(gcf,'PaperPosition',[0.5 0.5 13 12]);
set(gcf,'Position',[0.5 0.5 13 12]);
subplot(2,1,1);
plot(z,pdfK,'--');hold on;
plot(z,pdfE,'-.');
plot(z,pdfC,'-');
axis([-4 7 0 0.5])
set(gca,'XTick',-4:1:7);
set(gca,'YTick',0:0.1:0.5);
set(gca,'FontSize',8);
xlabel('X values','FontSize',10);
ylabel('pdf','FontSize',10);
title('(a)','FontSize',10);

%%%%%% Case b)

mE=2;vE=3;
pdfE=gausspdf(z,[mE vE]);
pdfC=pdfprod(z,pdfK,pdfE);

subplot(2,1,2);
plot(z,pdfK,'--');hold on;
plot(z,pdfE,'-.');
plot(z,pdfC,'-');
plot(z,gausspdf(z,[(vG*mE+vE*mG)/(vG+vE),vG*vE/(vG+vE)]),'--r');
axis([-4 7 0 0.5])
set(gca,'XTick',-4:1:7);
set(gca,'YTick',0:0.1:0.5);
set(gca,'FontSize',8);
xlabel('X values','FontSize',10);
ylabel('pdf','FontSize',10);
title('(b)','FontSize',10);

print -deps example12

save example12

