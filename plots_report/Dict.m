close all

s1 = centers(:,10);
s2 = centers(:,510);
s3 = centers(:,1010);
s4 = centers(:,1510);

figure
subplot(4,1,1); plot(s1);
title('Exemple de atomi din dictionar');
subplot(4,1,2); plot(s2);
subplot(4,1,3); plot(s3);
subplot(4,1,4); plot(s4);

% Save as fig
savefig('Dict.fig')

% Save as pdf
ps = get(gcf, 'Position');
ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
print(gcf, '-dpdf', 'Dict.pdf');
