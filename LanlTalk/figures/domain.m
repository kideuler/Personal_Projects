%% flower hole

g = 0.9;
colorg = [g g g];

square = [0 0; 1 0; 1 1; 0 1; 0 0];
n = 1000;
theta = [0:2*pi/n:2*pi];
x = (0.25 + 0.1*sin(5*theta)).*cos(theta)/3 + 0.5;
y = (0.25 + 0.1*sin(5*theta)).*sin(theta)/3 + 0.5;

hold on
axis square
fill(square(:,1), square(:,2), colorg)
fill(x, y, 'white')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
text(0.15, 0.85, '\Omega', 'FontSize', 14)
text(0.6117,0.5317,'\leftarrow \Gamma_1', 'FontSize', 14)
text(1,0.3,'\Gamma_2 \rightarrow', 'FontSize', 14, 'HorizontalAlignment','right')
filename = 'domain_square_flower';
saveas(gcf, [aesem_paper_root, '/figures/',filename],'epsc');
system(['epstopdf ', aesem_paper_root, '/figures/',filename,'.eps --output ', aesem_paper_root, '/figures/',filename,'.pdf']);
pause(0.01);

%% ellipse hole

g = 0.9;
colorg = [g g g];

square = [0 0; 1 0; 1 1; 0 1; 0 0];
n = 1000;
theta = [0:2*pi/n:2*pi];
x = 0.2*cos(theta) + 0.5;
y = 0.1*sin(theta) + 0.5;

hold on
axis square
fill(square(:,1), square(:,2), colorg)
fill(x, y, 'white')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
text(0.15, 0.85, '\Omega', 'FontSize', 14)
text(0.7,0.5,'\leftarrow \Gamma_1', 'FontSize', 14)
text(1,0.3,'\Gamma_2 \rightarrow', 'FontSize', 14, 'HorizontalAlignment','right')
filename = 'domain_square_ellipse';
saveas(gcf, [aesem_paper_root, '/figures/',filename],'epsc');
system(['epstopdf ', aesem_paper_root, '/figures/',filename,'.eps --output ', aesem_paper_root, '/figures/',filename,'.pdf']);
pause(0.01);