slvecx = [0, 0];
slvecy = [-1, 1];

plot(slvecx, slvecy)
hold on

slvecx = slvecx(2)-slvecx(1);
slvecy = slvecy(2)-slvecy(1);
% 
theta = linspace(-pi,pi);
xlosvec = cos(theta);
ylosvec = sin(theta);

% cosang = 1;
slvec_norm = sqrt(slvecx.^2 + slvecy.^2);
losvec_norm = sqrt(xlosvec.^2 + ylosvec.^2);

xlosvec_normed = xlosvec./losvec_norm;
ylosvec_normed = ylosvec./losvec_norm;

xslvec_normed = slvecx./slvec_norm;
yslvec_normed = slvecy./slvec_norm;


cosang = cos(deg2rad(90-rad2deg(acos(abs(xlosvec_normed*xslvec_normed + ylosvec_normed*yslvec_normed)))));
cosang = sin(acos(xlosvec_normed*xslvec_normed + ylosvec_normed*yslvec_normed));
cosang = sqrt(1-(xlosvec_normed*xslvec_normed + ylosvec_normed*yslvec_normed).^2);

for i=1:length(theta)
    plot([0, cosang(i)*xlosvec(i)], [0, cosang(i)*ylosvec(i)]);
end