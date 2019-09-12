function initgaus = get_gaussian_boundary(domainsize, sigma, pv)

initgaus = zeros(domainsize);
initgaus(1,:) = pv;
initgaus(:,1) = pv;
initgaus(end,:) = pv;
initgaus(:,end) = pv;

x = linspace(-2,2,domainsize(1));
y = linspace(-2,2,domainsize(1));
[X, Y] = meshgrid(x,y);
gaussian = exp(-(X.^2 + Y.^2)/sigma^2);
initgaus = conv2(initgaus, gaussian, 'same');

% figure()
% imagesc(initgaus)
% title('gaus')
