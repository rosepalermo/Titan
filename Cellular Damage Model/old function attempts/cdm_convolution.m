% setup the grid
x = linspace(-2,2,100);
y = linspace(-2,2,100);
[X,Y] = meshgrid(x,y);

% define the circular lake
% 0 is land
% -1 is lake
Z = double(X.^2 + Y.^2 > 1)-1;

% create the filter
filter_size = 11; % must be an odd number for the convolution
std_dev = .01; % standard deviation for the filter
filter_x = linspace(-3*std_dev,3*std_dev,filter_size);
filter_y = linspace(-3*std_dev,3*std_dev,filter_size);
[filter_X,filter_Y] = meshgrid(filter_x,filter_y);

filter = normpdf(sqrt(filter_X.^2 + filter_Y.^2),0,std_dev); % gaussian
filter = filter/sum(filter(:)); % normalized gaussian

% plot the filter
figure(1); clf
imagesc(filter_X, filter_Y, filter)
shading flat

% % simulate the erosion
figure(2);
threshold = -0.5;
for i=1:30
    % save the last Z
    Z_last = Z;
    
    % perform the convolution
    Z = conv2(Z,filter,'same');
    
    % erode the lake
    Z(Z < threshold) = -1; 
    
    % make sure the lake remains lake, this is equivalent to a Dirichlet boundary condition on the lake
    Z(Z_last < threshold) = -1; 
    
    % plot the damage
    figure(2)
    subplot(2,1,1);
    imagesc(X,Y,Z);
    shading flat
    axis square
    colorbar
    
    % plot the binary lake
    subplot(2,1,2);
    lake = double(Z > threshold);
    imagesc(X,Y,lake);
    shading flat
    axis square
    colorbar
    pause(0.03)
end
