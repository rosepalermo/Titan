
% make a movie from the figures
v = VideoWriter('uniformt1v1');
open(v);

for k = 1:100
    % Create a mat filename, and load it into a structure called matData.
%     matFileName = sprintf('4test%d.mat', k);
    figFileName = sprintf('uniformt1v1%d.fig', k);
    if exist(figFileName, 'file')
%         matData = load(matFileName);
        h = open(figFileName)
        h.Position = [0,0,1501,750]
    else
        fprintf('File %s does not exist.\n', figFileName);
    end
%     axis([900 1100 450 650])
%     axis square
    frame = getframe(gcf);
    writeVideo(v,frame);
    close all
end

close(v);

   