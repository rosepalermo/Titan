% make a movie from the figures
v = VideoWriter('4fetch5_2018_zoom1');
open(v);

for k = 1:92
    % Create a mat filename, and load it into a structure called matData.
    matFileName = sprintf('4test%d.mat', k);
    if exist(matFileName, 'file')
        matData = load(matFileName);
    else
        fprintf('File %s does not exist.\n', matFileName);
    end
    axis([900 1100 450 650])
    axis square
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);

   