% make a movie from the figures
v = VideoWriter('whats_the_bug_now');
open(v);

for k = 1:92
    % Create a mat filename, and load it into a structure called matData.
%     matFileName = sprintf('4test%d.mat', k);
    figFileName = sprintf('4test%d.fig', k);
    if exist(figFileName, 'file')
        matData = load(figFileName);
    else
        fprintf('File %s does not exist.\n', figFileName);
    end
    axis([900 1100 450 650])
    axis square
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);

   