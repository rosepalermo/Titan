function cells = get_intersecting_cells(start,final)
    
    if (start(1)-final(1) == 0)
        % perfectly vertical
        cells = zeros(2,abs(start(2)-final(2))+1);
        cells(1,:) = start(1);
        cells(2,:) = start(2):sign(final(2)-start(2)):final(2);
    elseif (start(2)-final(2) == 0)
        % perfectly horizontal
        cells = zeros(2,abs(start(1)-final(1))+1);
        cells(2,:) = start(2);
        cells(1,:) = start(1):sign(final(1)-start(1)):final(1);
    else
        % neither vertical nor horizontal
        
    end
    
end