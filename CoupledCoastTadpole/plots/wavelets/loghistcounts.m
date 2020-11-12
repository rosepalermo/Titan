function [N, edges, bin] = loghistcounts(X)
    logX = log(X);
    [N, logedges, bin] = histcounts(logX);
    edges = exp(logedges);
end