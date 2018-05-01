n = 10;
scores = zeros(n-1,1);
for i=0:n
    scores(i+1) = runTrial(i,1);
    drawnow
end

mean(scores)