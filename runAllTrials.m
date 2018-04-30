n = 20;
scores = zeros(n-1,1);
for i=1:n
    scores(i+1) = runTrial(i,0);
end

mean(scores)