function score = runTrial(n)

experimentalRun = n;
fprintf(['Loading the data file #' num2str(experimentalRun) ' \n']);
filename = ['data/run_' num2str(experimentalRun,'%03d') '.csv'];
experimentalData = csvread(filename);

internalState = estInitialize();

numDataPoints = size(experimentalData,1);
estimatedPosition_x = zeros(numDataPoints,1);
estimatedPosition_y = zeros(numDataPoints,1);
estimatedAngle = zeros(numDataPoints,1);

dt = experimentalData(2,1) - experimentalData(1,1);
for k = 1:numDataPoints
    t = experimentalData(k,1);
    gamma = experimentalData(k,2);
    omega = experimentalData(k,3);
    measx = experimentalData(k,4);
    measy = experimentalData(k,5);
    
    %run the estimator:
    [x, y, theta, internalState] = estRun(t, dt, internalState, gamma, omega, [measx, measy]);

    %keep track:
    estimatedPosition_x(k) = x;
    estimatedPosition_y(k) = y;
    estimatedAngle(k) = theta;
end    

% make sure the angle is in [-pi,pi]
estimatedAngle = mod(estimatedAngle+pi,2*pi)- pi;

posErr_x = estimatedPosition_x - experimentalData(:,6);
posErr_y = estimatedPosition_y - experimentalData(:,7);
angErr   = mod(estimatedAngle - experimentalData(:,8) + pi, 2*pi) - pi;

fprintf('Final error: \n');
fprintf(['   pos x = ' num2str(posErr_x(end)) ' m \n']);
fprintf(['   pos y = ' num2str(posErr_y(end)) ' m \n']);
fprintf(['   angle = ' num2str(angErr(end)) ' rad \n']);
score = norm([posErr_x(end); posErr_y(end); angErr(end)],1);
% our scalar score
fprintf(['Score: ' num2str(score) ' \n'])

end