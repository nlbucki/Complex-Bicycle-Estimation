function score = runTrial(n,show_plots)
tic;
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

toc;

%%
%==========================================================================
% Make some plots:
%==========================================================================
% Feel free to add additional plots, if you like.
fprintf('Generating plots \n')

if show_plots
    %% Figure 1
    figure;
    hold on;
    plot(experimentalData(:,4), experimentalData(:,5), 'rx');
    plot(estimatedPosition_x, estimatedPosition_y, 'b-');
    plot(experimentalData(:,6), experimentalData(:,7), 'k:.');
    hold off

    xlabel('x-position [m]');
    ylabel('y-position [m]');
    legend('meas','est','true');

    %% Figure 2
    figure;
    subplot(5,1,1);
    hold on;
    plot(experimentalData(:,1), experimentalData(:,6), 'k:.');
    plot(experimentalData(:,1), experimentalData(:,4), 'rx');
    plot(experimentalData(:,1), estimatedPosition_x, 'b-');
    ylabel('Position x [m]');
    legend('truth','meas','est');
    hold off;
    subplot(5,1,2);
    hold on;
    plot(experimentalData(:,1), experimentalData(:,7), 'k:.');
    plot(experimentalData(:,1), experimentalData(:,5), 'rx');
    plot(experimentalData(:,1), estimatedPosition_y, 'b-');
    ylabel('Position y [m]');
    hold off;
    subplot(5,1,3);
    hold on;
    plot(experimentalData(:,1), experimentalData(:,8), 'k:.');
    plot(experimentalData(:,1), estimatedAngle, 'b-');
    ylabel('Angle theta [rad]');
    hold off;
    subplot(5,1,4);
    plot(experimentalData(:,1), experimentalData(:,2), 'g-');
    ylabel('Steering angle gamma [rad]');
    subplot(5,1,5);
    plot(experimentalData(:,1), experimentalData(:,3), 'g-');
    ylabel('Pedal speed omega [rad/s]');
    xlabel('Time [s]');
end
end