function [x,y,theta,internalStateOut] = estRun(time, dt, internalStateIn, steeringAngle, pedalSpeed, measurement)
% In this function you implement your estimator. The function arguments
% are:
%  time: current time in [s]
%  dt: current time step [s]
%  internalStateIn: the estimator internal state, definition up to you.
%  steeringAngle: the steering angle of the bike, gamma, [rad]
%  pedalSpeed: the rotational speed of the pedal, omega, [rad/s]
%  measurement: the position measurement valid at the current time step
%
% Note: the measurement is a 2D vector, of x-y position measurement.
%  The measurement sensor may fail to return data, in which case the
%  measurement is given as NaN (not a number).
%
% The function has four outputs:
%  est_x: your current best estimate for the bicycle's x-position
%  est_y: your current best estimate for the bicycle's y-position
%  est_theta: your current best estimate for the bicycle's rotation theta
%  internalState: the estimator's internal state, in a format that can be understood by the next call to this function

% Example code only, you'll want to heavily modify this.
% this needs to correspond to your init function:

xm = internalStateIn.xm;
Pm = internalStateIn.Pm;

theta = xm(3);
omega = pedalSpeed;
gamma = steeringAngle;

B = 0.8;
r = 0.425;
R = [1.0881, 1.5315;
    1.5315, 2.9845];
Q = [0.01 0 0;
     0 0.01 0;
     0 0 0.1];
L = eye(3);
M = eye(2);

H = [1 0 -0.5*B*sin(theta);
     0 1 0.5*B*cos(theta)];
A = [1 0 -5*r*omega*sin(theta);
     0 1 5*r*omega*cos(theta);
     0 0 1];

xp = [5*r*omega*cos(theta);
      5*r*omega*sin(theta);
      5*r*omega/B*tan(gamma)]*dt + xm;
Pp = A*Pm*A' + L*Q*L';

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    K = Pp*H'*inv(H*Pp*H' + M*R*M');
    h = xp(1:2) + [0.5*B*cos(xp(3)); 0.5*B*sin(xp(3))];
    xm = xp + K*(measurement - h);
    Pm = (eye(3) - K*H)*Pp;
else
    Pm = Pp;
    xm = xp;
end

%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:

internalStateOut.xm = xm;
internalStateOut.Pm = Pm;
x = xm(1);
y = xm(2);
theta = xm(3);

end