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
Q = [0.1 0 0;
     0 0.1 0;
     0 0 0.1];

sigmas = zeros(3,6);
decomp = chol(3*Pm, 'lower');
for i = 1:3
    sigmas(:,i) = xm + decomp(:,i);
    sigmas(:,i+3) = xm - decomp(:,i);
end

for i = 1:6
    s3 = sigmas(3,i) + 5*r*omega/B*tan(gamma)*dt;
    sigmas(:,i) = [sigmas(1,i) + B*(sin(s3)-sin(theta))/tan(gamma);
                  sigmas(2,i) - B*(cos(s3)-cos(theta))/tan(gamma);
                  s3];
end

xp = mean(sigmas, 2);
Pp = zeros(3);
for i = 1:6
    Pp = Pp + 1/6*(sigmas(:,i) - xp)*(sigmas(:,i) - xp)';
end
Pp = Pp + Q;

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    meas_sigmas = zeros(2,6);
    for i = 1:6
        meas_sigmas(1,i) = sigmas(1,i) + 0.5*B*cos(sigmas(3,i));
        meas_sigmas(2,i) = sigmas(2,i) + 0.5*B*cos(sigmas(3,i));
    end
    z_hat = mean(meas_sigmas, 2);
    Pzz = zeros(2);
    for i = 1:6
        Pzz = Pzz + 1/6*(meas_sigmas(:,i) - z_hat)*(meas_sigmas(:,i) - z_hat)';
    end
    Pzz = Pzz + R;
    Pxz = zeros(3,2);
    for i = 1:6
        Pxz = Pxz + 1/6*(sigmas(:,i) - xp)*(meas_sigmas(:,i) - z_hat)';
    end
    K = Pxz*inv(Pzz);
    xm = xp + K*(measurement' - z_hat);
    Pm = Pp - K*Pzz*K';
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