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
Q = [0.5 0;
     0 0.1];

Pxi = blkdiag(Pm, Q, R);
xim = [xm; zeros(4,1)];
 
sigmas = zeros(7,14);
decomp = chol(7*Pxi, 'lower');
for i = 1:7
    sigmas(:,i) = xim + decomp(:,i);
    sigmas(:,i+7) = xim - decomp(:,i);
end

sigma_xp = zeros(3,14);
for i = 1:14
    s3 = sigmas(3,i) + (5*r*omega + sigmas(4,i))/B*tan(gamma + sigmas(5,i))*dt;
    sigma_xp(:,i) = [sigmas(1,i) + B*(sin(s3)-sin(theta))/tan(gamma + sigmas(5,i));
                  sigmas(2,i) - B*(cos(s3)-cos(theta))/tan(gamma + sigmas(5,i));
                  s3];
%     sigma_xp(:,i) = [sigmas(1,i) + 5*r*omega*cos(theta)*dt;
%                      sigmas(2,i) + 5*r*omega*sin(theta)*dt;
%                      sigmas(3,i) + (5*r*omega + sigmas(4,i))/B*tan(gamma + sigmas(5,i))*dt];
end

xp = mean(sigma_xp, 2);
Pp = zeros(3);
for i = 1:14
    Pp = Pp + 1/14*(sigma_xp(:,i) - xp)*(sigma_xp(:,i) - xp)';
end

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    meas_sigmas = zeros(2,14);
    for i = 1:14
        meas_sigmas(1,i) = sigma_xp(1,i) + 0.5*B*cos(sigma_xp(3,i)) + sigmas(6,i);
        meas_sigmas(2,i) = sigma_xp(2,i) + 0.5*B*cos(sigma_xp(3,i)) + sigmas(7,i);
    end
    z_hat = mean(meas_sigmas, 2);
    Pzz = zeros(2);
    for i = 1:14
        Pzz = Pzz + 1/14*(meas_sigmas(:,i) - z_hat)*(meas_sigmas(:,i) - z_hat)';
    end
    Pxz = zeros(3,2);
    for i = 1:14
        Pxz = Pxz + 1/14*(sigma_xp(:,i) - xp)*(meas_sigmas(:,i) - z_hat)';
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