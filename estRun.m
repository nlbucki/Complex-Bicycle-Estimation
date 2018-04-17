function [x,y,theta,internalStateOut] = estRun(~, dt, internalStateIn, steeringAngle, pedalSpeed, measurement)
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
 
% parameters for UKF
n = 7;
kappa = 5;
alpha = 0.25;
lambda = alpha^2 * (n+kappa) - n;
beta = 3;
 
xim = [xm; zeros(4,1)];
Pxi = blkdiag(Pm,Q,R);
 
sigmas = zeros(n,2*n+1);
decomp = chol((n+lambda)*Pxi, 'lower');
for i = 1:n
    sigmas(:,i) = xim + decomp(:,i);
    sigmas(:,i+n) = xim - decomp(:,i);
end
sigmas(:,end) = xim;

sigma_xp = zeros(3,2*n+1);
for i = 1:2*n+1
    [~, retx] = ode45(@(t,x) dynamics(t,x, omega, gamma, B, r), [0 dt], sigmas(1:5,i));
    sigma_xp(:,i) = retx(end,1:3)';
%     s3 = sigmas(3,i) + (5*r*omega + sigmas(4,i))/B*tan(gamma + sigmas(5,i))*dt;
%     sigma_xp(:,i) = [sigmas(1,i) + B*(sin(s3)-sin(theta))/tan(gamma + sigmas(5,i));
%                   sigmas(2,i) - B*(cos(s3)-cos(theta))/tan(gamma + sigmas(5,i));
%                   s3];
%     sigma_xp(:,i) = [sigmas(1,i) + 5*r*omega*cos(theta)*dt;
%                      sigmas(2,i) + 5*r*omega*sin(theta)*dt;
%                      sigmas(3,i) + (5*r*omega + sigmas(4,i))/B*tan(gamma + sigmas(5,i))*dt];
end

xp = zeros(3,1);
for i=1:2*n
    xp = xp + (1/(2*(n+lambda)))*sigma_xp(:,i);
end
xp = xp + sigma_xp(:,2*n+1)*lambda/(n+lambda);

Pp = zeros(3);
for i = 1:2*n
    Pp = Pp + (1/(2*(n+lambda)))*(sigma_xp(:,i) - xp)*(sigma_xp(:,i) - xp)';
end
Pp = Pp + (lambda/(n+lambda) + (1-alpha^2 + beta))*(sigma_xp(:,2*n+1) - xp)*(sigma_xp(:,2*n+1) - xp)';

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    meas_sigmas = zeros(2,2*n+1);
    for i = 1:2*n+1
        meas_sigmas(1,i) = sigma_xp(1,i) + 0.5*B*cos(sigma_xp(3,i)) + sigmas(6,i);
        meas_sigmas(2,i) = sigma_xp(2,i) + 0.5*B*cos(sigma_xp(3,i)) + sigmas(7,i);
    end
    z_hat = zeros(2,1);
    for i=1:2*n
        z_hat = z_hat + (1/(2*(n+lambda)))*meas_sigmas(:,i);
    end
    z_hat = z_hat + meas_sigmas(:,2*n+1)*lambda/(n+lambda);
    
    Pzz = zeros(2);
    for i = 1:2*n
        Pzz = Pzz + (1/(2*(n+lambda)))*(meas_sigmas(:,i) - z_hat)*(meas_sigmas(:,i) - z_hat)';
    end
    Pzz = Pzz + (lambda/(n+lambda) + (1-alpha^2 + beta))*(meas_sigmas(:,2*n+1) - z_hat)*(meas_sigmas(:,2*n+1) - z_hat)';
    Pxz = zeros(3,2);
    for i = 1:2*n
        Pxz = Pxz + (1/(2*(n+lambda)))*(sigma_xp(:,i) - xp)*(meas_sigmas(:,i) - z_hat)';
    end
    Pxz = Pxz + (lambda/(n+lambda) + (1-alpha^2 + beta))*(sigma_xp(:,2*n+1) - xp)*(meas_sigmas(:,2*n+1) - z_hat)';
    K = Pxz/Pzz;
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

function dxdt = dynamics(~,x,w,g,B,r)
% ODE for dynamics

dxdt = [5*r*(w+x(4))*cos(x(3)); 5*r*(w+x(4))*sin(x(3)); 5*r*(w+x(4))*tan(g+x(5))/B; 0; 0];

end