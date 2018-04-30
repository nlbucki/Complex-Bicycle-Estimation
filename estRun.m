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

internalStateOut = internalStateIn;
particles = internalStateIn.particles;
omega = pedalSpeed;
gamma = steeringAngle;
R = [1.0881, 1.5315;
     1.5315, 2.9845];

p_particles = zeros(size(particles));
for i=1:size(particles,2)
    v1 = normrnd(0,0.1);
    v2 = normrnd(0,0.1);
    p_particles(:,i) = particles(:,i) + ...
                     [5*particles(4,i)*(omega+v1)*cos(particles(3,i));
                     5*particles(4,i)*(omega+v1)*sin(particles(3,i));
                     5*particles(4,i)*(omega+v1)*tan(gamma+v2)/particles(5,i);
                     0; 0]*dt;
end

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    
    m_particles = zeros(2,size(p_particles,2));
    for i=1:size(p_particles,2)
        m_particles(:,i) = [p_particles(1,i) + 0.5*p_particles(5,i)*cos(p_particles(3,i)), ...
                            p_particles(2,i) + 0.5*p_particles(5,i)*sin(p_particles(3,i))];
    end
    probs = mvnpdf(m_particles', measurement, R);
    alpha = 1/sum(probs);
    probs = alpha*probs;
    cumprobs = cumsum(probs);
    f_particles = zeros(size(p_particles));
    for i=1:size(p_particles,2)
        X = rand;
        idx = find(cumprobs >= X, 1);
        f_particles(:,i) = p_particles(:,idx);
    end
else
    f_particles = p_particles;
end


%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:

internalStateOut.x = median(f_particles(1,:));
internalStateOut.y = median(f_particles(2,:));
internalStateOut.theta = median(f_particles(3,:));
internalStateOut.particles = f_particles;
x = internalStateOut.x;
y = internalStateOut.y;
theta = internalStateOut.theta;
end