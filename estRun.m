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

% Read in variables
particles = internalStateIn.particles;
omega = pedalSpeed;
gamma = steeringAngle;

% Variance of measurement noise as measured from calibration data file
R = [1.0881, 1.5315;
     1.5315, 2.9845];

% Propagate particles through the dynamics
p_particles = zeros(size(particles));
for i=1:size(particles,2)
    % Create process noise for each particle
    v1 = normrnd(0,0.1);
    v2 = normrnd(0,0.1);
    
    % First derivative of state vector
    particle_dot = [5*particles(4,i)*(omega + v1)*cos(particles(3,i));
                     5*particles(4,i)*(omega + v1)*sin(particles(3,i));
                     5*particles(4,i)*(omega + v1)*tan(gamma + v2)/particles(5,i);
                     0;
                     0];
    % Second derivative of state vector
    particle_ddot = [-(5*particles(4,i)*(omega + v1))^2*sin(particles(3,i))*tan(gamma + v2)/particles(5,i);
                     (5*particles(4,i)*(omega + v1))^2*cos(particles(3,i))*tan(gamma + v2)/particles(5,i);
                     0;
                     0;
                     0];
    % Second order dynamics approximation
    p_particles(:,i) = particles(:,i) + dt*particle_dot + dt^2/2*particle_ddot;
end

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % We have recieved a new measurement
    
    % Generate estimates of measurements we expect to recieve
    m_particles = zeros(2,size(p_particles,2));
    for i=1:size(p_particles,2)
        m_particles(:,i) = [p_particles(1,i) + 0.5*p_particles(5,i)*cos(p_particles(3,i)), ...
                            p_particles(2,i) + 0.5*p_particles(5,i)*sin(p_particles(3,i))];
    end
    
    % Check probability of getting each estimated measurment given the
    % actual measurement
    probs = mvnpdf(m_particles', measurement, R);
    if sum(probs) < 1e-9
        % Our estimator is too far from the measurement! Reset our
        % estimate to match the latest measurement
        disp('Estimator Reset!')
        for i=1:size(p_particles,2)
            p_particles(1:2,i) = [measurement(1) - 0.5*p_particles(5,i)*cos(p_particles(3,i));
                                  measurement(2) - 0.5*p_particles(5,i)*sin(p_particles(3,i))];
        end
        f_particles = p_particles;
        
        % Roughen particles after reset to prevent particle convergence
        sigma = zeros(size(f_particles,1),1);
        K = 0.01;
        for i=1:size(f_particles,1)
           % sigma_i = K*E_i*N^(-1/d)
           sigma(i) = K*abs(peak2peak(f_particles(i,:)))*size(f_particles,2)^(-1/size(f_particles,1));
        end
        delta = mvnrnd(zeros(5,1), diag(sigma), size(f_particles,2));
        f_particles = f_particles + delta';
    else
        % Estimate is reasonably valid, choose new particles
        alpha = 1/sum(probs);
        probs = alpha*probs;
        cumprobs = cumsum(probs);
        f_particles = zeros(size(p_particles));
        for i=1:size(p_particles,2)
            X = rand;
            idx = find(cumprobs >= X, 1);
            f_particles(:,i) = p_particles(:,idx);
        end
    end
else
    % We did not recieve a new measurement, so simply return the particles
    % from the prior update
    f_particles = p_particles;
end

%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:

% Use the median particle as estimate to reduce effect of outliers
internalStateOut.particles = f_particles;
x = median(f_particles(1,:));
y = median(f_particles(2,:));
theta = median(f_particles(3,:));

end