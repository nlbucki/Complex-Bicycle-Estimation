function internalState = estInitialize
% Fill in whatever initialization you'd like here. This function
% generates the internal state of the estimator at time 0. You may do
% whatever you like here, but you must return something that is in the
% format as may be used by your run() function.
%

% we make the interal state a structure, with the first three elements the
% positions x, y; the angle theta; and our favourite colour.

% note that there is *absolutely no prescribed format* for this internal state.
% You can put in it whatever you like. Probably, you'll want to keep the position
% and angle, and probably you'll remove the color.

N = 10;
particles = zeros(5,N);
for i=1:N
    particles(1,i) = unifrnd(-20,20);
    particles(2,i) = unifrnd(-20,20);
    particles(3,i) = unifrnd(-pi,pi);
    particles(4,i) = unifrnd(0.425*0.95,0.425*1.05);
    particles(5,i) = unifrnd(0.8*0.9,0.8*1.1);
end

internalState.x = 0;
internalState.y = 0;
internalState.theta = pi/4;
internalState.particles = particles; 
end


