function internalState = estInitialize

% Number of particles
N = 500;

% Initialize particles
particles = zeros(5,N);
for i=1:N
    % Start roughly at (0,0)
    particles(1,i) = unifrnd(-10,10);
    particles(2,i) = unifrnd(-10,10);
    
    % Start facing roughly northeast
    particles(3,i) = unifrnd(-pi/4, 3*pi/4);
    
    % Radius is +/- 5 percent of nominal value
    particles(4,i) = unifrnd(0.425*0.95,0.425*1.05);
    
    % Baseline is +/- 10 percent of nominal value
    particles(5,i) = unifrnd(0.8*0.9,0.8*1.1);
end

internalState.particles = particles; 

end


