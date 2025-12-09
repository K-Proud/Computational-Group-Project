% Number of gridpoints and domain size
N = 200; L = 10;
% position of kink
a=0;


% Spatial grid and grid spacing
x = linspace(-L,L,N);
dx = x(2)-x(1);

v=0.5;
gamma = 1/sqrt((1-v^2));

% Initial condition
p0 = tanh(gamma*(2*x-a)) ;




% Initialize the spatial parameters for simulation
T = 20; dt = 0.01; time = 0:dt:T;

% Preallocate the solution array
phi = zeros(N, length(time));
p1 = tanh(gamma*(2*x-a-v*dt));


phi(:, 1) = p0;
phi(:, 2) = p1; % Set initial condition

% Laplacian matrix
% Construct the Laplacian matrix using finite difference method
one_vec = ones(N,1);
Laplacian = spdiags([one_vec, -2*one_vec, one_vec],[1,0,-1],N,N);
% Set the boundary conditions
% Laplacian(1,1) = -2; Laplacian(N,N) = -2;
Laplacian = (1/dx^2)*Laplacian;

e = zeros(N,1);
e(1) = -tanh(L)/(dx)^2; e(N) = tanh(L)/(dx)^2;

% Time-stepping loop for the simulation
for n = 2:length(time)-1
    % Update the solution using a numerical scheme (e.g., explicit Euler method)
   phi(:,n+1) = 2*phi(:,n) - phi(:,n-1) + (dt^2)*(Laplacian*phi(:,n)-2*phi(:,n).*((phi(:,n).^2)-1) + e);
end


close all;
plot(phi(:,1))


figure;
saveToGif(x,phi,"animated.gif",[-L,L,-2,2],100)


function saveToGif(x,phi,name,axs,frames)
    plt = plot(x,phi(:,1));
    axis(axs)
    for i = 1:frames
        %exportgraphics(gcf,name,'Append',true);
        axis(axs)
        plt.YData = phi(:,round(i*end/frames)); % Update plot data for the current frame
        drawnow; % Update the figure window
        %pause(0.1); % Pause to control the frame rate
    end
end

