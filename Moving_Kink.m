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
p0 = tanh(gamma*(x-a)) ;


% Initialize the spatial parameters for simulation
T = 15; dt = 0.01; time = 0:dt:T;

% Preallocate the solution array
phi = zeros(N, length(time));
p1 = tanh(gamma*(x-a-v*dt));


phi(:, 1) = p0;
phi(:, 2) = p1; % Set initial condition

% Laplacian matrix
% Construct the Laplacian matrix using finite difference method
one_vec = ones(N,1);
Laplacian = spdiags([one_vec, -2*one_vec, one_vec],[1,0,-1],N,N);
% Set the boundary conditions
Laplacian = (1/dx^2)*Laplacian;

e = zeros(N,1);
e(1) = -tanh(L)/(dx)^2; e(N) = tanh(L)/(dx)^2;

% Time-stepping loop for the simulation
for n = 2:length(time)-1
    % Update the solution using a numerical scheme (e.g., explicit Euler method)
   phi(:,n+1) = 2*phi(:,n) - phi(:,n-1) + (dt^2)*(Laplacian*phi(:,n)-2*phi(:,n).*((phi(:,n).^2)-1) + e);
end


close all;


figure;
vid = VideoWriter('moving_kink_video.mp4', 'MPEG-4');
vid.FrameRate = 30;
open(vid)
saveToGif(x,phi,vid,[-L,L,-2,2],150,v,gamma,T)


figure;
plot(x,phi(:,1),LineWidth=2)
hold on;
plot(x,phi(:,round(end/2)),LineWidth=2)
plot(x,phi(:,end),LineWidth=2)
plot(x,tanh(gamma*x),'--')
plot(x,tanh(gamma*(x-v*T/2*ones(size(x)))),'--')
plot(x,tanh(gamma*(x-v*T*ones(size(x)))),'--')
legend('$t=0$','$t=T/2$','$t=T$','$\tanh(\gamma x)$','$\tanh(\gamma (x-vT/2))$','$\tanh(\gamma (x-vT))$','Interpreter','latex','Location','best','FontSize',20)
xlabel('$x$','Interpreter','latex','FontSize',24)
ylabel('$\phi$','Interpreter','latex','FontSize',24)

function saveToGif(x,phi,vid,axs,frames,v,gamma,T)
    plt1 = plot(x,phi(:,1),'LineWidth',2);
    hold on;
    plt2 = plot(x,tanh(gamma*x),'c--');
    axis(axs)
    xlabel("$x$",Interpreter='latex',FontSize=14)
    ylabel("$\phi$",Interpreter='latex',FontSize=14)
    title('Animation of Moving Kink','FontSize',20)
    legend('Numerical Approximation','Exact Solution',FontSize=14)
    for i = 1:frames
        frame= getframe(gcf);
        %exportgraphics(gcf,name,'Append',true);
        writeVideo(vid, frame.cdata);
        axis(axs)
        plt1.YData = phi(:,round(i*end/frames)); % Update plot data for the current frame
        plt2.YData = tanh(gamma*(x-v*round(i)*T/frames*ones(size(x)))); % Update plot data for the current frame
        drawnow; % Update the figure window
        %pause(0.1); % Pause to control the frame rate
    end
end

close(vid);
