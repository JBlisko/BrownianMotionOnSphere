clc
close all
clear all

%% Particle Properties
n = 100;                                 % Number of Particles

phi = acos(2*rand(n,1)-1)/2;
theta = 2 * pi * rand(n,1);
R = 5.3e-05;                       % Droplet Radius
p_r = 1 * 10^(-6);                  % Paricle Radius (m)
cl = 10^(-6);                    %Characteristic Length (um)

x = R*sin(phi).*cos(theta)/cl;
y = R*sin(phi).*sin(theta)/cl;
z = R*cos(phi)/cl;

r = [x y z];                            % Particle Coordinates Normalized
% r = [0 0 R; 0.02*R 0.02*R 0.9996*R; R 0 0];
C_0 = [0 0 0];                          % Water Droplet Center Position
I = eye(3);                             % Identity Matrix

x_rand = 2*rand(n,1)-1;
y_rand = 2*rand(n,1)-1;
z_rand = 2*rand(n,1)-1;

p_rand = [x_rand y_rand z_rand];
p_norm = zeros(n,3);
P_i = zeros(n,3);

for ii = 1:n
    dev_ten_i = I-([r(ii,:)*r(ii,1);r(ii,:)*r(ii,2);r(ii,:)*r(ii,3)]/sum(r(ii,:).^2));
    P_i(ii,:) = p_rand(ii,:)*dev_ten_i/norm(p_rand(ii,:)*dev_ten_i);
end

%% Particle Plotting
fig1 = figure;
[u,v,w] = sphere;
u = u(11:end,:);
v = v(11:end,:);
w = w(11:end,:);
s = surf(R*u/cl,R*v/cl,R*w/cl,'FaceAlpha',0.25,'FaceColor','c','edgecolor','none');
axis equal
xlim([-R*1.2/cl R*1.2/cl])
ylim([-R*1.2/cl R*1.2/cl])
zlim([0 R*1.2/cl])
view(0,90);
hold on
h = plot3(r(:,1),r(:,2),r(:,3),'r.','MarkerSize',30);
g = quiver3(r(:,1),r(:,2),r(:,3),P_i(:,1),P_i(:,2),P_i(:,3),0.2);
g.Color = 'black';

%% Particle Motion

d_W = randn(n,3);                   % Normally Distributed Walk

rho = 998;                          % Density of Water (kg/m^3)
eta = 10^-3;                        % Dynamic Viscosity of Fluid (Ns/m^2)
zeta = (6*pi*eta*p_r);              % Friction Coefficient
T = 293;                            % Temperature (K)
k_b = 1.38064852 * 10^(-23);        % Boltzmann Constant (um^2*kg/s^2*K)
D = (T*k_b) / zeta;                 % Diffusivity

gamma = 0.073;                      % Surface Tension (N/m)
%{
tau = 100*eta*cl/gamma;
dt = 0.1*tau;
dt_n = .001*tau;
%}
dt=10^-2;
k = 25*cl^2;                          % Surface Constraint Spring Constant (N/um)
k_ev = 1e-1;                        % Excluded Volume Spring Constant (N/m)

F_s_i = zeros(n,3);                 % Spring Force Initialize
F_ev = zeros(n,3);                  % Excluded Volume Initialize
F_hyd = zeros(n,3);                 % Hydrodynamic Force Initialize
F_ij_22 = zeros(n,3);               % Quadropole-Quadropole Initialize
T_ij_22 = zeros(n,3);               % Quadropole-Quadropole Torque Initialize
dist = zeros(n);

Q_i_0 = 0.245;                       % Fixed Monopolar Rise for i
Q_j_0 = 0.245;                       % Fixed Monopolar Rise for j
Q_i_2 = 0.245;                       % Fixed Quadropolar Rise for i
Q_j_2 = 0.245;                       % Fixed Quadropolar Rise for j

r_in = [x y z];                     % Test to see the change in position

for t = 1:300
    delete(h)
    delete(g)
    delete(s)
    
    d_W = randn(n,3);
    
    for i = 1:n
        temp_1 = [0 0 0];
        temp_2 = [0 0 0];
        temp_3 = [0 0 0];
        temp_4 = [0 0 0];
        
        r_n_i = norm(r(i,:)-C_0);
        F_s_i(i,:) = -k*(r_n_i-R/cl)*(r(i,:)-C_0)/r_n_i;
        
        for j = 1:n
            d = (R/p_r)*acos(dot((r(i,:)-C_0),(r(j,:)-C_0))/...
                (sqrt(sum((r(i,:)-C_0).^2)*sum((r(j,:)-C_0).^2))));
            dist(i,j) = R * acos(dot((r(i,:)-C_0),(r(j,:)-C_0))/...
                (sqrt(sum((r(i,:)-C_0).^2)*sum((r(j,:)-C_0).^2))));
            
            if j~=i
                % Quadropole-Quadropole Force
                dev_ten_i = I-([r(i,:)*r(i,1);r(i,:)*r(i,2);r(i,:)*r(i,3)]/sum(r(i,:).^2));
                d_ij = ((r(j,:)*dev_ten_i)-(r(i,:)*dev_ten_i))/norm((r(j,:)*dev_ten_i)-(r(i,:)*dev_ten_i));
                phi_i = atan2(dot(r(i,:)/r_n_i,cross(P_i(i,:),d_ij)),dot(P_i(i,:),d_ij));
                
                r_n_j = norm(r(j,:)-C_0);
                dev_ten_j = I-([r(j,:)*r(j,1);r(j,:)*r(j,2);r(j,:)*r(j,3)]/sum(r(j,:).^2));
                d_ji = ((r(j,:)*dev_ten_j)-(r(i,:)*dev_ten_j))/norm((r(j,:)*dev_ten_j)-(r(i,:)*dev_ten_j));
                phi_j = atan2(dot(r(j,:)/r_n_j,cross(P_i(j,:),d_ji)),dot(P_i(j,:),d_ji));
                
                F_ij_22(i,:) = (48*pi*gamma*p_r^4*Q_i_2*Q_j_2/(d^5))*...
                    (d_ij*cos(2*(phi_i+phi_j))-cross(r(i,:)/r_n_i,d_ij)*sin(2*(phi_i+phi_j)))/cl;
                
                temp_1 = temp_1 + F_ij_22(i,:);
                
                % Excluded Volume Force
                delta = (d-2*p_r/cl);
                
                if delta < 0
                    F_ev(i,:) = -.0025*k_ev*delta*d_ij*cl;
                    temp_2 = temp_2 + F_ev(i,:);
                end
                
                % Hydrodynamic Force
                F_hyd(i,:) = -6*pi*eta^2*(Q_j_0/d)^2*d_ij/rho;
                temp_3 = temp_3 + F_hyd(i,:);
                
                % Quadropolar Torque
                T_ij_22(i,:) = -24*pi*gamma*p_r^4*Q_i_2*Q_j_2*sin(2*(phi_i+phi_j))*r(i,:)/(d^4*r_n_i)/(cl^2);
                temp_4 = temp_4 + T_ij_22(i,:);
            end
        end
        
        F_ij_22(i,:) = temp_1;
        F_ev(i,:) = temp_2;
        F_hyd(i,:) = temp_3;
        T_ij_22(i,:) = temp_4;
    end
    
    % Change in Position
    d_r = sqrt(2*D*dt)*d_W/cl + dt*(F_s_i+F_ij_22+F_ev)/zeta;
    r = r + d_r;
    
    %{
    for k = 1:n
       if r(k,3)/R < 0
           r(k,3) = -r(k,3);
       end
    end
    %}
    
    % Change in Orientation
    w_i = cl^2*T_ij_22/(4*pi*eta*p_r^3);
    d_p = cross(w_i,P_i)*dt;
    for iii = 1:n
        P_i(iii,:) = (P_i(iii,:)+d_p(iii,:))/norm(P_i(iii,:)+d_p(iii,:));
    end
    
    % Updated Plot
    h = plot3(r(:,1),r(:,2),r(:,3),'r.','MarkerSize',30);
    hold on
    g = quiver3(r(:,1),r(:,2),r(:,3),P_i(:,1),P_i(:,2),P_i(:,3),0.2);
    g.Color = 'black';
    [u,v,w] = sphere;
    u = u(11:end,:);
    v = v(11:end,:);
    w = w(11:end,:);
    s = surf(R.*u/p_r,R.*v/p_r,R.*w/p_r,'FaceAlpha',0.25,'FaceColor','c','edgecolor','none');
    axis equal
    xlim([-R*1.2/p_r R*1.2/p_r])
    ylim([-R*1.2/p_r R*1.2/p_r])
    zlim([0 R*1.2/p_r])
    
    fr_get(t) = getframe(fig1, [10 10 520 400]);
end

gif_Writer = VideoWriter('Animation 2');
gif_Writer.FrameRate = 30;
open(gif_Writer);
writeVideo(gif_Writer,fr_get)
close(gif_Writer);

hold off

close all