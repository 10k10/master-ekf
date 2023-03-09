clc;
clear all;
close all;

Ix = 2.0; % Moments of inertia about x
Iy = 3.0; % Moments of inertia about y
Iz = 1.2; % Moments of inertia about z

% Constantes do sistema
kx = (Iy - Iz)/Ix;
ky = (Iz - Ix)/Iy;
kz = (Ix - Iy)/Iz;

% Parâmetros de simulação
Ts = 1e-3;
sim_time = 100;
sim_steps = sim_time/Ts;

% Estados:
% x0 = q0, x1 = q1, x2 = q3, x4 = q4
% x5 = wx, x6 = wy, x7 = wz

% Condições iniciais
q_0 = [0 1 0 0]';
omega_0 = [-0.1 -0.03 0.5]'; % rad/s

x0 = [q_0; omega_0];

x_m = zeros(7, sim_steps);
x_m(:,1) = x0;

% Simula o modelo
for k=2:sim_steps
    q = 0.5.*[x_m(1,k-1) -x_m(2,k-1) -x_m(3,k-1) -x_m(4,k-1);
              x_m(2,k-1)  x_m(1,k-1) -x_m(4,k-1)  x_m(3,k-1);
              x_m(3,k-1)  x_m(4,k-1)  x_m(1,k-1) -x_m(2,k-1);
              x_m(4,k-1) -x_m(3,k-1)  x_m(2,k-1)  x_m(1,k-1)]*[0; x_m(5,k-1); x_m(6,k-1); x_m(7,k-1)];
                    
    omega = [kx*x_m(6,k-1)*x_m(7,k-1);
             ky*x_m(5,k-1)*x_m(7,k-1);
             kz*x_m(5,k-1)*x_m(6,k-1)];
    
    x_m(:,k) = x_m(:,k-1) + Ts.*[q; omega];
end



% EKF Variables
x = zeros(7, sim_steps);
x0a = x0(1:4);
x0b = x0(5:7)*(1+0.2*(rand-0.5));
x(:,1) = [x0a; x0b];

ap = 0.7;
aq = 0.0005;
ar = 0.005;

P = ap.*eye(7);
Q = aq.*eye(7);
R = ar.*eye(4);

H = [eye(4) zeros(4,3)];
F = zeros(7);
Z = zeros(4,sim_steps);

q_ekf_error = zeros(1,sim_steps);
w_ekf_error = zeros(3,sim_steps);

dv = 0.02;

% Vetores de medidas no referencial inercial
n1 = [1 0 0]';
n2 = [0 1 0]';

% Pesos das medidas
a = [0.5 0.5];

% Vetor de quaternions estimados
q_bar = zeros(sim_steps,4);
q_error = zeros(sim_steps,1);

contador_acertos = 0;

tic

for i=2:sim_steps
    k = i;
    % Rotaciona as medidas inerciais conforme os quaternions do modelo
    qn = x_m(1:4,i)';
    qb1 = qrotate([0 n1'],qconj(qn)) + normrnd(0,dv/2,4,1);
    qb2 = qrotate([0 n2'],qconj(qn)) + normrnd(0,dv/2,4,1);
    
    % Remove 0 da frente após rotação
    b1 = qb1(2:4);
    b2 = qb2(2:4);
    
    %-----------------------------------
    %-------------- QUEST --------------
    %-----------------------------------
    
    q_bar(i,:) = QUEST(a,[b1 b2],[n1 n2]);
    %q_bar(i,:) = qMethod(a,[b1 b2],[n1 n2]);
    
    if qn(1) < 0
        q_bar(i,:) = -q_bar(i,:);
    end 
    
    q_bar(i,:) = [q_bar(i,4) q_bar(i,1) q_bar(i,2) q_bar(i,3)];
    
    %-----------------------------------
    %--------------- EKF ---------------
    %-----------------------------------
    
    q = 0.5.*[x(1,k-1) -x(2,k-1) -x(3,k-1) -x(4,k-1);
              x(2,k-1)  x(1,k-1) -x(4,k-1)  x(3,k-1);
              x(3,k-1)  x(4,k-1)  x(1,k-1) -x(2,k-1);
              x(4,k-1) -x(3,k-1)  x(2,k-1)  x(1,k-1)]*[0; x(5,k-1); x(6,k-1); x(7,k-1)];
                    
    omega = [kx*x(6,k-1)*x(7,k-1);
             ky*x(5,k-1)*x(7,k-1);
             kz*x(5,k-1)*x(6,k-1)];
    
    % Faz a predição dos estados usando fd(x[k-1],u[k-1])
    x(:,k) = x(:,k-1) + Ts.*[q; omega];
    
    % Computa o jacobiano
    
    q0 = x(1,k-1);
    q1 = x(2,k-1);
    q2 = x(3,k-1);
    q3 = x(4,k-1);
    wx = x(5,k-1);
    wy = x(6,k-1);
    wz = x(7,k-1);
    
    F = [0  -wx -wy -wz -q1   -q2  -q3;
         wx  0   wz -wy  q0   -q3   q2;
         wy -wz  0   wx  q3    q0  -q1;
         wz  wy -wx  0  -q2    q1   q0;
         0   0   0   0   0   kx*wz kx*wy;
         0   0   0   0 ky*wz   0   ky*wx;
         0   0   0   0 kz*wy kz*wx  0];
     
    F = eye(7) + F;
     
    P = F*P*F' + Q;
    Z(:,k) = qn' + normrnd(0,dv,4,1); % add noise
    %Z = Z/sqrt(Z(1)^2 + Z(2)^2 + Z(3)^2 + Z(4)^2); %normalize
    y = Z(:,k) - H*x(:,k);
    %y = H*x_m(:,k) - H*x(:,k);
    S = H*P*H' + R;
    K = P*H'*S^-1;
    x(:,k) = x(:,k) + K*y;
    P = (eye(7)-K*H)*P;
     
    q_m = x_m(1:4,k);
    q_y = Z(:,k);
    q_e = x(1:4,k);  
    
    q1 = qn;
    q2 = q_bar(i,:);

    % Calcula o ângulo principal (erro) entre o quaternion estimado e o real.
    q_error(i) = qangle(qmult(q1,qconj(q2)))*180/pi;
    q_ekf_error(i) = qangle(qmult(q1,qconj(q_e)))*180/pi;
    
    if q_error(i) < 0.1
        contador_acertos = contador_acertos+1;
    end
end

t2 = toc;

fprintf('Tempo de execução: %.4fs\r\n',t2);
fprintf('Fit: %.4f%%\r\n',contador_acertos/sim_steps*100);

figure(1)
subplot(2,1,1);
plot(x_m(1,:),'linewidth',2)
hold on
grid on
plot(x_m(2,:),'linewidth',2)
plot(x_m(3,:),'linewidth',2)
plot(x_m(4,:),'linewidth',2)
legend('$q_0$','$q_1$','$q_2$','$q_3$','Interpreter','latex','FontSize',16)

subplot(2,1,2);
plot(x_m(5,:),'linewidth',2)
hold on
grid on
plot(x_m(6,:),'linewidth',2)
plot(x_m(7,:),'linewidth',2)
legend('$\omega_x$','$\omega_y$','$\omega_z$','Interpreter','latex','FontSize',16)

figure(2)
subplot(2,2,1)
plot(q_bar(:,1),'linewidth',2)
hold on
grid on
plot(x_m(1,:),'--','linewidth',2)
plot(x(1,:),'linewidth',2)
ylim([-1 1])

subplot(2,2,2)
plot(q_bar(:,2),'linewidth',2)
hold on
grid on
plot(x_m(2,:),'--','linewidth',2)
plot(x(2,:),'linewidth',2)
ylim([-1 1])

subplot(2,2,3)
plot(q_bar(:,3),'linewidth',2)
hold on
grid on
plot(x_m(3,:),'--','linewidth',2)
plot(x(3,:),'linewidth',2)
ylim([-1 1])

subplot(2,2,4)
plot(q_bar(:,4),'linewidth',2)
hold on
grid on
plot(x_m(4,:),'--','linewidth',2)
plot(x(4,:),'linewidth',2)
ylim([-1 1])

figure(3)
plot(q_error,'linewidth',2)
hold on
grid on
plot(q_ekf_error,'linewidth',2)