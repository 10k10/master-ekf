clc;
clear all;
close all;

normalize_samples = 1; % true or false
max_samples = 1e6;

% Nomenclatura utilizada x_yvz:
% x: t = true,       m = measured
% y: s = sun sensor, m = magnetometer
% z: b = body,       n = inertial frame
t_svb = csvread('sim_data/true_sun.csv')';
m_svb = csvread('sim_data/meas_sun.csv')';
t_mvb = csvread('sim_data/true_mag.csv')';
m_mvb = csvread('sim_data/meas_mag.csv')';
qn = csvread('sim_data/qn.csv')';

samples = size(t_svb,2);

if normalize_samples
    for n=1:samples
        t_svb(:,n) = t_svb(:,n)/norm(t_svb(:,n));
        m_svb(:,n) = m_svb(:,n)/norm(m_svb(:,n));
        t_mvb(:,n) = t_mvb(:,n)/norm(t_mvb(:,n));
        m_mvb(:,n) = m_mvb(:,n)/norm(m_mvb(:,n));
    end
end

if samples > max_samples
    samples = max_samples;
    t_svb = t_svb(:,1:max_samples);
    m_svb = m_svb(:,1:max_samples);
    t_mvb = t_mvb(:,1:max_samples);
    m_mvb = m_mvb(:,1:max_samples);
    qn = qn(:,1:max_samples);
end

% Pesos das medidas
a = [0.5 0.5];

% Vetor de quaternions estimados
q_bar = zeros(4,samples);
q_error = zeros(samples,1);

contador_acertos = 0;
exec_time = 0;

for i=1:samples
    %-----------------------------------
    %--------------- DATA --------------
    %-----------------------------------
    
    b1 = m_svb(:,i);
    b2 = m_mvb(:,i);
    
    % Rotaciona as medidas inerciais conforme os quaternions do modelo
    qn1 = qrotate([0; b1],qconj(qn(:,i)));
    qn2 = qrotate([0; b2],qconj(qn(:,i)));
    
    % Remove 0 da frente após rotação
    n1 = qn1(2:4);
    n2 = qn2(2:4);
        
    %-----------------------------------
    %-------------- QUEST --------------
    %-----------------------------------   
    
    [q_bar(:,i),t] = QUEST(a,[b1 b2],[n1 n2],'QUEST-nr');
    exec_time = exec_time + t;
    
    q_bar(1,i) = -q_bar(1,i); % Misterioso, no modelo simulado não precisa
    
    %-----------------------------------
    %-------------- ERROR --------------
    %-----------------------------------
    
    q1 = qn(:,i);
    q2 = q_bar(:,i);
    
    qb = [0 1 0 0]';
    
    qb1 = qrotate(qb,q1);
    qb2 = qrotate(qb,q2);

    q_error(i) = real(acos(dot(qb1,qb2)/(norm(qb1)*norm(qb2)))*180/pi);
    
    
    if q_error(i) < 1e-3
        contador_acertos = contador_acertos+1;
    end
          
end

disp('Exec. time medio:')
disp(exec_time/samples)

disp('Porcentagem de acerto:')
disp(contador_acertos/samples*100)

figure(2)
subplot(2,2,1)
plot(q_bar(1,:),'linewidth',2)
hold on
grid on
plot(qn(1,:),'--','linewidth',2)
ylim([-1 1])

subplot(2,2,2)
plot(q_bar(2,:),'linewidth',2)
hold on
grid on
plot(qn(2,:),'--','linewidth',2)
ylim([-1 1])

subplot(2,2,3)
plot(q_bar(3,:),'linewidth',2)
hold on
grid on
plot(qn(3,:),'--','linewidth',2)
ylim([-1 1])

subplot(2,2,4)
plot(q_bar(4,:),'linewidth',2)
hold on
grid on
plot(qn(4,:),'--','linewidth',2)
ylim([-1 1])

figure(3)
plot(q_error,'linewidth',2)
grid on