% 2 cell FitzHugh Nagumo model

% TO RUN (at the terminal)
% ./run_... /usr/local/matlab2010b/ <arg1> <arg2> ... &

% function = FHnetwork(N,f,C, seed)

%if ischar(N);      N      = str2double(N);      end;
%...

% Replace random stream with one based on SEED
%RandStream.setDefaultStream(RandStream('mt19937ar','Seed',seed));

%Fname = sprintf('%04d_%02d_%03d_%02.2f_',N,f,C,..);
clc
close all
clear all;
global N I g E c stopper adj mat
N = 2;
stopper = 1;
tic
% fraction of the neurons that are excited
stim = 1; %fraction of stimulated neurons
stimmag=0.12; %magnitude of stimulation
I=zeros(N,1); %stimulation current
k=randperm(N);kk=k(1:ceil(stim*N));
I(kk)=stimmag
%I = 0.12*[1 zeros(1,N-1)];

%function [T,Y] = twocellFHmodel()

tspan = 50;

% system 1 - v1,w1 and s21
% Network of N neurons, out of which f*N are inhibitory and (1-f)*N are
% excitatory, and the connection density is C.


f = 0
C = 1


% first define some random adjacency matrix, which will be NxN, with a
% connection density of C.

iter = 1;

%adj = rand(N,N)>0.5;

% then define a random set of weights where the weight is necessarily 0 if
% the corresponding term in the adjacency matrix is 0.

%g = adj.*rand(N,N);

% for now keep the connection weight matrix same as the adjacency matrix

%g = [0 0.04; 0 0];
% f is the fraction of inhibitory neurons
fig_count = 1;
filepath = '/home/samyuktar/Documents/MATLAB/Images/';

% change the options for the ode solver
%for iter = 1:10
%for C = 0:0.1:1
adj=0.5+0.5*sign(C-rand(N));
for i=1:N,adj(i,i)=0;end;

%    adj  = [0 1 0 0 0 0;1 0 0 1 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 1 0 0];    
sprintf('actual C, %f',sum(sum(adj))/(N*N))
gsyn = 0.01;
g= gsyn*bsxfun(@rdivide,adj,sum(adj,1))

%g = [0,0.09;0,0];
stim;

%    for f = 0:0.1:1

%Synaptic reversal potentials
E=5*ones(1,N); %excitatory
k=randperm(N);kk=k(1:ceil(f*N));
E(kk)=-5
%E = E(randperm(N)); % randomly shuffle this so that which neuron is inhibitory and which is excitatory is not fixed
% write the FH equations of each neuron i which gets input from each of its
% N-1 neighbours (depending on the value of gji)


%[T,Y] = ode15s(@oscnetwork_opt2,tspan,zeros(1,N*(N+1)),options);

h = 0.001;                                              % step size
%    tfinal = 10;
tot_time = ceil(tspan/h);                                      %ceil is round up
%    initial values vector    
Y0 = zeros(N*(N+1),1);
% state vector: y = [ X dX Z dZ Phi dPhi ]'
%    dynamics = @(t,y) oscnetwork_opt(t,y);
%    t = zeros(tot_time,1);
%   T = 1:tot_time;
T = 1:h:tspan;
dynamics = @(t,y) oscnetwork_opt_vec_out_fast(t,y);

%    dynamics =  oscnetwork_opt(t,Y);
% result matrices
Y = zeros(tot_time,N*(N+1));
% set initial condition and integrate
Y(1,:) = Y0;
    
for i= 1:tot_time
          % update time and dynamics
      [ T(i+1), Y(i+1,:) ] = integrator_rk4(h,dynamics,T(i),Y(i,:)');
end
    
    


%[T,Y] = rk4code2(N,adj,tspan);
%[T,Y] = ode45(@oscnetwork_opt,tspan,zeros(1,(2*N)+sum(sum(adj))),options);
size(Y)
% Plot voltage vs time
% figure(1)
grid on


list  = 1:N;
col=hsv(N);

%voltage_std_stt = [];

% plot all the plots
figure(fig_count);
count = 0;






for i = 1:N
    count = count+1;
    grid on
    if E(count)<0
        subplot(2,1,1)
        ax1 = plot(T,Y(:,i),'color',col(count,:));
        %voltage_std_stt = [voltage_std_stt Y(end,i)];
        hold on

    elseif E(count)>0
        subplot(2,1,2) 
        ax2 = plot(T,Y(:,i),'color',col(count,:));
        %voltage_std_stt = [voltage_std_stt Y(end,i)];
        hold on

    end
%    legendInfo{count} = ['Voltage Cell No: ' num2str(count)]; 
    xlabel('Time(s)');
    ylabel('Voltage (mV)');
    title('Voltage vs Time-for N cells');       
    %legendInfo{count} = ['Voltage Cell No: ' num2str(count)]; 


end
hold off
%legendInfo;
%legend(legendInfo);
%legend(ax2,2, legendInfo2);

%saveas(fig_count,[filepath,'N',int2str(N),'f',num2str(floor(f)),'p',num2str((f-floor(f))*10),'C',num2str(floor(C)),'p',num2str((C-floor(C))*10),'fig',int2str(fig_count),'tspan',num2str(tspan),'iter',num2str(iter)],'fig');
%saveas(fig_count,[filepath,'N',int2str(N),'f',num2str(floor(f)),'p',num2str((f-floor(f))*10),'C',num2str(floor(C)),'p',num2str((C-floor(C))*10),'fig',int2str(fig_count),'tspan',num2str(tspan),'iter',num2str(iter)],'png');
fig_count = fig_count+1;

%    end
%end


thresh = -c;

%for thresh = -0.1:0.1:0.5
%thresh = 0.1; %threshold to decide if the cell is active or not
%{
figure(fig_count);
fe = zeros(1,length(T));
fi = zeros(1,length(T));
for i = 1:length(T)
    count_act_inh = 0;
    count_act_exc = 0;
    for j = 1:N
        if (Y(i,j) > thresh) && E(j)<0 
            count_act_inh = count_act_inh + 1;
        elseif (Y(i,j) > thresh) && E(j)>0
            count_act_exc = count_act_exc + 1;
        end
    end
    fe(i) = count_act_exc/(N*(1-f));
    fi(i) = count_act_inh/(N*f);
end
fe;
fi;
subplot(2,1,1);
plot(T,fe,'g');
hold on;
plot(T,fi,'r');
hold off;
time_step = T(4)-T(3);
s= sprintf('time step, %f',time_step);
disp(s);
xlabel('Time(s)');ylabel('fraction of active inhibitory and excitatory neurons');
title(['Fraction of active inhibitory and excitatory neurons vs time',int2str(thresh)]);
legend('fraction-excitatory','fraction-inhibitory');

subplot(2,1,2);
plot(fe,fi,'k');
xlabel('fraction of active excitatory neurons');
ylabel('fraction of active inhibitory neurons');
title('Phase plot of Fraction of active excitatory to inhibitory neurons, to total number N of neurons');
%}
%    txtf = ['f',num2str(floor(f)),'p',num2str((f-floor(f)))]
%saveas(fig_count,[filepath,'phaseN',int2str(N),'f',num2str(floor(f)),'p',num2str((f-floor(f))*10),'C',num2str(floor(C)),'p',num2str((C-floor(C))*10),'fig',int2str(fig_count),'tspan',num2str(tspan),'iter',num2str(iter)],'fig');
%saveas(fig_count,[filepath,'phaseN',int2str(N),'f',num2str(floor(f)),'p',num2str((f-floor(f))*10),'C',num2str(floor(C)),'p',num2str((C-floor(C))*10),'fig',int2str(fig_count),'tspan',num2str(tspan),'iter',num2str(iter)],'png');
%    saveas(fig_count,[filepath,'fefiN100p5t30',int2str(fig_count)],'png');
fig_count = fig_count+1;
%end
    
    
%end

%save(Fname,'Y','...')

toc
