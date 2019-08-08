% This code computes the spontaneous curvature at the edges of a helicoid
% and its bending energy for various inner ramp radii.
% It requires the accompagning comsol file (hel00.mph), Comsol Multiphysics
% version 5.3a or higher with the module "Chemical engineering", and must 
% be ran with within a "Comsol with Matlab" instance.

% This code was developed for the research article:
% Chabanon M. & Rangamani P., "?Geometric coupling of helicoidal ramps and
% curvature-inducing proteins in organelle membranes", Journal of the Royal
% Society Interface (2019).

% Code author: Morgan Chabanon

clear all;  %close all
tic

name_model='hel00'; %name of comsol model
model = mphload(name_model);  %loads the comsol model


%R0=[0 0.02 0.04 0.06 0.08 0.1]; % values of the inner ramp radius
R0=linspace(0,0.1,50);

%The helicoid pitch is 2*pi*n_turns*p
p = 0.05; %typically between 0.01 and 0.1. 
n_turns = 1; %number of turns

c0=1; c1=5;  % Spontaneous curvature imposed at the bottom and top of 
%inner ramps respectively. A linear gradient from C0 to C1 is imposed at 
%the outer ramp

% Invariant
model.param.set('n', n_turns); % number of turns
model.param.set('p', p); % p=H/(2*pi*n)
model.param.set('c0', c0); % Spontaneous curvature at lower boundary
model.param.set('c1', c1); % Difference between upper and lower boundary C
model.param.set('L', 1); % Helicoid diameter

% Initialization
N=numel(R0);
n_data = 101; %number of elements on the edges +1 (imposed in comsol)
r_top = NaN*ones(N,n_data);
c_top = r_top;
r_bot = r_top;
c_bot = r_top;
z_cen = r_top;
c_cen = r_top;
cp_top = NaN*ones(N,1);
cp_bot = cp_top;
I=NaN*ones(1,n_data);
W = NaN*ones(N,1);
A = W;

% Parameteric study
for i = 1:N  %takes of the order of 5 to 15 sec per iteration on a desktop computer

    r0 = R0(i);   
    model.param.set('r0', r0);
    
    model.study('std1').run; %run comsol model
    
    %Evaluates C along the top boundary along r
    data = mpheval(model,'c','selection',3,'edim','edge');
    [r_top(i,:) , I] = sort(data.p(1,:),2);
    c_top(i,:) = data.d1(I);
    
    %Evaluates C along the bottom boundary along r
    data = mpheval(model,'c','selection',2,'edim','edge');
    [r_bot(i,:) , I] = sort(data.p(1,:),2);
    c_bot(i,:) = data.d1(I);
    
    %Evaluates C along the iner ramp
    data = mpheval(model,'c','selection',4,'edim','edge');
    [z_cen(i,:) , I] = sort(data.p(3,:),2);
    c_cen(i,:) = data.d1(I);
    
    %Evaluates C at the top point of the inner ramp
    data = mpheval(model,'c','selection',4,'edim','point');
    cp_top(i) = data.d1;
    
    %Evaluates C at the bottom point of the inner ramp
    data = mpheval(model,'c','selection',3,'edim','point');
    cp_bot(i) = data.d1;
    
    %Evaluates the bending energy by integratin over the surface
    W(i)=mphint2(model,'c^2-0.9*K' ,'surface'); % W/(k c_0^2 L^2)
    A(i)=mphint2(model,'1','surface'); % Area

    disp('-----')
    disp(['r_0 = ' num2str(r0)]);
    disp([num2str(floor(i/N*100)) '% progress']);

end

toc
clear data i model ans

% Saves data in .mat format
filename = strcat(name_model, '_', string(datetime, 'HHmmss'));
save(filename)


% Visualization of the spontaneous curvature at the upper and lower points
% of the inner ramp
figure
plot(R0,cp_top, '^', R0,cp_bot, 'v')
xlabel('r_0/L') ; ylabel('C(r_0)/C_0')
legend('top', 'bottom')
set(gca,'FontSize', 16, 'position', [0.18 0.20 0.75 0.75])
axis square ; grid minor
axis([0 0.1 -80*[1 -1]])

% Visualization of the bending energy
figure
semilogy(R0, W./A, 'o','color', [1 1 1]*0.6)
axis([0 0.1 1 1e4])
xlabel('r_0/L') ; ylabel('W_b/(kC_0^2\Omega)')
set(gca,'FontSize', 16, 'position', [0.18 0.20 0.75 0.75])
axis square ; grid minor

   

