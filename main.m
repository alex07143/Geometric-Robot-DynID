%% %% Sample code for human inertial parameter identification experiment using different regularizers
% 2019 Taeyoon Lee
close all
clear all
clc

addpath(genpath('Functions'))
addpath(genpath('Data'))

%% Loading data & inititalizations
% [Name]         [Description]                                                                     [Size]
%  A              Regressor matrix                                                            N_sample*(N_body*10)
%  b              Observation vector                                                          N_sample*1
%  Phi_prior      Column-wise concatenation of vectorized prior inertial parameters                 10*N_body
%  Sigma          Observation variance                                                        N_sample*N_sample
%  gamma          regularization factor                                                               1
%  Q              Bounding ellipsoidal region => [x;1]'* Q * [x;1]' >=0                (cell: N_body) * (matrix: 4*4)

Data       = load('Train_data_human.mat');

N_body     = size(Data.Phi_prior,2);
A          = Data.A;
b          = Data.b;
Phi_prior  = Data.Phi_prior;
Q          = Data.Q;

sigma      = diag([0.119 0.216 1].^2);
Sigma      = kron(eye(length(b)/3), sigma);
gamma      = 1e-2 * trace(A'*inv(Sigma)*A);

cvx_setup;
cvx_solver mosek    % MOSEK solver is highly recommended for improved numerical stability and convergence speed over SDPT3 and SEDUMI.
                    % Using CVX with MOSEK requires a CVX Professional license. 
                    % Academic users can obtain a free CVX Professional license by submitting an Academic License Request.
                    % See < http://cvxr.com/cvx/doc/mosek.html > for more information.  

%% Identification with different regularizers (using provable comparative analysis method)

% Step1: Solve 'regularized' least squares formulation with Euclidean regularizer.
[Params_Euc, LS_error_Euc]             = ID_Euclidean('regularized',A,b,Phi_prior,Sigma,gamma,Q);

% Step2: Store least square error from the previous experiment
c = LS_error_Euc;

% Step3: Solve 'point-to-set' identification formulation using other convex regularizers with least square error bound "c" 
[Params_pullback, LS_error_pullback]   = ID_ConstPullback('point-to-set',A,b,Phi_prior,Sigma,c,Q);
[Params_entropic, LS_error_entropic]   = ID_Entropic('point-to-set',A,b,Phi_prior,Sigma,c,Q);

% Step4: Check that least square errors on the training samples are identical !!
LS_error = [ LS_error_Euc LS_error_pullback LS_error_entropic ];
assert(abs(max(LS_error)-min(LS_error))/mean(LS_error) < 1e-5);

%% Visualization of the identified parameters with equivalent uniform ellipsoids
Body_Frames = load('Body_frames_human.mat');
trans = 0.5; color = 0.8 * rand(N_body,3);

figure(1); clf; hold on;

subplot(1,3,1); hold on;
title('Euclidean','FontSize',16);
for i = 1 : N_body
    plot_Inertia_ellipsoid(Body_Frames.T_sb(:,:,i), Params_Euc.Phi(:,i), trans, color(i,:));
end
view([-200 20]);
axis([-0.6 0.7 -0.3 0.4 -0.1 1.9])
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
lighting phong;
camlight;

subplot(1,3,2); hold on;
title('Const. pullback','FontSize',16);
for i = 1 : N_body
    plot_Inertia_ellipsoid(Body_Frames.T_sb(:,:,i), Params_pullback.Phi(:,i), trans, color(i,:));
end
view([-200 20]);
axis([-0.6 0.7 -0.3 0.4 -0.1 1.9])
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
lighting phong;
camlight;

subplot(1,3,3); hold on;
title('Entropic','FontSize',16);
for i = 1 : N_body
    plot_Inertia_ellipsoid(Body_Frames.T_sb(:,:,i), Params_entropic.Phi(:,i), trans, color(i,:));
end
view([-200 20]);
axis([-0.6 0.7 -0.3 0.4 -0.1 1.9])
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
lighting phong;
camlight;
