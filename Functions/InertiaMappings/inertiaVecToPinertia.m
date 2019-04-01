%% %% Tranform vectorized inertial parameters to pseudo inertia matrix representation
% 2019 Taeyoon Lee

%% Inputs
% [Name]      [Description]                                                                     [Size]
%  phi         Vectorized inertial parameters                                                    10*1

%% Outputs
% [Name]      [Description]                                                                     [Size]
%   P          pseudo inertia matrix                                                             4*4

%% Implementation
function [ P ] = inertiaVecToPinertia( phi )

m    = phi(1);
h    = [phi(2);
        phi(3);
        phi(4)];
Irot = [phi(5)   phi(8)  phi(10);
        phi(8)   phi(6)  phi(9);
        phi(10)  phi(9)  phi(7)];

P = [(1/2)*trace(Irot)*eye(3)-Irot, h;
                    h'            , m];

end
