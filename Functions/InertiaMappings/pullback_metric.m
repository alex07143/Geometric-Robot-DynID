%% %% Construct pullback form of the affine-invariant Riemannian metric defined on P(4) to R^10 coordinate
% 2019 Taeyoon Lee

%% Inputs
% [Name]      [Description]                                                                     [Size]
%  phi         Vectorized inertial parameters                                                    10*1

%% Outputs
% [Name]      [Description]                                                                     [Size]
%   g          matrix representation of the affine-invariant Riem. metric evaluated on phi       10*10
%              (e.x. ds^2 = dphi'* g * dphi)

%% Implementation
function [ g ] = pullback_metric( phi )

g = zeros(10,10);

P_inv = inv(inertiaVecToPinertia(phi));

for i = 1 : 10
    for j = 1 : 10
        e_i = zeros(10,1); e_i(i) = 1;
        e_j = zeros(10,1); e_j(j) = 1;
        V_i = inertiaVecToPinertia(e_i);
        V_j = inertiaVecToPinertia(e_j);
        
        g(i,j) = trace( P_inv * V_i * P_inv * V_j );
    end
end

end