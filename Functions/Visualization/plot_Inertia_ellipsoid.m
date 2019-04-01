%% %% Plot the equivalent uniform ellipsoid of the given inertial parameters
% 2019 Taeyoon Lee

%% Inputs
% [Name]      [Description]                                                                     [Size]
%  T_sb        Transformation matrix from stationary frame {s} to body attached frame {b}        4*4
%  phi_b       Vectorized inertial parameters expressed in body attached frame {b}              10*1
%  trans       Face alpha of the ellipsoid                                                        1
%  color       Face color of the ellipsoid in RGB                                                1*3

%% Implementation
function plot_Inertia_ellipsoid(T_sb, phi_b, trans, color)

P_b = inertiaVecToPinertia(phi_b);
m   = P_b(4,4);
p_b = P_b(1:3,4)/m;
SigmaC_b = P_b(1:3,1:3) - m * p_b * p_b';
[R, S, ~] = eig(SigmaC_b);
s = [S(1,1);S(2,2);S(3,3)];
if det(R)<0
    R_temp = R;
    R(:,1) = R_temp(:,2);
    R(:,2) = R_temp(:,1);
    s(1) = S(2,2);
    s(2) = S(1,1);
end
s(find(s <0)) = eps;
radii = sqrt(5*s/m);

T_bc = [   R  , p_b;
         0,0,0,   1];
T_sc = T_sb * T_bc;

[xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3),200);
R_sc = T_sc(1:3,1:3); p_sc = T_sc(1:3,4);
a = kron(R_sc(:,1),xc); b = kron(R_sc(:,2),yc); c = kron(R_sc(:,3),zc);
data = a+b+c; n = size(data,2);
x = data(1:n,:)+p_sc(1); y = data(n+1:2*n,:)+p_sc(2); z = data(2*n+1:end,:)+p_sc(3);


surf(x,y,z,'EdgeColor','none','FaceColor',color,'FaceAlpha',trans);
axis equal; 
draw_SE3(T_sb);

end