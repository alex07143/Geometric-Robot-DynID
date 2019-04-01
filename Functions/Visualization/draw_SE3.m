%% %% Draw x-y-z axes in space
% 2019 Taeyoon Lee

%% Inputs
% [Name]      [Description]                                                                     [Size]
%  T_sb        Transformation matrix from stationary frame {s} to the frame {b}                  4*4

%% Implementation
function draw_SE3(T_sb)
p=T_sb(1:3,4);
ax=p+T_sb(1:3,1)*0.05;
ay=p+T_sb(1:3,2)*0.05;
az=p+T_sb(1:3,3)*0.05;
hold on;
plot3([p(1),ax(1)], [p(2),ax(2)], [p(3),ax(3)],'r','LineWidth',1.5);
plot3([p(1),ay(1)], [p(2),ay(2)], [p(3),ay(3)],'g','LineWidth',1.5);
plot3([p(1),az(1)], [p(2),az(2)], [p(3),az(3)],'b','LineWidth',1.5);
end