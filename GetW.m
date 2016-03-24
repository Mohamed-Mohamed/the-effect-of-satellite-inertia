function [ F ] = GetW(x0 )
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg

load('I.mat');
F =[ x0(2)*(Ixz*x0(1) + Iyz*x0(2) + Izz*x0(3)) - x0(3)*(Ixy*x0(1) + Iyy*x0(2) + Iyz*x0(3))-A1(1);
       x0(3)*(Ixx*x0(1) + Ixy*x0(2) + Ixz*x0(3)) - x0(1)*(Ixz*x0(1) + Iyz*x0(2) + Izz*x0(3))-A1(2);
       x0(1)*(Ixy*x0(1) + Iyy*x0(2) + Iyz*x0(3)) - x0(2)*(Ixx*x0(1) + Ixy*x0(2) + Ixz*x0(3))-A1(3)];
end

