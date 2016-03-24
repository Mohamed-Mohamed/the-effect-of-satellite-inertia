function [ IP ] = Ip ( I, Start, End, m, state )
% This function is used to get parallel moment of inertia matrix
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% inputs
% I   : the inertia matrix (3x3)
% Start : point of I matrix (1x3)
% End : point of IP matrix (1x3)
% m : mass (1x1)
% state : if P from IG to IP if G from IP to IG
%% outputs
% IP : the parallel inertia matrix (3x3)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------
r=End-Start;
if state == 'P'
    IP(1,1)=I(1,1)+m*(r(2)^2+r(3)^2);
    IP(2,2)=I(2,2)+m*(r(1)^2+r(3)^2);
    IP(3,3)=I(3,3)+m*(r(1)^2+r(2)^2);
    IP(1,2)=I(1,2)-m*r(1)*r(2);
    IP(2,1)=I(2,1)-m*r(1)*r(2);
    IP(1,3)=I(1,3)-m*r(1)*r(3);
    IP(3,1)=I(3,1)-m*r(1)*r(3);
    IP(2,3)=I(2,3)-m*r(2)*r(3);
    IP(3,2)=I(3,2)-m*r(2)*r(3);
elseif state == 'G'
    IP(1,1)=I(1,1)-m*(r(2)^2+r(3)^2);
    IP(2,2)=I(2,2)-m*(r(1)^2+r(3)^2);
    IP(3,3)=I(3,3)-m*(r(1)^2+r(2)^2);
    IP(1,2)=I(1,2)+m*r(1)*r(2);
    IP(2,1)=I(2,1)+m*r(1)*r(2);
    IP(1,3)=I(1,3)+m*r(1)*r(3);
    IP(3,1)=I(3,1)+m*r(1)*r(3);
    IP(2,3)=I(2,3)+m*r(2)*r(3);
    IP(3,2)=I(3,2)+m*r(2)*r(3);
end
end