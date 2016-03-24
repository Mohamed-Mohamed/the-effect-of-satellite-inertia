%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg


% this program is used for showing a satalite about the earth
close all; clear all; clc;
%% constants
G=6.67428E-11; % gravitational constant
hz=15000; % simulation frequancy
step_sim=1; % simulation step
%% Intial condition
M=[5.972*10^24,1000];       % [M1,M2]
R1=[0;0;0];                 % position of M(1)
R2=[8000e3;0;0];                 % position of M(2)
I=[0,0,0];              % location of initial axis
V1=[0;0;0];                 % velocity of M(1)
V2=[0;3;7.5]*1e3;                 % velocity of M(2)
% inertia values about the center of M2 (satellite)
Ixx=1000;
Iyy=1000;
Izz=1000;
Ixy=0;
Ixz=0;
Iyz=0;
% the distance between the center of M2 (satellite) relative to axis at M2 (satellite)
xG2=0; yG2=0; zG2=0.1;
% image (inter your URL of the image)
image_file = 'D:\4th year of Aerospace\1st\Orbital Mechanics\AER-427, Orbital Mechanics, Mohamed Mohamed Elsayed,SCE 2, BN 13  By MATLAB\week 11\satalite about earth with tracking/earth.jpg';
%% RK4 parameter
tf=24*3600*0+24*3600/8*1.13;   % final time of soution
dt=0.1*0+100;            % time step
X0=[R1;R2;V1;V2];
B=[0;0;0;0;0;0;0;0;0;0;0;0];
sol(1:12,1)=X0;
order=12;
%% RK4 solution
for n=1:length(0:dt:tf)
    b=G*M(2)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    c=-G*M(1)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    A=[0,0,0,0,0,0,1,0,0,0,0,0; ...
        0,0,0,0,0,0,0,1,0,0,0,0; ...
        0,0,0,0,0,0,0,0,1,0,0,0; ...
        0,0,0,0,0,0,0,0,0,1,0,0; ...
        0,0,0,0,0,0,0,0,0,0,1,0; ...
        0,0,0,0,0,0,0,0,0,0,0,1;...
        -b,0,0,b,0,0,0,0,0,0,0,0; ...
        0,-b,0,0,b,0,0,0,0,0,0,0; ...
        0,0,-b,0,0,b,0,0,0,0,0,0; ...
        -c,0,0,c,0,0,0,0,0,0,0,0; ...
        0,-c,0,0,c,0,0,0,0,0,0,0; ...
        0,0,-c,0,0,c,0,0,0,0,0,0 ];
    [ XX ] = RK4( A,B,sol(1:12,n),dt,n*dt,(n+1)*dt,order );
    sol(1:12,n+1)=XX(1:12,2);
end
R1_x=sol(1,:);
R1_y=sol(2,:);
R1_z=sol(3,:);
R2_x=sol(4,:);
R2_y=sol(5,:);
R2_z=sol(6,:);
V1_x=sol(7,:);
V1_y=sol(8,:);
V1_z=sol(9,:);
V2_x=sol(10,:);
V2_y=sol(11,:);
V2_z=sol(12,:);
%% projected path of the satalite on the earth
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
for P=1:length(R1_x)
    points = intersectLineSphere([0,0,0,R2_x(P)-R1_x(P)+(P-1)*erot*dt,R2_y(P)-R1_y(P)+(P-1)*erot*dt,R2_z(P)-R1_z(P)+(P-1)*erot*dt], [0,0,0,6400e3]);
    XX1(P)=points(2,1);
    YY1(P)=points(2,2);
    ZZ1(P)=points(2,3);
end
%% center of masses parameters
r=[R2_x;R2_y;R2_z]-[R1_x;R1_y;R1_z];                                            % the distance betweem M1 & M2
Rc=(M(1)*[R1_x;R1_y;R1_z]+M(2)*[R2_x;R2_y;R2_z])/sum(M);                        % location of center of masses
Vc=(M(1)*[V1_x;V1_y;V1_z]+M(2)*[V2_x;V2_y;V2_z])/sum(M);                        % Vc = constant
Ac=(M(1)*G*M(2)*r.^3.*r-M(1)*G*M(2)*r.^3.*r)/sum(M);    % for check Ac = 0
for H=1:length(Rc(1,:))
    % acceleration of M1
    acc1(1:3,H)=(-1)^(1)*G*M(1)/(norm([R2_x(H);R2_y(H);R2_z(H)]'-[R1_x(H);R1_y(H);R1_z(H)]'))^3*([R2_x(H);R2_y(H);R2_z(H)]-[R1_x(H);R1_y(H);R1_z(H)]); 
    % acceleration of M2
    acc2(1:3,H)=(-1)^(2)*G*M(2)/(norm([R2_x(H);R2_y(H);R2_z(H)]'-[R1_x(H);R1_y(H);R1_z(H)]'))^3*([R2_x(H);R2_y(H);R2_z(H)]-[R1_x(H);R1_y(H);R1_z(H)]); 
end
%% velocity and acceleration magnituides
for h=1:length(Rc(1,:))
        MagV1(1,h)=norm([V1_x;V1_y;V1_z]);  % velocity magnituide of M1
        MagV2(1,h)=norm([V2_x;V2_y;V2_z]);  % velocity magnituide of M2
        MagA1(1,h)=norm(acc1(1:3,h));  % acceleration magnituide of M1
        MagA2(1,h)=norm(acc2(1:3,h));  % acceleration magnituide of M2
end
%% M2 (satellite) inertia 
% inertia tensor about C.G
IG=[Ixx, Ixy, Ixz; Ixy, Iyy, Iyz; Ixz, Iyz,Izz];
% inertia tensor about its axes
I0= Ip ( IG, [0, 0, 0], [xG2, yG2, zG2], M(2), 'P' );

% moment due to gravitational force of M2 (satellite)
for k=1:length(Rc(1,:))
    MF2(:,k)=cross([xG2, yG2, zG2],[M(2)*acc2(:,k)]')';
end
% angular acceleration (alpha) and angular velocity (w)
alpha(:,1)=[0;0;0];
alpha(:,2)=[0;0;0];
x0=[xG2 yG2 zG2];
for m=1:length(Rc(1,:))
    A1=MF2(:,m)-I0*alpha(:,m);
    save('I.mat', 'Ixx', 'Iyy', 'Izz', 'Ixy', 'Ixz', 'Iyz', 'A1');
    w(:,m)=[fsolve(@GetW,x0)]';
    clc;
    if m >=2
        alpha(:,m+1)=(w(:,m)-w(:,m-1))/dt;
    end
end
% Get phi, theta and psi of M2 (satellite) about its initial axes
for PP=1:length(w(1,:))
    PHI(PP)=trapz(w(1,1:PP))*180/pi;
    while abs(PHI(PP)) >= 360
        if PHI(PP)>0
            PHI(PP)=PHI(PP)-360;
        elseif PHI(PP)<0
            PHI(PP)=PHI(PP)+360;
        end
    end
    THETA(PP)=trapz(w(2,1:PP))*180/pi;
    while abs(THETA(PP)) >= 360
        if THETA(PP)>0
            THETA(PP)=THETA(PP)-360;
        elseif THETA(PP)<0
            THETA(PP)=THETA(PP)+360;
        end
    end
    PSI(PP)=trapz(w(3,1:PP))*180/pi;
    while abs(PSI(PP)) >= 360
        if PSI(PP)>0
            PSI(PP)=PSI(PP)-360;
        elseif PSI(PP)<0
            PSI(PP)=PSI(PP)+360;
        end
    end
end
%% unit vectors
% ut
for t=1:length(V2_x)
    utx2(t)=V2_x(t)/norm([V2_x(t) V2_y(t) V2_z(t)]);
    uty2(t)=V2_y(t)/norm([V2_x(t) V2_y(t) V2_z(t)]);
    utz2(t)=V2_z(t)/norm([V2_x(t) V2_y(t) V2_z(t)]);
end
% ub
for t=1:length(V2_x)
    VcrossAx2(t)=V2_y(t)*acc1(3,t)-V2_z(t)*acc1(2,t);
    VcrossAy2(t)=V2_z(t)*acc1(1,t)-V2_x(t)*acc1(3,t);
    VcrossAz2(t)=V2_x(t)*acc1(2,t)-V2_y(t)*acc1(1,t);
    MagVcrossA2(t)=sqrt(VcrossAx2(t).^2+VcrossAy2(t).^2+VcrossAz2(t).^2);
    ubx2(t)=VcrossAx2(t)/MagVcrossA2(t);
    uby2(t)=VcrossAy2(t)/MagVcrossA2(t);
    ubz2(t)=VcrossAz2(t)/MagVcrossA2(t);
end
% un
for t=1:length(V2_x)
    unx2(t)=uby2(t)*utz2(t)-ubz2(t)*uty2(t);
    uny2(t)=ubz2(t)*utx2(t)-ubx2(t)*utz2(t);
    unz2(t)=ubx2(t)*uty2(t)-uby2(t)*utx2(t);
end
%% Equatorial plane
    theta_vec=linspace(0,2*pi,30);
    r_vec=0:1e6:1.8e7;
    [theta_mat, r_mat]=meshgrid(theta_vec, r_vec);
    [x_mat, y_mat]=pol2cart(theta_mat, r_mat);
    z_mat=r_mat*0;
%% satellite tracking
figure(1);
% Texturemap the globe
% Load Earth image for texture map
cdata = imread(image_file);
J = imrotate(cdata,0,'bilinear');
V=1;
%--------------------------------------------------------------------------------------------------------------------------------------------------------
for p=1:step_sim:length(V2_y)
%     figure(V);
    set(0,'defaultfigureposition',[40 55 1300 600])
    cla;
    % Options
    space_color = 'k';
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
    % Earth texture image
    % Anything imread() will handle, but needs to be a 2:1 unprojected globe
    % Mean spherical earth
    erad    = 6371008.7714; % equatorial radius (meters)
    prad    = 6371008.7714; % polar radius (meters)
    %GMST0 = []; % Don't set up rotatable globe (ECEF)
    GMST0 = (4.89496121282306 - (p-1)*erot*dt); % Set up a rotatable globe at J2000.0
    set(gcf,'Color','w');
    % M2 trajectory
    plot3(I(1)*ones(1,length(R1_x))+R2_x-R1_x,I(2)*ones(1,length(R1_y))+R2_y-R1_y,I(3)*ones(1,length(R1_z))+R2_z-R1_z,'g','LineWidth',2);
    hold all;
    grid on;
    xlabel('X','Fontsize',18);
    ylabel('Y','Fontsize',18);
    zlabel('Z','Fontsize',18);
    title('Satellite about Earth','Fontsize',18);
    xlim auto;
    ylim auto;
    zlim auto;
    view(-25, 38);
    % M2 start
    [ model_handle ] = WireCubes( [R2_x(1,p)-R1_x(1,p),R2_y(1,p)-R1_y(1,p),R2_z(1,p)-R1_z(1,p)],1e6,[.9,0.2,0.5],1 );
    % Create wireframe globe
    % Create a 3D meshgrid of the sphere points using the ellipsoid function
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    globe = surf(y+R1_x(1,p), -x+R1_y(1,p), -z+R1_z(1,p), 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    if ~isempty(GMST0)
        hgx = hgtransform;
        set(hgx,'Matrix', makehgtform('zrotate',GMST0));
        set(globe,'Parent',hgx);
    end
    % Set image as color data (cdata) property, and set face color to indicate
    % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    %Projected path of the satalite on the earth
    plot3(XX1,YY1,ZZ1,'r','LineWidth',2);
    % Equatorial plane
    if (R2_z(1,p)-R1_z(1,p)) >=0
        color='yellow';
    else
        color=[0.5,0.5,0.5];
    end
    surf(x_mat, y_mat, z_mat, 'FaceAlpha', 0.3, 'EdgeColor', color, 'FaceColor', color, 'EdgeAlpha', 0.2);
    
    vectarrow([R2_x(p) R2_y(p) R2_z(p)],[R2_x(p) R2_y(p) R2_z(p)]+[utx2(p) uty2(p) utz2(p)]*5e6,2,'b')
    vectarrow([R2_x(p) R2_y(p) R2_z(p)],[R2_x(p) R2_y(p) R2_z(p)]+[unx2(p) uny2(p) unz2(p)]*5e6,2,'r')
    vectarrow([R2_x(p) R2_y(p) R2_z(p)],[R2_x(p) R2_y(p) R2_z(p)]+[ubx2(p) uby2(p) ubz2(p)]*5e6,2,[0.3,0.7,0.1])
    
    [ Q ] = RM ( PHI(p), THETA(p), PSI(p), '123');
    VV1=Q*[utx2(p); uty2(p); utz2(p)];
    VV2=Q*[unx2(p); uny2(p); unz2(p)];
    VV3=Q*[ubx2(p); uby2(p); ubz2(p)];
    
    vectarrow([R2_x(p) R2_y(p) R2_z(p)],[R2_x(p) R2_y(p) R2_z(p)]+VV1'*5e6,2,'c')
    vectarrow([R2_x(p) R2_y(p) R2_z(p)],[R2_x(p) R2_y(p) R2_z(p)]+VV2'*5e6,2,'m')
    vectarrow([R2_x(p) R2_y(p) R2_z(p)],[R2_x(p) R2_y(p) R2_z(p)]+VV3'*5e6,2,'k')
    legend('Satellite trajectory','Satellite','Projected path of the satalite on the earth','Equatorial plane','x','y','z','x''','y''','z''','location','northeast');
    vectarrow([0 0 0],[2 0 0]*1e7,2,'b')
    vectarrow([0 0 0],[0 2 0]*1e7,2,'b')
    vectarrow([0 0 0],[0 0 2]*1e7,2,'b')
    
    [ Q ] = RM ((p-1)*erot*dt*180/pi, 0, 0, '313');
    vv1=Q*[2;0;0]*1e7;
    vv2=Q*[0;2;0]*1e7;
    vv3=Q*[0;0;2]*1e7;
    
    vectarrow([0 0 0],vv1',2,'r')
    vectarrow([0 0 0],vv2',2,'r')
    vectarrow([0 0 0],vv3',2,'r')

%     img= getframe(gcf);
%     imwrite(img.cdata, [num2str(V), '.png']);
%     V=V+1;
    pause(1/hz)
%     close;
end
%% 
figure(2);
set(gcf,'color','w')
set(0,'defaultfigureposition',[40 55 1300 600])
subplot(3,1,1);
plot(0:dt:tf+dt,PHI,'b','linewidth',2);
ylabel('\phi (degree)','Fontsize',18);
grid on;
xlim([0, tf+dt]);
title(['Satellite angles change' ' (Ixx_G = '  num2str(Ixx) ', Iyy_G = '  num2str(Iyy) ', Izz_G = '  num2str(Izz) ', C.G = [' num2str(xG2) ', ' num2str(yG2) ', ' num2str(zG2) '])'],'Fontsize',18);
subplot(3,1,2);
plot(0:dt:tf+dt,THETA,'r','linewidth',2);
grid on;
xlim([0, tf+dt]);
ylabel('\theta (degree)','Fontsize',18);
subplot(3,1,3);
plot(0:dt:tf+dt,PSI,'m','linewidth',2);
grid on;
xlim([0, tf+dt]);
 xlabel('Time (sec)','Fontsize',18);
 ylabel('\psi (degree)','Fontsize',18);