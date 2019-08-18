%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                     SHOULDER MODEL 3D                           %%%%%
%%%%%                                                                 %%%%%
%%%%% Author: Humberto J De las Casas	                              %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% All bodies in global coordinates:
%%%%% Axis on Schematic (Anterior View): Left "Z", Up "Y", Out "X".
ground=[0 0 0];
clavicle=[0.006325 0.00693 0.025465];
scapula=[-0.01433 0.02007 0.135535]+clavicle;
humerus=[-0.00955 -0.034 0.009]+scapula;
ulna=[0.0061 -0.2904 -0.0123]+humerus;
radius=[0.0004 -0.011503 0.019999]+ulna;
proximal_row=[0.018 -0.242 0.025]+radius;
hand=[0.003992 -0.015054 0.002327]+proximal_row;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Bodies used for the project: "ground" and "humerus". PLOTS.
[x,y,z]=ellipsoid(ground(1),ground(2)-0.5,ground(3),0.1,0.5,0.12);
figure; surf(x,y,z,'FaceColor','k'); alpha 0.15; axis equal; view(90,0)
[x,y,z]=sphere; hold on; xlabel('X'); ylabel('Y'); zlabel('Z'); camroll(90)
surf(x*0.03,y*0.03,z*0.03,'FaceColor','b')
[x,y,z]=sphere;
surf(x*0.03+humerus(1),y*0.03+humerus(2),z*0.03+humerus(3),'FaceColor','r')
surf(x*0.03+humerus(1),y*0.03+humerus(2)-0.55,z*0.03+humerus(3),...
     'FaceColor','y'); legend('Body','Ground','Humerus','Wrist')
plot3([ground(1) humerus(1)],[ground(2) humerus(2)],...
      [ground(3) humerus(3)],'k','LineWidth',10)
plot3([humerus(1) humerus(1)],[humerus(2) humerus(2)-0.55],...
      [humerus(3) humerus(3)],'k','LineWidth',10)
plot3([ground(1) 0.30],[ground(2) 0],[ground(3) 0],'r','LineWidth',2)
plot3([ground(1) 0],[ground(2) 0.30],[ground(3) 0],'y','LineWidth',2)
plot3([ground(1) 0],[ground(2) 0],[ground(3) 0.30],'g','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Muscles Attachments (10) - INITIAL POSITIONS.
DELT1_A=[0.009 -0.119 0.006]+humerus;           % Deltoid - 1.
DELT1_B=[-0.014 0.011 0.08]+clavicle;
DELT2_A=[0.005 -0.136 0.006]+humerus;           % Deltoid - 2.
DELT2_B=[-0.011 0 0.006]+scapula;
DELT3_A=[-0.056 0.001 -0.025]+scapula;          % Deltoid - 3.
DELT3_B=[0.002 -0.076 0.01]+humerus;
SUPSP_A=[0.003 0.011 0.026]+humerus;            % Supraspinatus.
SUPSP_B=[-0.044 -0.015 -0.059]+scapula;
INFSP_A=[-0.009 0.005 0.024]+humerus;           % Infraspinatus.
INFSP_B=[-0.074 -0.055 -0.048]+scapula;
SUBSC_A=[0.014 0.008 -0.013]+humerus;           % Subscapularis.
SUBSC_B=[-0.072 -0.039 -0.065]+scapula;
TMIN_A=[-0.001 -0.013 0.022]+humerus;           % Teres Minor.
TMIN_B=[-0.096 -0.081 -0.053]+scapula;
TMAJ_A=[0.01 -0.054 -0.006]+humerus;            % Teres Major.
TMAJ_B=[-0.105 -0.103 -0.058]+scapula;
CORB_A=[0.012 -0.041 -0.027]+scapula;           % Coracobra Chialis.
CORB_B=[0.007 -0.15 -0.008]+humerus;
TRIlong_A=[-0.046 -0.041 -0.014]+scapula;       % Triceps Brachi Long Head.
TRIlong_B=[-0.022 0.01 -0.001]+ulna;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Ploting the Muscles Attachments.
plot3([DELT1_A(1) DELT1_B(1)],[DELT1_A(2) DELT1_B(2)],...
      [DELT1_A(3) DELT1_B(3)],'c','LineWidth',2);
plot3([DELT2_A(1) DELT2_B(1)],[DELT2_A(2) DELT2_B(2)],...
      [DELT2_A(3) DELT2_B(3)],'c','LineWidth',2);
plot3([DELT3_A(1) DELT3_B(1)],[DELT3_A(2) DELT3_B(2)],...
      [DELT3_A(3) DELT3_B(3)],'c','LineWidth',2);
plot3([SUPSP_A(1) SUPSP_B(1)],[SUPSP_A(2) SUPSP_B(2)],...
      [SUPSP_A(3) SUPSP_B(3)],'c','LineWidth',2);
plot3([INFSP_A(1) INFSP_B(1)],[INFSP_A(2) INFSP_B(2)],...
      [INFSP_A(3) INFSP_B(3)],'c','LineWidth',2);
plot3([SUBSC_A(1) SUBSC_B(1)],[SUBSC_A(2) SUBSC_B(2)],...
      [SUBSC_A(3) SUBSC_B(3)],'c','LineWidth',2);
plot3([TMIN_A(1) TMIN_B(1)],[TMIN_A(2) TMIN_B(2)],...
      [TMIN_A(3) TMIN_B(3)],'c','LineWidth',2);
plot3([TMAJ_A(1) TMAJ_B(1)],[TMAJ_A(2) TMAJ_B(2)],...
      [TMAJ_A(3) TMAJ_B(3)],'c','LineWidth',2);
plot3([CORB_A(1) CORB_B(1)],[CORB_A(2) CORB_B(2)],...
      [CORB_A(3) CORB_B(3)],'c','LineWidth',2);
plot3([TRIlong_A(1) TRIlong_B(1)],[TRIlong_A(2) TRIlong_B(2)],...
      [TRIlong_A(3) TRIlong_B(3)],'c','LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Symbolic Variables. 
syms q1 q2 q3 q1dot q2dot q3dot q1ddot q2ddot q3ddot g
qdot=[q1dot;q2dot;q3dot];qddot=[q1ddot;q2ddot;q3ddot];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Frames transformation on function of q=[q1 q2 q3].  
Q_00=eye(4)+[zeros(4,3) [ground(1);ground(2);ground(3);0]];
Q_Humerus_0=(eye(4)+[zeros(4,3) [humerus(1);humerus(2);humerus(3);0]]);
Q_Humerus_1=Q_Humerus_0*[cos(q1) -sin(q1) 0 0;sin(q1) cos(q1) 0 0 ...
        ;0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1];
Q_Humerus_2=Q_Humerus_1*...
        [cos(q2) -sin(q2) 0 0;sin(q2) cos(q2) ...
        0 0;0 0 1 0;0 0 0 1]*[0 -1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1]*...
        [1 0 0 0;0 0 -1 0;0 1 0 0; 0 0 0 1];
Q_Humerus_3=Q_Humerus_2*[cos(q3) -sin(q3) ...
        0 0;sin(q3) cos(q3) 0 0;0 0 1 0;0 0 0 1];

Q_Ulna_1=Q_Humerus_1*(eye(4)+[zeros(4,3) [ulna(1)-humerus(1);...
        -ulna(3)+humerus(3);ulna(2)-humerus(2);0]]);
Q_Ulna_2=Q_Humerus_2*(eye(4)+[zeros(4,3) [-ulna(3)+humerus(3);...
        ulna(2)-humerus(2);ulna(1)-humerus(1);0]]);
Q_Ulna_3=Q_Humerus_3*(eye(4)+[zeros(4,3) [-ulna(3)+humerus(3);...
        ulna(2)-humerus(2);ulna(1)-humerus(1);0]]);

Q_Radius_1=Q_Ulna_1*(eye(4)+[zeros(4,3) [radius(1)-ulna(1);...
        -radius(3)+ulna(3);radius(2)-ulna(2);0]]);
Q_Radius_2=Q_Ulna_2*(eye(4)+[zeros(4,3) [-radius(3)+ulna(3);...
        radius(2)-ulna(2);radius(1)-ulna(1);0]]);
Q_Radius_3=Q_Ulna_3*(eye(4)+[zeros(4,3) [-radius(3)+ulna(3);...
        radius(2)-ulna(2);radius(1)-ulna(1);0]]);

Q_PR_1=Q_Radius_1*(eye(4)+[zeros(4,3) [proximal_row(1)-radius(1);...
        -proximal_row(3)+radius(3);proximal_row(2)-radius(2);0]]);
Q_PR_2=Q_Radius_2*(eye(4)+[zeros(4,3) [-proximal_row(3)+radius(3);...
        proximal_row(2)-radius(2);proximal_row(1)-radius(1);0]]);
Q_PR_3=Q_Radius_3*(eye(4)+[zeros(4,3) [-proximal_row(3)+radius(3);...
        proximal_row(2)-radius(2);proximal_row(1)-radius(1);0]]);

Q_Hand_1=Q_PR_1*(eye(4)+[zeros(4,3) [hand(1)-proximal_row(1);...
        -hand(3)+proximal_row(3);hand(2)-proximal_row(2);0]]);
Q_Hand_2=Q_PR_2*(eye(4)+[zeros(4,3) [-hand(3)+proximal_row(3);...
        hand(2)-proximal_row(2);hand(1)-proximal_row(1);0]]);
Q_Hand_3=Q_PR_3*(eye(4)+[zeros(4,3) [-hand(3)+proximal_row(3);...
        hand(2)-proximal_row(2);hand(1)-proximal_row(1);0]]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian respect to the center of mass.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Physical parameters - Humerus:
Humerus_mass=1.99757;
Humerus_mass_center_1=Q_Humerus_1*(eye(4)+[zeros(4,3) ...
        [0.018064;0.012746;-0.140141;0]]);
Humerus_mass_center_2=Q_Humerus_2*(eye(4)+[zeros(4,3) ...
        [0.012746;-0.140141;0.018064;0]]);    
Humerus_mass_center_3=Q_Humerus_3*(eye(4)+[zeros(4,3) ...
        [0.012746;-0.140141;0.018064;0]]);
Humerus_I=[0.01227763 -3.4741E-4 -2.325E-4;-3.4741E-4 0.00255133 0.0012293; ...
        -2.325E-4 0.0012293 0.01257888];                                   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Physical parameters - Ulna, Radius and Hand:
Ulna_mass=1.1053;
Radius_mass=0.23359;
Hand_mass=0.5819;
Ulna_mass_center_1=Q_Ulna_1*(eye(4)+[zeros(4,3) ...
        [0.00971783;-0.024286;-0.0959509;0]]);
Ulna_mass_center_2=Q_Ulna_2*(eye(4)+[zeros(4,3) ...
        [-0.024286;-0.0959509;0.00971783;0]]);
Ulna_mass_center_3=Q_Ulna_3*(eye(4)+[zeros(4,3) ...
        [-0.024286;-0.0959509;0.00971783;0]]);
Radius_mass_center_1=Q_Radius_1*(eye(4)+[zeros(4,3) ...
        [0.0336341;-0.0156;-0.181559;0]]);
Radius_mass_center_2=Q_Radius_2*(eye(4)+[zeros(4,3) ...
        [-0.0156;-0.181559;0.0336341;0]]);
Radius_mass_center_3=Q_Radius_3*(eye(4)+[zeros(4,3) ...
        [-0.0156;-0.181559;0.0336341;0]]);
Hand_mass_center_1=Q_Hand_1*(eye(4)+[zeros(4,3) ...
        [-0.00301314;0.00112205;-0.0424993;0]]);
Hand_mass_center_2=Q_Hand_2*(eye(4)+[zeros(4,3) ...
        [0.00112205;-0.0424993;-0.00301314;0]]);
Hand_mass_center_3=Q_Hand_3*(eye(4)+[zeros(4,3) ...
        [0.00112205;-0.0424993;-0.00301314;0]]);
Ulna_I=[0.00541309 3.1686E-4 -7.615E-5;3.1686E-4 ...
        0.00115318 0.00109169;-7.615E-5 0.00109169 0.00494361];
Radius_I=[4.3855E-4 3.014E-5 -4.24E-6; 3.014E-5 8.859E-5 ...
        6.418E-5; -4.24E-6 6.418E-5 4.0258E-4];
Hand_I=[1.1E-4 9.0E-7 -2.0E-7;9.0E-7 6.0E-5 1.2E-5; ...
        -2.0E-7 1.2E-5 1.5E-4];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Due to the model only has 3 DOFs on the humerus socket-ball joint, 
%%%%% the Z vectors direction are the same for all frames.
Z_0=Q_Humerus_0(1:3,3);
Z_1=Q_Humerus_1(1:3,3);
Z_2=Q_Humerus_2(1:3,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian - Humerus:
J_Humerus=cross(Z_0,Humerus_mass_center_1(1:3,4)-Q_Humerus_0(1:3,4));
J_Humerus(1:3,2)=cross(Z_1,Humerus_mass_center_2(1:3,4)-Q_Humerus_0(1:3,4));
J_Humerus(1:3,3)=cross(Z_2,Humerus_mass_center_3(1:3,4)-Q_Humerus_0(1:3,4));    
J_Humerus(4:6,1:3)=[Z_0 Z_1 Z_2];
J_Humerus=simplify(J_Humerus);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian - Ulna:
J_Ulna=cross(Z_0,Ulna_mass_center_1(1:3,4)-Q_Humerus_0(1:3,4));
J_Ulna(1:3,2)=cross(Z_1,Ulna_mass_center_2(1:3,4)-Q_Humerus_0(1:3,4));
J_Ulna(1:3,3)=cross(Z_2,Ulna_mass_center_3(1:3,4)-Q_Humerus_0(1:3,4));    
J_Ulna(4:6,1:3)=[Z_0 Z_1 Z_2];
J_Ulna=simplify(J_Ulna);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian - Radius:
J_Radius=cross(Z_0,Radius_mass_center_1(1:3,4)-Q_Humerus_0(1:3,4));
J_Radius(1:3,2)=cross(Z_1,Radius_mass_center_2(1:3,4)-Q_Humerus_0(1:3,4));
J_Radius(1:3,3)=cross(Z_2,Radius_mass_center_3(1:3,4)-Q_Humerus_0(1:3,4));    
J_Radius(4:6,1:3)=[Z_0 Z_1 Z_2];
J_Radius=simplify(J_Radius);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian - Hand:
J_Hand=cross(Z_0,Hand_mass_center_1(1:3,4)-Q_Humerus_0(1:3,4));
J_Hand(1:3,2)=cross(Z_1,Hand_mass_center_2(1:3,4)-Q_Humerus_0(1:3,4));
J_Hand(1:3,3)=cross(Z_2,Hand_mass_center_3(1:3,4)-Q_Humerus_0(1:3,4));    
J_Hand(4:6,1:3)=[Z_0 Z_1 Z_2];
J_Hand=simplify(J_Hand);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Mass Matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Mass - Humerus:
R1=Q_Humerus_1(1:3,1:3);
R2=Q_Humerus_2(1:3,1:3);
R3=Q_Humerus_3(1:3,1:3);
M_Humerus=simplify(...
            J_Humerus(4:6,1:3).'*R1*Humerus_I*R1.'*J_Humerus(4:6,1:3)+...
            J_Humerus(4:6,1:3).'*R2*Humerus_I*R2.'*J_Humerus(4:6,1:3)+...
            J_Humerus(4:6,1:3).'*R3*Humerus_I*R3.'*J_Humerus(4:6,1:3)+...
            Humerus_mass*J_Humerus(1:3,1:3).'*J_Humerus(1:3,1:3));
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Mass - Ulna:
R1=Q_Ulna_1(1:3,1:3);
R2=Q_Ulna_2(1:3,1:3);
R3=Q_Ulna_3(1:3,1:3);
M_Ulna=simplify(...
            J_Ulna(4:6,1:3).'*R1*Ulna_I*R1.'*J_Ulna(4:6,1:3)+...
            J_Ulna(4:6,1:3).'*R2*Ulna_I*R2.'*J_Ulna(4:6,1:3)+...
            J_Ulna(4:6,1:3).'*R3*Ulna_I*R3.'*J_Ulna(4:6,1:3)+...
            Ulna_mass*J_Ulna(1:3,1:3).'*J_Ulna(1:3,1:3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Mass - Radius:
R1=Q_Radius_1(1:3,1:3);
R2=Q_Radius_2(1:3,1:3);
R3=Q_Radius_3(1:3,1:3);
M_Radius=simplify(...
            J_Radius(4:6,1:3).'*R1*Radius_I*R1.'*J_Radius(4:6,1:3)+...
            J_Radius(4:6,1:3).'*R2*Radius_I*R2.'*J_Radius(4:6,1:3)+...
            J_Radius(4:6,1:3).'*R3*Radius_I*R3.'*J_Radius(4:6,1:3)+...
            Radius_mass*J_Radius(1:3,1:3).'*J_Radius(1:3,1:3));
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Mass - Hand:
R1=Q_Hand_1(1:3,1:3);
R2=Q_Hand_2(1:3,1:3);
R3=Q_Hand_3(1:3,1:3);
M_Hand=simplify(...
            J_Hand(4:6,1:3).'*R1*Hand_I*R1.'*J_Hand(4:6,1:3)+...
            J_Hand(4:6,1:3).'*R2*Hand_I*R2.'*J_Hand(4:6,1:3)+...
            J_Hand(4:6,1:3).'*R3*Hand_I*R3.'*J_Hand(4:6,1:3)+...
            Hand_mass*J_Hand(1:3,1:3).'*J_Hand(1:3,1:3));
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Mass Matrix.
M_q=M_Humerus+M_Ulna+M_Radius+M_Hand;
isequal(M_q,M_q.');                     % It is symmetric!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Coriolis Matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
q=[q1;q2;q3];
for i=1:3
    for j=1:3
        for k=1:3
c(i,j,k)=((diff(M_q(k,j),q(i)))+(diff(M_q(k,i),q(j)))...
        -(diff(M_q(i,j),q(k))))/2;
        end
    end
end
for k=1:3
    for j=1:3
C_qqdot(k,j)=c(1,j,k)*q1dot+c(2,j,k)*q2dot+c(3,j,k)*q3dot;
    end
end
for i=1:3
    for j=1:3
M_DOT(i,j)=diff(M_q(i,j),q1)*q1dot+diff(M_q(i,j),q2)*...
           q2dot+diff(M_q(i,j),q3)*q3dot;
    end
end
N_qqdot=simplify(M_DOT-2*C_qqdot);
isequal(N_qqdot,-N_qqdot.');                % It is skew-symmetric!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Gravity Matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
g_v=[0;0;-g];
P_Humerus=-Humerus_mass*g_v.'*Humerus_mass_center_3(1:3,4);
P_Ulna=-Ulna_mass*g_v.'*Ulna_mass_center_3(1:3,4);
P_Radius=-Radius_mass*g_v.'*Radius_mass_center_3(1:3,4);
P_Hand=-Hand_mass*g_v.'*Hand_mass_center_3(1:3,4);
P_Total=P_Humerus+P_Ulna+P_Radius+P_Hand;

g_q(1)=simplify(diff(P_Total,q1));
g_q(2,1)=simplify(diff(P_Total,q2));
g_q(3,1)=simplify(diff(P_Total,q3));   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PART II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Muscles Attachments (10) - IN FCN OF q=[q1 q2 q3].
DELT1_A=Q_Humerus_3*(eye(4)+[zeros(4,3) [0.009;-0.119;0.006;0]]);
DELT1_B=[-0.014 0.011 0.08]+clavicle;
DELT2_A=Q_Humerus_3*(eye(4)+[zeros(4,3) [0.005;-0.136;0.006;0]]);
DELT2_B=[-0.011 0 0.006]+scapula;
DELT3_A=[-0.056 0.001 -0.025]+scapula;
DELT3_B=Q_Humerus_3*(eye(4)+[zeros(4,3) [0.002;-0.076;0.01;0]]);
SUPSP_A=Q_Humerus_3*(eye(4)+[zeros(4,3) [0.003;0.011;0.026;0]]);
SUPSP_B=[-0.044 -0.015 -0.059]+scapula;
INFSP_A=Q_Humerus_3*(eye(4)+[zeros(4,3) [-0.009;0.005;0.024;0]]);
INFSP_B=[-0.074 -0.055 -0.048]+scapula;
SUBSC_A=Q_Humerus_3*(eye(4)+[zeros(4,3) [0.014;0.008;-0.013;0]]);
SUBSC_B=[-0.072 -0.039 -0.065]+scapula;
TMIN_A=Q_Humerus_3*(eye(4)+[zeros(4,3) [-0.001;-0.013;0.022;0]]);
TMIN_B=[-0.096 -0.081 -0.053]+scapula;
TMAJ_A=Q_Humerus_3*(eye(4)+[zeros(4,3) [0.01;-0.054;-0.006;0]]);
TMAJ_B=[-0.105 -0.103 -0.058]+scapula;
CORB_A=[0.012 -0.041 -0.027]+scapula;
CORB_B=Q_Humerus_3*(eye(4)+[zeros(4,3) [0.007;-0.15;-0.008;0]]);
TRIlong_A=[-0.046 -0.041 -0.014]+scapula;
TRIlong_B=Q_Ulna_3*(eye(4)+[zeros(4,3) [-0.022;0.01;-0.001;0]]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Computing muscle lengths:
L_DELT1=simplify(norm(DELT1_A(1:3,4)-DELT1_B.'));
L_DELT2=simplify(norm(DELT2_A(1:3,4)-DELT2_B.'));
L_DELT3=simplify(norm(DELT3_A.'-DELT3_B(1:3,4)));
L_SUPSP=simplify(norm(SUPSP_A(1:3,4)-SUPSP_B.'));
L_INFSP=simplify(norm(INFSP_A(1:3,4)-INFSP_B.'));
L_SUBSC=simplify(norm(SUBSC_A(1:3,4)-SUBSC_B.'));
L_TMIN=simplify(norm(TMIN_A(1:3,4)-TMIN_B.'));
L_TMAJ=simplify(norm(TMAJ_A(1:3,4)-TMAJ_B.'));
L_CORB=simplify(norm(CORB_A.'-CORB_B(1:3,4)));
L_TRIlong=simplify(norm(TRIlong_A.'-TRIlong_B(1:3,4)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Computing unit vector of muscles:
Unit_DELT1=simplify((DELT1_A(1:3,4)-DELT1_B.')/L_DELT1);
Unit_DELT2=simplify((DELT2_A(1:3,4)-DELT2_B.')/L_DELT2);
Unit_DELT3=simplify((DELT3_A.'-DELT3_B(1:3,4))/L_DELT3);
Unit_SUPSP=simplify((SUPSP_A(1:3,4)-SUPSP_B.')/L_SUPSP);
Unit_INFSP=simplify((INFSP_A(1:3,4)-INFSP_B.')/L_INFSP);
Unit_SUBSC=simplify((SUBSC_A(1:3,4)-SUBSC_B.')/L_SUBSC);
Unit_TMIN=simplify((TMIN_A(1:3,4)-TMIN_B.')/L_TMIN);
Unit_TMAJ=simplify((TMAJ_A(1:3,4)-TMAJ_B.')/L_TMAJ);
Unit_CORB=simplify((CORB_A.'-CORB_B(1:3,4))/L_CORB);
Unit_TRIlong=simplify((+TRIlong_A.'-TRIlong_B(1:3,4))/L_TRIlong);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Optimal lengths:
L_opt_DELT1=0.0976;
L_opt_DELT2=0.1078;
L_opt_DELT3=0.1367;
L_opt_SUPSP=0.0682;
L_opt_INFSP=0.0755;
L_opt_SUBSC=0.0873;
L_opt_TMIN=0.0741;
L_opt_TMAJ=0.1624;
L_opt_CORB=0.0932;
L_opt_TRIlong=0.134;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initial lengths:
q1=0;q2=0;q3=0;
L_0_DELT1=eval(L_DELT1);                                                        
L_0_DELT2=eval(L_DELT2);
L_0_DELT3=eval(L_DELT3);
L_0_SUPSP=eval(L_SUPSP);
L_0_INFSP=eval(L_INFSP);
L_0_SUBSC=eval(L_SUBSC);
L_0_TMIN=eval(L_TMIN);
L_0_TMAJ=eval(L_TMAJ);
L_0_CORB=eval(L_CORB);
L_0_TRIlong=eval(L_TRIlong);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Split percentage between segment and tendon lengths:
SP_DELT1=L_0_DELT1;
SP_DELT2=L_0_DELT2;
SP_DELT3=L_0_DELT3;
SP_SUPSP=L_0_SUPSP;
SP_INFSP=L_0_INFSP;
SP_SUBSC=L_0_SUBSC;
SP_TMIN=L_0_TMIN;
SP_TMAJ=L_0_TMAJ;
SP_CORB=L_0_CORB;
SP_TRIlong=L_0_TRIlong;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Dynamic Equation for initial values of q=[0 0 0].
DE_0=simplify(eval(M_q*qddot+C_qqdot*qdot+g_q));






