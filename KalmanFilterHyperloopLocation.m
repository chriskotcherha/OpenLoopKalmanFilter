function [predictedState, predictedCovariance] = KalmanFilterHyperloop( prevState, prevCovariance, IMUData, sensorData, execution )

sensorPositions = zeros(7,3);
thicknessOfRail = 1;


p1 = sensorPositions(1,:)'; %Position of the sensor with respect to the center of mass of the pod
p2 = sensorPositions(2,:)'; %Position of the sensor with respect to the center of mass of the pod
p3 = sensorPositions(3,:)'; %Position of the sensor with respect to the center of mass of the pod
tck = thicknessOfRail/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing state, kalman gain etc%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaT=1/4000; %One time step

%Currently using IMU as control to make the state prediction and so
%initializing that

%%%%%%%%%%%%%%%% ----------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xkk=prevState; %This is the previous state and in the case of k=1 the initial state
uk=IMUData; %Control vector, currently IMU being used
Pkk=prevCovariance;

%
axL=uk(1);  %accelaration in x from IMU
ayL=uk(2);  %acceleration in y from IMU
azL=uk(3);  %acceleration in z from IMU
omegaX=uk(4);   %angular velocity in x from IMU
omegaY=uk(5);   %angular velocity in y from IMU
omegaZ=uk(6);   %angular velocity in z from IMU
%

%quaternions from the state
q1=xkk(7);  
q2=xkk(8);
q3=xkk(9);
q0=xkk(10);

%rotation matrix from the quaternions
Rot=[1-2*q2^2-2*q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2);...
     2*(q1*q2+q0*q3) 1-2*q1^2-2*q3^2 2*(q2*q3-q0*q1);...
     2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) 1-2*q1^2-2*q2^2;];


 %Fk matrix used for predicting the next state and the error in the
 %next state
 Fk = [1 0 0 deltaT 0 0 0 0 0 0;...
       0 1 0 0 deltaT 0 0 0 0 0;...
       0 0 1 0 0 deltaT 0 0 0 0;...
       0 0 0 1 0 0 0 0 0 0;...
       0 0 0 0 1 0 0 0 0 0;...
       0 0 0 0 0 1 0 0 0 0;...
       0 0 0 0 0 0 1 0 0 0;...
       0 0 0 0 0 0 0 1 0 0;...
       0 0 0 0 0 0 0 0 1 0;...
       0 0 0 0 0 0 0 0 0 1;];



 %Bk matrix, multiplied to the control vector (currently IMU) in order
 %to also predict the next state of the pod
 Bk=[ 0 0 0 0 0 0;...
      0 0 0 0 0 0;...
      0 0 0 0 0 0;...
      Rot(1,1) Rot(1,2) Rot(1,3) 0 0 0;...
      Rot(2,1) Rot(2,2) Rot(2,3) 0 0 0;...
      Rot(3,1) Rot(3,2) (Rot(3,3)-g/(azL)) 0 0 0;...
      0 0 0 q0/2 -q3/2 q2/2;...
      0 0 0 q3/2 q0/2 -q1/2;...
      0 0 0 -q2/2 q1/2 q0/2;...
      0 0 0 -q1/2 -q2/2 -q3/2;];

%This is the portion for Wk and this is the same as patricks E80 code for rn. We will be expirimentally determining this next semester%

%OLDuk=Control(:,k);
OmegaMatrix=[0 omegaZ -omegaY omegaX;...
            -omegaZ 0 omegaX omegaY;...
            omegaY -omegaX 0 omegaZ;...
            -omegaX -omegaY -omegaZ 0;];
% dqdt=(1/2).*(OmegaMatrix*xkk(7:10));
% dq1dt=dqdt(1);
% dq2dt=dqdt(2);
% dq3dt=dqdt(3);
% dq0dt=dqdt(4);
% domegadt=(NEXTuk(4:6)-OLDuk(4:6))./(2*deltaT);
% domegaXdt=domegadt(1);
% domegaYdt=domegadt(2);
% domegaZdt=domegadt(3);


% wk=(0.5*deltaT^2).*([(Rot*[axL;ayL;azL;])' ([1-4*q2*dq2dt-4*q3*dq3dt 2*(q1*dq2dt+q2*dq1dt-q0*dq3dt-q3*dq0dt) 2*(q1*dq3dt+q3*dq1dt+q0*dq2dt+q2*dq0dt);... %([1-4*q2*dq2dt-4*q3*dq3dt 2*(q1*dq2dt+q2*dq1dt+q0*dq3dt+q3*dq0dt) 2*(-q0*dq2dt-q2*dq0dt+q1*dq3dt+q3*dq1dt);%2*(-q0*dq3dt-q3*dq0dt+q1*dq2dt+q2*dq1dt) 1-4*q1*dq1dt-4*q3*dq3dt 2*(q2*dq3dt+q3*dq2dt+q0*dq1dt+q1*dq0dt); 2*(q1*dq3dt+q3*dq1dt+q0*dq2dt+q2*dq0dt) 2*(-q0*dq1dt-q1*dq0dt+q2*dq3dt+q3*dq2dt) 1-4*q2*dq2dt-4*q1*dq1dt;]*[axL;ayL;azL;]+Rot*([NEXTaxL-OLDaxL;NEXTayL-OLDayL;NEXTazL-OLDazL;]./(2*deltaT)))'  
%                  2*(q1*dq2dt+q2*dq1dt+q0*dq3dt+q3*dq0dt) 1-4*q1*dq1dt-4*q3*dq3dt 2*(q2*dq3dt+q3*dq2dt-q0*dq1dt-q1*dq0dt);...
%                  2*(q1*dq3dt+q3*dq1dt-q0*dq2dt-q2*dq0dt) 2*(q2*dq3dt+q3*dq2dt+q0*dq1dt+q1*dq0dt) 1-4*q2*dq2dt-4*q1*dq1dt;]*[axL;ayL;azL;]+Rot*([NEXTaxL-OLDaxL;NEXTayL-OLDayL;NEXTazL-OLDazL;]./(2*deltaT)))' ((1/2)*(OmegaMatrix*dqdt+...
%       [0 domegaZdt -domegaYdt domegaXdt;...
%       -domegaZdt 0 domegaXdt domegaYdt;...
%       domegaYdt -domegaXdt 0 domegaZdt;...
%       -domegaXdt -domegaYdt -domegaZdt 0;]*xkk(7:10)))'
%     ]').^2;
% Wk=(wk*wk')*(k~=1)+diag([wk])*(k==1);%diag(wk);
%Experimental determination bit ends%



%The derivative of the rotation matrix in terms of quaternions










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREDICTION STEP

xkp1k=xkk+deltaT*[  xkk(4:6);...
                       (Rot*uk(1:3))-[0;0;9.81;];...
                       (1/2).*(OmegaMatrix*xkk(7:10));];  %prediction step of the state 

Pkp1k=Fk*Pkk*Fk'+Bk*Qk*Bk'; %prediction step of the error, needs to be experimentally determined 


normQuat=sqrt(sum((xkp1k(7:10)).^2));
xkp1k(7:10)=xkp1k(7:10)./normQuat;
normQuat=sqrt(sum((xkp1k(7:10)).^2));


%%%%%%%%%%%%%%%%%%%----------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRECTIVE STEPS

% CORRECTIVE STEP 1, XZ LASER SCANNER
if execution(1)==1
        
    z = sensorData(1,:);
    
    rz1=xkp1k(3);
    q11=xkp1k(7);  
    q12=xkp1k(8);
    q13=xkp1k(9);
    q10=xkp1k(10);
    Rot1=[1-2*q12^2-2*q13^2 2*(q11*q12-q10*q13) 2*(q11*q13+q10*q12);...
     2*(q11*q12+q10*q13) 1-2*q11^2-2*q13^2 2*(q12*q13-q10*q11);...
     2*(q11*q13-q10*q12) 2*(q12*q13+q10*q11) 1-2*q11^2-2*q12^2;];
    sz1=(Rot1(3,:)*p1+rz1);
    sq1=sqrt(Rot1(2,2)^+Rot1(1,2)^2);
    
    dd1dq1=(sq1*([2*q13 2*q10 -4*q11]*p1)-sz1*(Rot1(2,2)*(-4*q11) + Rot1(1,2)*2*q12)/sq1)/(sq1^2);
    dd1dq2=(sq1*([-2*q10 2*q13 -4*q12]*p1)-sz1*(Rot1(2,2)*0 + Rot1(1,2)*2*q11)/sq1)/(sq1^2);
    dd1dq3=(sq1*([2*q11 2*q12 0]*p1)-sz1*(Rot1(2,2)*(-4*q13) + Rot1(1,2)*-2*q10)/sq1)/(sq1^2);
    dd1dq0=(sq1*([-2*q12 2*q11 0]*p1)-sz1*(Rot1(2,2)*0 + Rot1(1,2)*-2*q13)/sq1)/(sq1^2);
    
    dtheta1dq1=(-1/sqrt(1-((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1)^2))*(sq1*((-4*q11)*Rot1(2,1)+2*q12*Rot1(1,2)-0*Rot1(2,2)-2*q12*Rot1(1,1)) - (Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))*(Rot1(2,2)*(-4*q11) + Rot1(1,2)*2*q12)/sq1)/(sq1^2);
    dtheta1dq2=(-1/sqrt(1-((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1)^2))*(sq1*(0*Rot1(2,1)+2*q11*Rot1(1,2)-q12*-4*Rot1(2,2)-2*q11*Rot1(1,1)) - (Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))*(Rot1(2,2)*0 + Rot1(1,2)*2*q11)/sq1)/(sq1^2);
    dtheta1dq3=(-1/sqrt(1-((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1)^2))*(sq1*(-4*q13*Rot1(2,1)+2*q10*Rot1(1,2)-q13*-4*Rot1(2,2)-q10*-2*Rot1(1,1)) - (Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))*(Rot1(2,2)*(-4*q13) + Rot1(1,2)*-2*q10)/sq1)/(sq1^2);
    dtheta1dq0=(-1/sqrt(1-((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1)^2))*(sq1*(0*Rot1(2,1)+2*q13*Rot1(1,2)-0*Rot1(2,2)-q13*-2*Rot1(1,1)) - (Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))*(Rot1(2,2)*0 + Rot1(1,2)*-2*q13)/sq1)/(sq1^2);
    
    
    H1kp1=[0 0 1/sq1 0 0 0 dd1dq1 dd1dq2 dd1dq3 dd1dq0;...
            0 0 0 0 0 0 dtheta1dq1 dtheta1dq2 dtheta1dq3 dtheta1dq0;];
    S1kp1=XZScannerCovariance; %Experimentally determined
    K1kp1=Pkp1k*H1kp1'/(H1kp1*Pkp1k*H1kp1'+S1kp1);

    
    h1kp1=[sz1/sq1;...
           acos((Rot1(1,2)*Rot1(2,1)-Rot1(1,1)*Rot1(2,2))/sq1);];
        
    x1kp1kp1=xkp1k+k1kp1*(z1kp1-h1kp1);
    P1kp1kp1=(eye(10,10)-K1kp1*H1kp1)*Pkp1k;
else
    x1kp1kp1=xkp1k;
    P1kp1kp1=Pkp1k;
end


% CORRECTIVE STEP 2, YZ LASER SCANNER
if execution(2)==1
        
    rz2=x1kp1kp1(3);
    q21=x1kp1kp1(7);  
    q22=x1kp1kp1(8);
    q23=x1kp1kp1(9);
    q20=x1kp1kp1(10);
    Rot2=[1-2*q22^2-2*q23^2 2*(q21*q22-q20*q23) 2*(q21*q23+q20*q22);...
     2*(q21*q22+q20*q23) 1-2*q21^2-2*q23^2 2*(q22*q23-q20*q21);...
     2*(q21*q23-q20*q22) 2*(q22*q23+q20*q21) 1-2*q21^2-2*q22^2;];
    sz2=(Rot2(3,:)*p2+rz2);
    sq2=sqrt(Rot2(2,1)^+Rot2(1,1)^2);
    
    
    dd2dq1=(sq2*([2*q23 2*q20 -4*q21]*p1)-sz2*(Rot2(2,1)*(2*q22) + Rot2(1,1)*0)/sq2)/(sq2^2);
    dd2dq2=(sq2*([-2*q20 2*q23 -4*q22]*p1)-sz2*(Rot2(2,1)*2*q21 + Rot2(1,1)*-4*q22)/sq2)/(sq2^2);
    dd2dq3=(sq2*([2*q21 2*q22 0]*p1)-sz2*(Rot2(2,1)*(2*q20) + Rot2(1,1)*-4*q23)/sq2)/(sq2^2);
    dd2dq0=(sq2*([-2*q22 2*q21 0]*p1)-sz2*(Rot2(2,1)*2*q23 + Rot2(1,1)*0)/sq2)/(sq2^2);
    
    dtheta2dq1=(1/sqrt(1-((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2)^2))*(sq2*((-4*q21)*Rot2(2,1)+2*q22*Rot2(1,2)-0*Rot2(2,2)-2*q22*Rot2(1,1)) - (Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))*(Rot2(2,2)*(Rot2(2,1)*(2*q22) + Rot2(1,1)*0)/sq2)/(sq2^2);
    dtheta2dq2=(1/sqrt(1-((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2)^2))*(sq2*(0*Rot2(2,1)+2*q21*Rot2(1,2)-q22*-4*Rot2(2,2)-2*q21*Rot2(1,1)) - (Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))*(Rot2(2,1)*2*q21 + Rot2(1,1)*-4*q22)/sq2)/(sq2^2);
    dtheta2dq3=(1/sqrt(1-((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2)^2))*(sq2*(-4*q23*Rot2(2,1)+2*q20*Rot2(1,2)-q23*-4*Rot2(2,2)-q20*-2*Rot2(1,1)) - (Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))*(Rot2(2,1)*(2*q20) + Rot2(1,1)*-4*q23)/sq2)/(sq2^2);
    dtheta2dq0=(1/sqrt(1-((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2)^2))*(sq2*(0*Rot2(2,1)+2*q23*Rot2(1,2)-0*Rot2(2,2)-q23*-2*Rot2(1,1)) - (Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))*(Rot2(2,1)*2*q23 + Rot2(1,1)*0)/sq2)/(sq2^2);
    
    
    H2kp1=[0 0 1/sq2 0 0 0 dd2dq1 dd2dq2 dd2dq3 dd2dq0;...
            0 0 0 0 0 0 dtheta2dq1 dtheta2dq2 dtheta2dq3 dtheta2dq0;];
    S2kp1=YZScannerCovariance; %Experimentally determined
    K2kp1=P1kp1kp1*H2kp1'/(H2kp1*P1kp1kp1*H2kp1'+S2kp1);

    
    h2kp1=[sz2/sq2;...
           acos((Rot2(1,2)*Rot2(2,1)-Rot2(1,1)*Rot2(2,2))/sq2);];
        
    x2kp1kp1=x1kp1kp1+k2kp1*(z2kp1-h2kp1);
    P2kp1kp1=(eye(10,10)-K2kp1*H2kp1)*P1kp1kp1;
else
    x2kp1kp1=x1kp1kp1;
    P2kp1kp1=Pkp1kp1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = Rot*p1+c;
s2 = Rot*p2+c;
s3 = Rot*p3+c;


d1 = s1(3,1)/(sqrt(Rot(2,2)^2 + Rot(1,2)^2));
theta1 = acos(((Rot(1,2)*Rot(2,1))-(Rot(1,1)*Rot(2,2)))/(sqrt(Rot(2,2)^2 + Rot(1,2)^2)));

d2 = s2(3,1)/(sqrt(Rot(2,1)^2 + Rot(1,1)^2));
theta2 = acos(((Rot(1,1)*Rot(2,2))-(Rot(2,1)*Rot(1,2)))/(sqrt(Rot(2,1)^2 + Rot(1,1)^2)));

d3 = (tck - s3(2,1))/(sqrt(Rot(3,3)^2 + Rot(1,3)^2));
theta3 = acos(((Rot(1,1)*Rot(3,3))-(Rot(1,3)*Rot(3,1)))/(sqrt(Rot(3,3)^2 + Rot(1,3)^2)));

%UPDATE 1

%% Initializing the H matrix for the Update step %%%%%%%%%%%%%%%%%%%%%%%%%%%

H1kp1=[ 1 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0;];
S1kp1= 0; %covariance matrix of error/uncertainty in relevant measurement outputs





% if rcond(H1kp1*Pkp1k*H1kp1'+S1kp1)<0.00000000001
%    k
 %   H1kp1*Pkp1k*H1kp1'+S1kp1
%end


K1kp1=Pkp1k*(H1kp1')/(H1kp1*Pkp1k*H1kp1'+S1kp1);
%Pkp1k
counter = counter+1;

%     %%%%%%%%%%%%%%%%%%%%%%Update 1%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Counter1==c
    x1kp1kp1=xkp1k+K1kp1*(z1kp1-H1kp1*xkp1k);
    P1kp1kp1=(eye(10,10)-K1kp1*H1kp1)*Pkp1k;
else
    x1kp1kp1=xkp1k;
    P1kp1kp1=Pkp1k;
end
if ((normQuat>1.000000000001) || (normQuat<0.999999999999))
    if ((normQuat>1.00000001) || (normQuat<0.99999999))
        fprintf('**************************************WARNING, quaternions norm not equal 1!!*************************************');
        k
        x1kp1kp1(7:10)
        break
    end
end







predictedState = x2kp1kp1;
predictedCovariance= P2kp1kp1;
    
    
    
    
    
    
    
    

