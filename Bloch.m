function [Mxf,Myf,Mzf] = Bloch(M0, T1, T2, tdt,tda, TSL, sla, dw0, it, ip, nca, wNC,phi, T1r, T2r,alphaTd)

dt = tdt(end)/length(tdt);
rfrot = tda.*dt; % Rotation in radians at each time point

M(:,1) = M0'; %magnetization at the beginning

Mz0 = M(3,1);

%% Tip down pulse
A = @(t) [exp(-t/T2) 0 0;0 exp(-t/T2) 0;0 0 exp(-t/T1)];
B = @(t) [0; 0; Mz0*(1-exp(-t/T1))];

for td = 1:length(tdt)
    M(:,td+1) = A(dt/2)*M(:,td)+B(dt/2);
    M(:,td+1) = Rot('x',rfrot(1),'lh')*M(:,td+1);	% RF Rotation.
    M(:,td+1) = A(dt/2)*M(:,td+1)+B(dt/2);
end

%% SL time

A = @(t) [exp(-t/T2r) 0 0;0 exp(-t/T1r) 0;0 0 exp(-t/T2r)];
B = @(t) [0; 0; 0]; % Mz0*cos(theta)*(1-exp(-t/T1r)) 

% First SL pulse
samples = round(TSL/dt);
for sl1 = 1:samples/2
    nca_t = nca*sin(wNC*sl1*dt + phi); %BNC amplitude
    theta = atan(sla(sl1)/(dw0 + nca_t));
    weff = sign(dw0 + nca_t + 1E-20).*sqrt(sla(sl1).^2 + (dw0 + nca_t)^2);
    
    M(:,td+sl1+1) = A(dt/2)*M(:,td+sl1)+B(dt/2);
    M(:,td+sl1+1) = Rot('x',theta,'lh')*Rot('z',weff*dt,'lh')*Rot('x',-theta,'lh')* M(:,td+sl1+1);	% SL Rotation.
    M(:,td+sl1+1) = A(dt/2)*M(:,td+sl1+1)+B(dt/2);
end

A = @(t) [exp(-t/T2) 0 0;0 exp(-t/T2) 0;0 0 exp(-t/T1)];
B = @(t) [0; 0; Mz0*(1-exp(-t/T1))];
% Inversion pulse
if it ~= 0
    ips = round(it/dt);
    ipa = -2*alphaTd/it;
    for is = 1:ips
        M(:,td+sl1+is+1) = A(dt/2)*M(:,td+sl1+is)+B(dt/2);
        M(:,td+sl1+is+1) = Rot('y',ipa*dt,'lh')* M(:,td+sl1+is+1);	% 180 RF Rotation.
        M(:,td+sl1+is+1) = A(dt/2)*M(:,td+sl1+is+1)+B(dt/2);
    end
else
    is = 0;
end

A = @(t) [exp(-t/T2r) 0 0;0 exp(-t/T1r) 0;0 0 exp(-t/T2r)];
B = @(t) [0; 0; 0];%Mz0*cos(theta)*(1-exp(-t/T1r)) 
% Second SL pulse
for sl2 = 1:samples/2
    nca_t = nca*sin(wNC*(sl1+sl2+is)*dt + phi);
    theta = atan(sla(sl1+sl2)/(dw0 + nca_t));
    weff = sign(dw0 + nca_t + 1E-20).*sqrt(sla(sl1+sl2).^2 + (dw0 + nca_t)^2);

    M(:,td+sl1+is+sl2+1) = A(dt/2)*M(:,td+sl1+is+sl2)+B(dt/2);
    M(:,td+sl1+is+sl2+1) = Rot('x',theta,'lh')*Rot('z',weff*dt,'lh')*Rot('x',-theta,'lh')* M(:,td+sl1+is+sl2+1);	% RF Rotation, notice theta has the oposite sign
    M(:,td+sl1+is+sl2+1) = A(dt/2)*M(:,td+sl1+is+sl2+1)+B(dt/2);
end

A = @(t) [exp(-t/T2) 0 0;0 exp(-t/T2) 0;0 0 exp(-t/T1)];
B = @(t) [0; 0; Mz0*(1-exp(-t/T1))];

%% Tip up pulse

A = @(t) [exp(-t/T2) 0 0;0 exp(-t/T2) 0;0 0 exp(-t/T1)];
B = @(t) [0; 0; Mz0*(1-exp(-t/T1))];

for tu = 1:length(tdt)
    M(:,td+sl1+is+sl2+tu+1) = A(dt/2)*M(:,td+sl1+is+sl2+tu)+B(dt/2);
    M(:,td+sl1+is+sl2+tu+1) = Rot('x',ip*rfrot(1),'lh')* M(:,td+sl1+is+sl2+tu+1);	% RF Rotation.
    M(:,td+sl1+is+sl2+tu+1) = A(dt/2)*M(:,td+sl1+is+sl2+tu+1)+B(dt/2);
end


Mxf = M(1,:);
Myf = M(2,:);
Mzf = M(3,:);
