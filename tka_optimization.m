%% Final Project - George Durrant


%% Objective Function

% f(x) = w1*ovm1 + w2*ovm2
% x = tp, rf_fron, rt_fron, rf_sag, rt_sag
% x0 is initial design variables values
x0 = [5, 10, -20, 15, -50];

options = optimset('LargeScale', 'off');
[finalx, fval] = fmincon('obfun',x0,[],[],[],[],[],[],'confun',options);
[cin,ceq] = confun(finalx);


%% Find strain for tp=5mm, tp=8mm, tp=11mm

% INITIALIZE
L = 80; % mm, length
W = 40; % mm, width
H = 40; % mm, bone height
LL = 20; % mm, P1 location
LR = 60; % mm, P2 location
v1_tp1L = zeros(1, L+1); 
v1_tp2L = zeros(1, L+1);
v1_tp3L = zeros(1, L+1);
v2_tp1L = zeros(1, L+1);
v2_tp2L = zeros(1, L+1);
v2_tp3L = zeros(1, L+1);
v1_tp1R = zeros(1, L+1); 
v1_tp2R = zeros(1, L+1);
v1_tp3R = zeros(1, L+1);
v2_tp1R = zeros(1, L+1);
v2_tp2R = zeros(1, L+1);
v2_tp3R = zeros(1, L+1);
e1_tp1 = zeros(1, L+1);
e1_tp2 = zeros(1, L+1);
e1_tp3 = zeros(1, L+1);
e2_tp1 = zeros(1, L+1);
e2_tp2 = zeros(1, L+1);
e2_tp3 = zeros(1, L+1);
v3L = zeros(1, L+1);
v3R = zeros(1, L+1);
e3 = zeros(1, L+1);

% LOADING CASES
P1L = 1500;
P1R = 1500;
P2L = 2500;
P2R = 500;

% FIND LAMBDA, K
[lambda_tp1, k_tp1] = find_lambda_k(5);
[lambda_tp2, k_tp2] = find_lambda_k(8);
[lambda_tp3, k_tp3] = find_lambda_k(11);

% tp < 4.1 for no liftoff
[l3, k3] = find_lambda_k(4.1);

% FIND DEFLECTION
for i=1:L+1
    if (i-1)<LL
        x = i-1;
        a = LL;
        b = LR;
    else
        x = L-(i-1);
        a = LR;
        b = LL;
    end
    v1_tp1L(i) = calcDeflection(x, a, b, L, P1L, lambda_tp1, k_tp1);
    v1_tp2L(i) = calcDeflection(x, a, b, L, P1L, lambda_tp2, k_tp2);
    v1_tp3L(i) = calcDeflection(x, a, b, L, P1L, lambda_tp3, k_tp3);
    v2_tp1L(i) = calcDeflection(x, a, b, L, P2L, lambda_tp1, k_tp1);
    v2_tp2L(i) = calcDeflection(x, a, b, L, P2L, lambda_tp2, k_tp2);
    v2_tp3L(i) = calcDeflection(x, a, b, L, P2L, lambda_tp3, k_tp3);
    
    v3L(i) = calcDeflection(x, a, b, L, P2L, l3, k3);
end
for i=1:L+1
    if (i-1)<LR
        x = i-1;
        a = LR;
        b = LL;
    else
        x = L-(i-1);
        a = LL;
        b = LR;
    end
    v1_tp1R(i) = calcDeflection(x, a, b, L, P1R, lambda_tp1, k_tp1);
    v1_tp2R(i) = calcDeflection(x, a, b, L, P1R, lambda_tp2, k_tp2);
    v1_tp3R(i) = calcDeflection(x, a, b, L, P1R, lambda_tp3, k_tp3);
    v2_tp1R(i) = calcDeflection(x, a, b, L, P2R, lambda_tp1, k_tp1);
    v2_tp2R(i) = calcDeflection(x, a, b, L, P2R, lambda_tp2, k_tp2);
    v2_tp3R(i) = calcDeflection(x, a, b, L, P2R, lambda_tp3, k_tp3);
    
    v3R(i) = calcDeflection(x, a, b, L, P2R, l3, k3);
end

% FIND STRAIN
for i=1:L+1
    e1_tp1(i) = -(v1_tp1L(i)+v1_tp1R(i))/H*100;
    e1_tp2(i) = -(v1_tp2L(i)+v1_tp2R(i))/H*100;
    e1_tp3(i) = -(v1_tp3L(i)+v1_tp3R(i))/H*100;
    e2_tp1(i) = -(v2_tp1L(i)+v2_tp1R(i))/H*100;
    e2_tp2(i) = -(v2_tp2L(i)+v2_tp2R(i))/H*100;
    e2_tp3(i) = -(v2_tp3L(i)+v2_tp3R(i))/H*100;
    
    e3(i) = -(v3L(i)+v3R(i))/H*100;
end


%% Plot strain

lin = linspace(0, L, L+1);
figure
plot(lin, e1_tp1); hold on; plot(lin, e1_tp2); plot(lin, e1_tp3);
title('Trabecular Bone Strain, Habitual Loading')
xlabel('Length (mm)'); ylabel('Strain (%)');
legend({'tp=5mm', 'tp=8mm', 'tp=11mm'}, 'Location', 'northwest')

figure
plot(lin, e2_tp1); hold on; plot(lin, e2_tp2); plot(lin, e2_tp3);
title('Trabecular Bone Strain, Extreme Loading')
xlabel('Length (mm)'); ylabel('Strain (%)');
legend({'tp=5mm', 'tp=8mm', 'tp=11mm'}, 'Location', 'northwest')


%% Von Mises stress vs Rt_frontal

% result 4) plot max von mises stress
% with x axis as 14-21 mm for Rtfrontal 
% hold everything else constant, use 3 values of tp
% tp = 5, 8, 11

tp_vals = [5, 8, 11];
% Rt_front_vals = linspace(14, 21, 70);
Rt_front_vals = linspace(-15, -21, 60);
Rf_front = finalx(2);
Rt_sag = finalx(5);
Rf_sag = finalx(4);

ovm1_tp1 = zeros(1, 60);
ovm1_tp2 = zeros(1, 60);
ovm1_tp3 = zeros(1, 60);
ovm2_tp1 = zeros(1, 60);
ovm2_tp2 = zeros(1, 60);
ovm2_tp3 = zeros(1, 60);

for ii=1:3
    
    for jj=1:60
        
        % initialize
        tp = tp_vals(ii);
        R1 = Rf_front;
        R1p = Rf_sag;
%         R2 = 0-Rt_front_vals(jj);
        R2 = Rt_front_vals(jj);
        R2p = Rt_sag;
        psi = 0;
        E1 = 210e3;
        E2 = 1e3;
        v = 0.3;

        % loading cases
        P1L = 1500;
        P1R = 1500;
        P2L = 2500;
        P2R = 500;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND HERTZ CONTACT STRESSES AT POINT A %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [S1, a1] = modHertz(P1L, R1, R1p, R2, R2p, psi, tp, E1, E2, v);
        [S2L, a2L] = modHertz(P2L, R1, R1p, R2, R2p, psi, tp, E1, E2, v);
        [S2R, a2R] = modHertz(P2R, R1, R1p, R2, R2p, psi, tp, E1, E2, v);

        % loading case 1
        or_A1 = S1(1);
        ot_A1 = S1(2);
        oz_A1 = S1(3);
        Pmax1 = -oz_A1;

        % loading case 2, left
        or_A2L = S2L(1);
        ot_A2L = S2L(2);
        oz_A2L = S2L(3);
        Pmax2L = -oz_A2L;

        % loading case 2, right
        or_A2R = S2R(1);
        ot_A2R = S2R(2);
        oz_A2R = S2R(3);
        Pmax2R = -oz_A2R;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND HERTZ CONTACT STRESSES AT POINT B %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % loading case 1
        or_B1 = 0.31*or_A1;
        ot_B1 = 0.31*ot_A1;
        oz_B1 = 2*(0.31*-Pmax1)+or_B1;

        % loading case 2, left
        or_B2L = 0.31*or_A2L;
        ot_B2L = 0.31*ot_A2L;
        oz_B2L = 2*(0.31*-Pmax2L)+or_B2L;

        % loading case 2, right
        or_B2R = 0.31*or_A2R;
        ot_B2R = 0.31*ot_A2R;
        oz_B2R = 2*(0.31*-Pmax2R)+or_B2R;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND HERTZ CONTACT STRESSES AT POINT C %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % loading case 1
        or_C1 = (1/3)*(1-2*v)*Pmax1;
        ot_C1 = -(1/3)*(1-2*v)*Pmax1;
        oz_C1 = 0;

        % loading case 2, left
        or_C2L = (1/3)*(1-2*v)*Pmax2L;
        ot_C2L = -(1/3)*(1-2*v)*Pmax2L;
        oz_C2L = 0;

        % loading case 2, right
        or_C2R = (1/3)*(1-2*v)*Pmax2R;
        ot_C2R = -(1/3)*(1-2*v)*Pmax2R;
        oz_C2R = 0;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND BENDING STRESS USING BOEF %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % initialize
        L = 80; % mm, length
        W = 40; % mm, width
        H = 40; % mm, bone height
        LL = 20; % mm, P1 location
        LR = 60; % mm, P2 location
        Epe = 1e3; % MPa, plastic
        Eco = 210e3; % MPa, metal
        Etra = 500; % MPa, bone
        M1L = zeros(1, L+1); 
        M1R = zeros(1, L+1);
        M2L = zeros(1, L+1);
        M2R = zeros(1, L+1);
        M1tot = zeros(1, L+1);
        M2tot = zeros(1, L+1);
        % find lambda, k, sigEI, NAt
        [lambda, k, sigEI, NAt] = find_lambda_k(tp);
        % calculate moment 
        for i=1:L+1
            % when x<a
            if (i-1)<LL
                x = i-1;
                a = LL;
                b = LR;
            % when a<x<b
            else
                x = L-(i-1);
                a = LR;
                b = LL;
            end
            M1L(i) = calcMoment(x, a, b, L, P1L, lambda);
            M2L(i) = calcMoment(x, a, b, L, P2L, lambda);
        end
        for i=1:L+1
            if (i-1)<LR
                x = i-1;
                a = LR;
                b = LL;
            else
                x = L-(i-1);
                a = LL;
                b = LR;
            end
            M1R(i) = calcMoment(x, a, b, L, P1R, lambda);
            M2R(i) = calcMoment(x, a, b, L, P2R, lambda);
        end
        for i=1:L+1
            M1tot(i) = M1L(i)+M1R(i);
            M2tot(i) = M2L(i)+M2R(i);
        end
        % calculate bending stress
        ob1 = zeros(1, L+1);
        ob2 = zeros(1, L+1);
        for i=1:L+1
            ob1(i) = Epe*M1tot(i)*NAt/sigEI;
            ob2(i) = Epe*M2tot(i)*NAt/sigEI;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND BENDING STRESS AT POINT A %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % loading case 1
        M_A1 = M1tot(LL+1);
        ob_A1 = Epe*M_A1*NAt/sigEI;

        % loading case 2, left
        M_A2L = M2tot(LL+1);
        ob_A2L = Epe*M_A2L*NAt/sigEI;

        % loading case 2, right
        M_A2R = M2tot(LL+1);
        ob_A2R = Epe*M_A2R*NAt/sigEI;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND BENDING STRESS AT POINT B %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % loading case 1
        M_B1 = M1tot(LL+1);
        NAtB1 = NAt+0.51*a1;
        ob_B1 = Epe*M_B1*NAtB1/sigEI;

        % loading case 2, left
        M_B2L = M2tot(LL+1);
        NAtB2L = NAt+0.51*a2L;
        ob_B2L = Epe*M_B2L*NAtB2L/sigEI;

        % loading case 2, right
        M_B2R = M2tot(LR+1);
        NAtB2R = NAt+0.51*a2R;
        ob_B2R = Epe*M_B2R*NAtB2R/sigEI;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND BENDING STRESS AT POINT C %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % loading case 1
        xp1_mina = LL-a1;
        Mp1L_mina = calcMoment(xp1_mina, LL, LR, L, P1L, lambda);
        Mp1R_mina = calcMoment(xp1_mina, LR, LL, L, P1R, lambda);
        Mp1_mina = Mp1L_mina+Mp1R_mina;

        xp1_plusa = LL+a1;
        Mp1L_plusa = calcMoment(L-xp1_plusa, LR, LL, L, P1L, lambda);
        Mp1R_plusa = calcMoment(xp1_plusa, LR, LL, L, P1R, lambda);
        Mp1_plusa = Mp1L_plusa+Mp1R_plusa;

        M_C1 = max(Mp1_mina, Mp1_plusa);
        ob_C1 = Epe*M_C1*NAt/sigEI;

        % loading case 2, left
        xp2_minaL = LL-a2L;
        Mp2L_minaL = calcMoment(xp2_minaL, LL, LR, L, P2L, lambda);
        Mp2R_minaL = calcMoment(xp2_minaL, LR, LL, L, P2R, lambda);
        Mp2_minaL = Mp2L_minaL+Mp2R_minaL;

        xp2_plusaL = LL+a2L;
        Mp2L_plusaL = calcMoment(L-xp2_plusaL, LR, LL, L, P2L, lambda);
        Mp2R_plusaL = calcMoment(xp2_plusaL, LR, LL, L, P2R, lambda);
        Mp2_plusaL = Mp2L_plusaL+Mp2R_plusaL;

        M_C2L = max(Mp2_minaL, Mp2_plusaL);
        ob_C2L = Epe*M_C2L*NAt/sigEI;

        % loading case 2, right
        xp2_minaR = LR-a2R;
        Mp2L_minaR = calcMoment(L-xp2_minaR, LR, LL, L, P2L, lambda);
        Mp2R_minaR = calcMoment(xp2_minaR, LR, LL, L, P2R, lambda);
        Mp2_minaR = Mp2L_minaR+Mp2R_minaR;

        xp2_plusaR = LR+a2R;
        Mp2L_plusaR = calcMoment(L-xp2_plusaR, LR, LL, L, P2L, lambda);
        Mp2R_plusaR = calcMoment(L-xp2_plusaR, LL, LR, L, P2R, lambda);
        Mp2_plusaR = Mp2L_plusaR+Mp2R_plusaR;

        M_C2R = max(Mp2_minaR, Mp2_plusaR);
        ob_C2R = Epe*M_C2R*NAt/sigEI;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ADD BENDING STRESSES TO RADIAL HERTZ STRESSES %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % point A
        orb_A1 = or_A1+ob_A1;
        orb_A2L = or_A2L+ob_A2L;
        orb_A2R = or_A2R+ob_A2R;

        % point B 
        orb_B1 = or_B1+ob_B1;
        orb_B2L = or_B2L+ob_B2L;
        orb_B2R = or_B2R+ob_B2R;

        % point C
        orb_C1 = or_C1+ob_C1;
        orb_C2L = or_C2L+ob_C2L;
        orb_C2R = or_C2R+ob_C2R;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND VON MISES STRESSES %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        % loading case 1
        ovm1A = sqrt((1/2)*( (orb_A1-ot_A1)^2 + (ot_A1-oz_A1)^2 + (oz_A1-orb_A1)^2));
        ovm1B = sqrt((1/2)*( (orb_B1-ot_B1)^2 + (ot_B1-oz_B1)^2 + (oz_B1-orb_B1)^2));
        ovm1C = sqrt((1/2)*( (orb_C1-ot_C1)^2 + (ot_C1-oz_C1)^2 + (oz_C1-orb_C1)^2));

        % loading case 2, left
        ovm2AL = sqrt((1/2)*( (orb_A2L-ot_A2L)^2 + (ot_A2L-oz_A2L)^2 + (oz_A2L-orb_A2L)^2));
        ovm2BL = sqrt((1/2)*( (orb_B2L-ot_B2L)^2 + (ot_B2L-oz_B2L)^2 + (oz_B2L-orb_B2L)^2));
        ovm2CL = sqrt((1/2)*( (orb_C2L-ot_C2L)^2 + (ot_C2L-oz_C2L)^2 + (oz_C2L-orb_C2L)^2));

        % loading case 2, right
        ovm2AR = sqrt((1/2)*( (orb_A2R-ot_A2R)^2 + (ot_A2R-oz_A2R)^2 + (oz_A2R-orb_A2R)^2));
        ovm2BR = sqrt((1/2)*( (orb_B2R-ot_B2R)^2 + (ot_B2R-oz_B2R)^2 + (oz_B2R-orb_B2R)^2));
        ovm2CR = sqrt((1/2)*( (orb_C2R-ot_C2R)^2 + (ot_C2R-oz_C2R)^2 + (oz_C2R-orb_C2R)^2));

        % max von mises stresses
        
        if ii==1  
            ovm1_tp1(jj) = max([ovm1A, ovm1B, ovm1C]);
            ovm2_tp1(jj) = max([ovm2AL, ovm2BL, ovm2CL, ovm2AR, ovm2BR, ovm2CR]);
        elseif ii==2
            ovm1_tp2(jj) = max([ovm1A, ovm1B, ovm1C]);
            ovm2_tp2(jj) = max([ovm2AL, ovm2BL, ovm2CL, ovm2AR, ovm2BR, ovm2CR]);
        else
            ovm1_tp3(jj) = max([ovm1A, ovm1B, ovm1C]);
            ovm2_tp3(jj) = max([ovm2AL, ovm2BL, ovm2CL, ovm2AR, ovm2BR, ovm2CR]);
        end
    end
end

%% Plot Von Mises stress vs Rt_frontal
Rtlin = linspace(15, 21, 60);
figure
plot(Rtlin, ovm1_tp1); hold on;
plot(Rtlin, ovm1_tp2); plot(Rtlin, ovm1_tp3);
title('Von Mises stress vs R_tfrontal, Habitual Loading')
xlabel('R_tfrontal (mm)'); ylabel('Von Mises stress (MPa');
legend({'tp=5mm', 'tp=8mm', 'tp=11mm'}, 'Location', 'northwest')

figure 
plot(Rtlin, ovm2_tp1); hold on;
plot(Rtlin, ovm2_tp2); plot(Rtlin, ovm2_tp3);
title('Von Mises stress vs R_tfrontal, Extreme Loading')
xlabel('R_tfrontal (mm)'); ylabel('Von Mises stress (MPa');
legend({'tp=5mm', 'tp=8mm', 'tp=11mm'}, 'Location', 'northwest')


%% Metal Stresses

NAtm = NAt+tp;

% loading case 1
M_1m = M1tot(LL+1);
ob_1m = Eco*M_1m*NAtm/sigEI;

% loading case 2, left
M_2Lm = M2tot(LL+1);
ob_2Lm = Eco*M_2Lm*NAtm/sigEI;

% loading case 2, right
M_2Rm = M2tot(LR+1);
ob_2Rm = Eco*M_2Rm*NAtm/sigEI;

%% Bone Stresses

NAtb = NAt-NAt+(12-tp)/2;

% loading case 1
M_1b = M1tot(LL+1);
ob_1b = Etra*M_1b*NAtb/sigEI;

% loading case 2, left
M_2Lb = M2tot(LL+1);
ob_2Lb = Etra*M_2Lb*NAtb/sigEI;

% loading case 2, right
M_2Rb = M2tot(LR+1);
ob_2Rb = Etra*M_2Rb*NAtb/sigEI;

