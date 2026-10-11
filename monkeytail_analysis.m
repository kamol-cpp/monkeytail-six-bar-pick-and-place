function results = monkeytail_analysis(overrides)
%MONKEYTAIL_ANALYSIS Kinematics and inverse dynamics of the MonkeyTail six-bar.
%   SI units
%   Coordinates are relative to O4, with x right, y up, and CCW positive.
%   Six links: ground; crank O2-A; coupler A-B; rocker O4-B-C;
%   sliding bar C-E; and a massless swivel block pinned to ground at O5.
%   Forces
%     F_A = force ON CRANK at A from coupler;
%     F_B = force ON COUPLER at B from rocker;
%     F_C = force ON ROCKER at C from sliding bar;
%     R_O2, R_O4 = ground forces on crank and rocker;
%     R_O5 = signed scalar force on sliding bar along nT = perp(eT);
%     tau = required crank-shaft torque, positive CCW.

    cfg = defaultSettings();
    if nargin == 1 && (ischar(overrides) || isstring(overrides))
        assert(strcmpi(overrides,'defaults'), ...
            'MonkeyTail:Input','The only text option is ''defaults''.');
        results = cfg;
        return
    end
    if nargin == 1
        assert(isstruct(overrides) && isscalar(overrides), ...
            'MonkeyTail:Input','Overrides must be a scalar struct.');
        fields = fieldnames(overrides);
        for k = 1:numel(fields)
            assert(isfield(cfg,fields{k}),'MonkeyTail:Input', ...
                'Unknown setting: %s.',fields{k});
            cfg.(fields{k}) = overrides.(fields{k});
        end
    end

    p = modelParameters(cfg);
    theta2 = (0:cfg.numSamples-1) * (2*pi/cfg.numSamples);
    motion = kinematics(theta2,p);
    dynamics = inverseDynamics(motion,p);
    validation = validateModel(motion,dynamics,p);
    extrema = turningPoints(motion,p);
    summary = forceSummary(theta2,dynamics);

    results = struct('settings',cfg,'model',p,'motion',motion, ...
        'dynamics',dynamics,'validation',validation, ...
        'turningPoints',extrema,'forceSummary',summary);
    if cfg.verbose
        printSummary(results);
    end
    results.figures = {};
    if cfg.makePlots
        results.figures{1} = plotResults(results);
    end
    if nargout == 0
        clear results
    end
end

function cfg = defaultSettings()
    cfg.rpm = 20;
    cfg.numSamples = 1440;
    cfg.groundLength_in = 8.76;
    cfg.groundAngle_rad = pi/2.8851;
    cfg.crankRadius_in = 2.63;
    cfg.couplerLength_in = 8.34;
    cfg.rockerJointRadius_in = 3.7;
    cfg.rockerTipRadius_in = 6;
    cfg.sliderSupport_in = [-6;6];
    cfg.sliderLength_in = 7;
    cfg.linkWidths_in = [1,1,1,1.5];
    cfg.thickness_mm = 6;
    cfg.density_lbm_in3 = 0.0433;
    cfg.cadVolumes_in3 = [0.8,2.15,1.59];
    cfg.gripper_g = 40;
    cfg.payload_g = 100;
    cfg.gripperOffset_in = 0.5;
    cfg.gripperSize_in = [2.5,2.5,1];
    cfg.assemblyBranch = 1;
    cfg.gravity_m_s2 = 9.80665;
    cfg.motorRatedTorque_kgfcm = 11;
    cfg.makePlots = true;
    cfg.figureVisible = 'on';
    cfg.verbose = true;
end

function p = modelParameters(c)
    checkNumbers(c.rpm,1,0,true);
    checkNumbers(c.numSamples,1,360,true);
    assert(mod(c.numSamples,1)==0,'numSamples must be an integer.');
    positive = [c.groundLength_in,c.crankRadius_in,c.couplerLength_in, ...
        c.rockerJointRadius_in,c.rockerTipRadius_in,c.sliderLength_in, ...
        c.thickness_mm,c.density_lbm_in3,c.gravity_m_s2];
    checkNumbers(positive,9,0,false);
    checkNumbers(c.linkWidths_in,4,0,false);
    checkNumbers(c.cadVolumes_in3,3,0,false);
    checkNumbers(c.gripperSize_in,3,0,false);
    checkNumbers(c.sliderSupport_in,2,-Inf,true);
    checkNumbers(c.groundAngle_rad,1,-Inf,true);
    checkNumbers(c.gripperOffset_in,1,0,true);
    checkNumbers(c.gripper_g,1,0,true);
    checkNumbers(c.payload_g,1,0,true);
    assert(ismember(c.assemblyBranch,[-1,1]) && isscalar(c.assemblyBranch), ...
        'MonkeyTail:Branch','assemblyBranch must be +1 or -1.');
    assert(strcmp(c.figureVisible,'on') || strcmp(c.figureVisible,'off'), ...
        'MonkeyTail:Input','figureVisible must be on or off.');
    checkNumbers(c.motorRatedTorque_kgfcm,1,0,false);
    assert(islogical(c.makePlots) && isscalar(c.makePlots), ...
        'makePlots must be true or false.');
    assert(islogical(c.verbose) && isscalar(c.verbose), ...
        'verbose must be true or false.');
    assert(c.rockerJointRadius_in <= c.rockerTipRadius_in, ...
        'MonkeyTail:Geometry','B must lie between O4 and C.');

    inch = 0.0254;
    lbm = 0.45359237;
    p.inch = inch;
    p.lbm = lbm;
    p.g = c.gravity_m_s2;
    p.omega2 = c.rpm*2*pi/60;
    p.alpha2 = 0;
    p.l2 = c.crankRadius_in*inch;
    p.l3 = c.couplerLength_in*inch;
    p.rB = c.rockerJointRadius_in*inch;
    p.rC = c.rockerTipRadius_in*inch;
    p.barLength = c.sliderLength_in*inch;
    p.gripOffset = c.gripperOffset_in*inch;
    p.O4 = [0;0];
    p.O2 = c.groundLength_in*inch* ...
        [cos(c.groundAngle_rad);sin(c.groundAngle_rad)];
    p.O5 = c.sliderSupport_in(:)*inch;
    p.branch = c.assemblyBranch;
    masses = c.density_lbm_in3*c.cadVolumes_in3(:)*lbm;
    p.m2 = masses(1); p.m3 = masses(2); p.m4 = masses(3);
    p.mBar = c.density_lbm_in3*c.sliderLength_in* ...
        c.linkWidths_in(4)*(c.thickness_mm/25.4)*lbm;
    p.mGripPayload = (c.gripper_g+c.payload_g)/1000;
    p.mT = p.mBar+p.mGripPayload;
    widths = c.linkWidths_in*inch;
    p.I2 = p.m2*(p.l2^2+widths(1)^2)/12;
    p.I3 = p.m3*(p.l3^2+widths(2)^2)/12;
    p.I4 = p.m4*(p.rC^2+widths(3)^2)/12;
    p.sT = (p.mBar*p.barLength/2- ...
        p.mGripPayload*p.gripOffset)/p.mT;
    gripSize = c.gripperSize_in*inch;
    p.IT = p.mBar*(p.barLength^2+widths(4)^2)/12 + ...
        p.mBar*(p.barLength/2-p.sT)^2 + ...
        p.mGripPayload*(gripSize(1)^2+gripSize(2)^2)/12 + ...
        p.mGripPayload*(-p.gripOffset-p.sT)^2;
    p.lengthScale = max([p.l2,p.l3,p.rC,p.barLength,norm(p.O2)]);
    p.forceScale = (p.m2+p.m3+p.m4+p.mT)*p.g;
    p.motorRatedTorque_Nm = c.motorRatedTorque_kgfcm*0.0980665;
end

function s = kinematics(theta2,p)
    n = numel(theta2);
    s.theta2 = theta2;
    e2 = [cos(theta2);sin(theta2)];
    s.A = p.O2+p.l2*e2;
    d = vectorLength(s.A-p.O4);
    assert(all(d > 1e-12*p.lengthScale),'MonkeyTail:Geometry', ...
        'Coincident circle centers; assembly is indeterminate.');
    u = (s.A-p.O4)./d;
    along = (p.rB^2-p.l3^2+d.^2)./(2*d);
    h2 = p.rB^2-along.^2;
    assert(all(h2 > 1e-12*p.lengthScale^2), ...
        'MonkeyTail:Geometry', ...
        'Four-bar cannot close or reaches a true singular configuration.');
    s.B = p.O4+along.*u+p.branch*sqrt(h2).*perp(u);
    e4 = (s.B-p.O4)/p.rB;
    e3 = (s.A-s.B)/p.l3; % theta3 is B -> A, NOT A -> B.
    s.C = p.O4+p.rC*e4;
    s.theta3 = unwrap(atan2(e3(2,:),e3(1,:)));
    s.theta4 = unwrap(atan2(e4(2,:),e4(1,:)));
    s.vA = p.l2*p.omega2*perp(e2);
    s.aA = p.l2*(p.alpha2*perp(e2)-p.omega2^2*e2);
    s.omega3 = zeros(1,n); s.omega4 = zeros(1,n);
    s.alpha3 = zeros(1,n); s.alpha4 = zeros(1,n);
    s.kinematicRcond = zeros(1,n);
    for k = 1:n
        J = [p.l3*perp(e3(:,k)),p.rB*perp(e4(:,k))];
        s.kinematicRcond(k) = rcond(J/p.lengthScale);
        assert(s.kinematicRcond(k)>1e-10,'MonkeyTail:Singularity', ...
            'Four-bar Jacobian is nearly singular at %.6g deg.', ...
            rad2deg(theta2(k)));
        omega = J\s.vA(:,k);
        alpha = J\(s.aA(:,k)+p.l3*omega(1)^2*e3(:,k)+ ...
            p.rB*omega(2)^2*e4(:,k));
        s.omega3(k) = omega(1); s.omega4(k) = omega(2);
        s.alpha3(k) = alpha(1); s.alpha4(k) = alpha(2);
    end
    s.vB = p.rB*s.omega4.*perp(e4);
    s.aB = p.rB*(s.alpha4.*perp(e4)-s.omega4.^2.*e4);
    s.vC = p.rC*s.omega4.*perp(e4);
    s.aC = p.rC*(s.alpha4.*perp(e4)-s.omega4.^2.*e4);
    s.sliderDistance = vectorLength(p.O5-s.C);
    assert(all(s.sliderDistance>1e-10*p.lengthScale), ...
        'MonkeyTail:Geometry','C coincides with O5; bar angle is undefined.');
    s.eT = (p.O5-s.C)./s.sliderDistance;
    s.nT = perp(s.eT);
    s.thetaT = unwrap(atan2(s.eT(2,:),s.eT(1,:)));
    s.sliderSpeed = -sum(s.vC.*s.eT,1);
    s.omegaT = -sum(s.vC.*s.nT,1)./s.sliderDistance;
    s.sliderAcceleration = -sum(s.aC.*s.eT,1)+ ...
        s.sliderDistance.*s.omegaT.^2;
    s.alphaT = (-sum(s.aC.*s.nT,1)- ...
        2*s.sliderSpeed.*s.omegaT)./s.sliderDistance;
    s.E = s.C+p.barLength*s.eT;
    s.P = s.C-p.gripOffset*s.eT;
    angularT = s.alphaT.*s.nT-s.omegaT.^2.*s.eT;
    s.vP = s.vC-p.gripOffset*s.omegaT.*s.nT;
    s.aP = s.aC-p.gripOffset*angularT;
    s.G2 = p.O2+p.l2/2*e2;
    s.G3 = (s.A+s.B)/2;
    s.G4 = p.O4+p.rC/2*e4;
    s.GT = s.C+p.sT*s.eT;
    s.vG2 = p.l2/2*p.omega2*perp(e2);
    s.aG2 = p.l2/2*(p.alpha2*perp(e2)-p.omega2^2*e2);
    s.vG3 = (s.vA+s.vB)/2; s.aG3 = (s.aA+s.aB)/2;
    s.vG4 = s.vC/2; s.aG4 = s.aC/2;
    s.vGT = s.vC+p.sT*s.omegaT.*s.nT;
    s.aGT = s.aC+p.sT*angularT;
    s.transmissionAngle_deg = rad2deg(acos(min(1, ...
        abs(sum(e3.*e4,1)))));
    s.idealRockerTorqueRatio = nan(1,n);
    if p.omega2 > 0
        moving = abs(s.omega4)>1e-8*p.omega2;
        s.idealRockerTorqueRatio(moving) = ...
            p.omega2./s.omega4(moving);
    end
    s.sliderEngagementMargin = p.barLength-s.sliderDistance;
end

function f = inverseDynamics(s,p)
    n = numel(s.theta2);
    values = zeros(12,n);
    f.scaledResidual_N = zeros(1,n);
    f.dynamicsRcond = zeros(1,n);
    momentRows = [3,6,9,12];
    for k = 1:n
        A = zeros(12); b = zeros(12,1);
        % Crank
        A(1:2,1:2) = eye(2); A(1:2,7:8) = eye(2);
        A(3,1:2) = momentRow(s.A(:,k)-s.G2(:,k));
        A(3,7:8) = momentRow(p.O2-s.G2(:,k)); A(3,12) = 1;
        b(1:2) = p.m2*(s.aG2(:,k)+[0;p.g]);
        b(3) = p.I2*p.alpha2;
        % Coupler
        A(4:5,1:2) = -eye(2); A(4:5,3:4) = eye(2);
        A(6,1:2) = -momentRow(s.A(:,k)-s.G3(:,k));
        A(6,3:4) = momentRow(s.B(:,k)-s.G3(:,k));
        b(4:5) = p.m3*(s.aG3(:,k)+[0;p.g]);
        b(6) = p.I3*s.alpha3(k);
        % Rocker
        A(7:8,3:4) = -eye(2); A(7:8,5:6) = eye(2);
        A(7:8,9:10) = eye(2);
        A(9,3:4) = -momentRow(s.B(:,k)-s.G4(:,k));
        A(9,5:6) = momentRow(s.C(:,k)-s.G4(:,k));
        A(9,9:10) = momentRow(p.O4-s.G4(:,k));
        b(7:8) = p.m4*(s.aG4(:,k)+[0;p.g]);
        b(9) = p.I4*s.alpha4(k);
        % Bar+gripper+payload
        A(10:11,5:6) = -eye(2); A(10:11,11) = s.nT(:,k);
        A(12,5:6) = -momentRow(s.C(:,k)-s.GT(:,k));
        A(12,11) = cross2(p.O5-s.GT(:,k),s.nT(:,k));
        b(10:11) = p.mT*(s.aGT(:,k)+[0;p.g]);
        b(12) = p.IT*s.alphaT(k);
        As = A;
        As(momentRows,:) = A(momentRows,:)/p.lengthScale;
        As(:,12) = As(:,12)*p.lengthScale;
        bs = b;
        bs(momentRows) = b(momentRows)/p.lengthScale;
        f.dynamicsRcond(k) = rcond(As);
        assert(f.dynamicsRcond(k)>1e-10,'MonkeyTail:Singularity', ...
            'Dynamics matrix is nearly singular at %.6g deg.', ...
            rad2deg(s.theta2(k)));
        z = As\bs;
        values(:,k) = z;
        values(12,k) = z(12)*p.lengthScale;
        f.scaledResidual_N(k) = norm(As*z-bs,inf);
    end
    assert(all(isfinite(values(:))),'MonkeyTail:Nonfinite', ...
        'Nonfinite force or torque result.');
    f.F_A = values(1:2,:); f.F_B = values(3:4,:);
    f.F_C = values(5:6,:); f.R_O2 = values(7:8,:);
    f.R_O4 = values(9:10,:); f.R_O5 = values(11,:);
    f.R_O5_vector = s.nT.*f.R_O5;
    f.tau_Nm = values(12,:);
    f.power_W = f.tau_Nm*p.omega2;
    f.tauRMS_Nm = sqrt(mean(f.tau_Nm.^2));
    f.tauPeak_Nm = max(abs(f.tau_Nm));
    f.idealRatedTorqueRatio = p.motorRatedTorque_Nm/ ...
        max(f.tauPeak_Nm,realmin);
end

function v = validateModel(s,f,p)
    closure = s.A-s.B-p.l3*[cos(s.theta3);sin(s.theta3)];
    sliderClosure = s.C+s.sliderDistance.*s.eT-p.O5;
    v.loopClosure_m = max([vectorLength(closure),vectorLength(sliderClosure)]);
    v.reactionTangency_N = max(abs(sum(f.R_O5_vector.*s.eT,1)));
    v.newtonEulerResidual_N = max(f.scaledResidual_N);
    v.minimumKinematicRcond = min(s.kinematicRcond);
    v.minimumDynamicsRcond = min(f.dynamicsRcond);
    v.minimumSliderMargin_m = min(s.sliderEngagementMargin);
    assert(v.minimumSliderMargin_m>0,'MonkeyTail:Travel', ...
        'The finite bar leaves the swivel: minimum margin %.6g m.', ...
        v.minimumSliderMargin_m);
    v.maximumPoseStep_m = max(vectorLength(diff(s.B,1,2)));
    assert(v.maximumPoseStep_m<0.1*p.lengthScale, ...
        'MonkeyTail:Branch','Discontinuous or inadequately sampled assembly.');
    masses = [p.m2,p.m3,p.m4,p.mT];
    inertias = [p.I2,p.I3,p.I4,p.IT];
    pos = {s.G2,s.G3,s.G4,s.GT};
    vel = {s.vG2,s.vG3,s.vG4,s.vGT};
    acc = {s.aG2,s.aG3,s.aG4,s.aGT};
    omegas = [repmat(p.omega2,size(s.theta2)); ...
        s.omega3;s.omega4;s.omegaT];
    alphas = [repmat(p.alpha2,size(s.theta2)); ...
        s.alpha3;s.alpha4;s.alphaT];
    energy = zeros(size(s.theta2)); energyRate = energy;
    fdVelocityError = 0; fdAccelerationError = 0; fdPowerError = 0;
    if p.omega2>0
        dt = 2*pi/(numel(s.theta2)*p.omega2);
    end
    for j = 1:4
        energy = energy + masses(j)/2*sum(vel{j}.^2,1) + ...
            inertias(j)/2*omegas(j,:).^2 + masses(j)*p.g*pos{j}(2,:);
        energyRate = energyRate + masses(j)*sum(vel{j}.*acc{j},1) + ...
            inertias(j)*omegas(j,:).*alphas(j,:) + ...
            masses(j)*p.g*vel{j}(2,:);
        if p.omega2>0
            [dv,da] = periodicDerivatives(pos{j},dt);
            fdVelocityError = max(fdVelocityError,max(vectorLength(dv-vel{j})));
            fdAccelerationError = max(fdAccelerationError, ...
                max(vectorLength(da-acc{j})));
        end
    end
    v.energy_J = energy;
    v.energyRate_W = energyRate;
    v.powerBalanceError_W = max(abs(f.power_W-energyRate));
    if p.omega2>0
        de = periodicDerivatives(energy,dt);
        fdPowerError = max(abs(de-f.power_W));
    end
    v.finiteDifferenceVelocityError_m_s = fdVelocityError;
    v.finiteDifferenceAccelerationError_m_s2 = fdAccelerationError;
    v.finiteDifferencePowerError_W = fdPowerError;
    v.derivativeCheckApplicable = p.omega2>0;
    assert(v.loopClosure_m<1e-10*p.lengthScale, ...
        'MonkeyTail:Validation','Loop closure failed.');
    assert(v.newtonEulerResidual_N<1e-9*p.forceScale, ...
        'MonkeyTail:Validation','Newton-Euler residual failed.');
    assert(v.reactionTangency_N<1e-9*p.forceScale, ...
        'MonkeyTail:Validation','Slider reaction is not perpendicular.');
    powerScale = max([p.forceScale*p.lengthScale*p.omega2, ...
        max(abs(f.power_W)),1e-12]);
    assert(v.powerBalanceError_W<1e-9*powerScale, ...
        'MonkeyTail:Validation','Analytical power balance failed.');
    if p.omega2>0
        assert(fdVelocityError<1e-4*p.lengthScale*p.omega2 && ...
            fdAccelerationError<1e-3*p.lengthScale*p.omega2^2 && ...
            fdPowerError<1e-3*powerScale,'MonkeyTail:Validation', ...
            'Finite-difference check failed; increase numSamples.');
    end
    v.passed = true;
end

function [first,second] = periodicDerivatives(x,dt)
    xm2 = circshift(x,2,2); xm1 = circshift(x,1,2);
    xp1 = circshift(x,-1,2); xp2 = circshift(x,-2,2);
    first = (xm2-8*xm1+8*xp1-xp2)/(12*dt);
    second = (-xp2+16*xp1-30*x+16*xm1-xm2)/(12*dt^2);
end

function t = turningPoints(s,p)
    q = p; q.omega2 = 1;
    sample = kinematics(s.theta2,q);
    angles = [s.theta2,2*pi];
    speeds = [sample.omega4,sample.omega4(1)];
    roots = zeros(1,numel(s.theta2));
    numRoots = 0;
    for k = 1:numel(s.theta2)
        if speeds(k)==0
            numRoots = numRoots+1;
            roots(numRoots) = angles(k);
        elseif speeds(k)*speeds(k+1)<0
            numRoots = numRoots+1;
            roots(numRoots) = fzero(@(a) rockerSpeed(a,q), ...
                angles(k:k+1)); 
        end
    end
    roots = roots(1:numRoots);
    roots = sort(mod(roots,2*pi));
    if ~isempty(roots)
        roots = roots([true,diff(roots)>1e-9]);
    end
    poses = kinematics(roots,p);
    t.CrankAngle_deg = rad2deg(roots(:));
    t.C_x_in = poses.C(1,:)'/p.inch;
    t.C_y_in = poses.C(2,:)'/p.inch;
    t.P_x_in = poses.P(1,:)'/p.inch;
    t.P_y_in = poses.P(2,:)'/p.inch;
    t.BarAngle_deg = rad2deg(poses.thetaT(:));
end

function w = rockerSpeed(a,p)
    q = kinematics(a,p); w = q.omega4;
end

function t = forceSummary(theta,f)
    forces = [vectorLength(f.F_A);vectorLength(f.F_B);vectorLength(f.F_C); ...
        vectorLength(f.R_O2);vectorLength(f.R_O4);abs(f.R_O5)];
    [peak,index] = max(forces,[],2);
    t.Force = {'F_A';'F_B';'F_C';'R_O2';'R_O4';'R_O5'};
    t.Peak_N = peak;
    t.Peak_lbf = peak/4.4482216152605;
    t.CrankAngle_deg = rad2deg(theta(index(:)))';
end

function printSummary(r)
    [peakForce,index] = max(r.forceSummary.Peak_N);
    fprintf('Speed: %.1f rpm; payload: %.1f g\n', ...
        r.settings.rpm,r.settings.payload_g);
    fprintf('Peak torque: %.4f N*m; RMS torque: %.4f N*m\n', ...
        r.dynamics.tauPeak_Nm,r.dynamics.tauRMS_Nm);
    fprintf('Maximum reaction force: %.3f N (%s)\n', ...
        peakForce,r.forceSummary.Force{index});
    fprintf('Model checks passed.\n');
end

function fig = plotResults(r)
    s = r.motion; f = r.dynamics; p = r.model;
    a = rad2deg(s.theta2);
    fig = figure('Name','MonkeyTail analysis','Color','w', ...
        'Visible',r.settings.figureVisible,'Position',[80,80,1250,850]);
    subplot(2,3,1); plot(s.C(1,:)/p.inch,s.C(2,:)/p.inch,'LineWidth',1.5);
    hold on; plot(s.P(1,:)/p.inch,s.P(2,:)/p.inch,'LineWidth',1.5);
    axis equal; grid on; xlabel('x [in]'); ylabel('y [in]');
    legend('C: attachment','P: payload COM','Location','best'); title('Trajectories');
    subplot(2,3,2); plot(a,rad2deg([s.theta4;s.thetaT])','LineWidth',1.2);
    grid on; xlabel('Crank angle [deg]'); ylabel('Angle [deg]');
    legend('Rocker','Sliding bar','Location','best'); title('Orientation');
    subplot(2,3,3); plot(a,[vectorLength(s.vC);vectorLength(s.vP)]'/p.inch,'LineWidth',1.2);
    grid on; xlabel('Crank angle [deg]'); ylabel('Speed [in/s]');
    legend('C','P','Location','best'); title('Attachment and payload speed');
    subplot(2,3,4); plot(a,[vectorLength(s.aC);vectorLength(s.aP)]'/p.inch,'LineWidth',1.2);
    grid on; xlabel('Crank angle [deg]'); ylabel('Acceleration [in/s^2]');
    legend('C','P','Location','best'); title('Attachment and payload acceleration');
    subplot(2,3,5); plot(a,[vectorLength(f.F_A);vectorLength(f.F_B);vectorLength(f.F_C); ...
        vectorLength(f.R_O2);vectorLength(f.R_O4);abs(f.R_O5)]','LineWidth',1.2);
    grid on; xlabel('Crank angle [deg]'); ylabel('Force magnitude [N]');
    legend('F_A','F_B','F_C','R_{O2}','R_{O4}','|R_{O5}|','Location','best');
    title('Joint and ground reactions');
    subplot(2,3,6); plot(a,f.tau_Nm,'LineWidth',1.5);
    hold on; plot([a(1),a(end)],[0,0],':k');
    grid on; xlabel('Crank angle [deg]'); ylabel('Crank torque [N m]');
    title('Signed torque (negative requires braking)');
    sgtitle(sprintf('MonkeyTail | %.3g rpm | %.3g g payload', ...
        r.settings.rpm,r.settings.payload_g));
end

function y = perp(x)
    y = [-x(2,:);x(1,:)];
end

function row = momentRow(r)
    row = [-r(2),r(1)];
end

function z = cross2(a,b)
    z = a(1,:).*b(2,:)-a(2,:).*b(1,:);
end

function y = vectorLength(x)
    y = sqrt(sum(x.^2,1));
end

function checkNumbers(x,count,minimum,allowEqual)
    assert(isnumeric(x) && isreal(x) && isvector(x) && ...
        numel(x)==count && all(isfinite(x(:))), ...
        'Expected %d finite real numeric values.', count);
    if allowEqual
        assert(all(x(:)>=minimum),'Values must be at least %g.',minimum);
    else
        assert(all(x(:)>minimum),'Values must be greater than %g.',minimum);
    end
end
