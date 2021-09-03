%% Copyright Information
% MIT License
% 
% Copyright (c) 2021 Open Drilling Project
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% MADSim (last updated with MATLAB R2021a, no additional toolboxes should
%          be required to run)
%   M   - Mechanics
%   A   - And
%   D   - Dynamics
%   Sim - Simulator
%
%   This is source code based on the BHA modeling developed in-house at
%   Scientific Drilling International, Inc (SDI). The provided model is a
%   simplified version of SDI's proprietary code, and is offered as a means
%   to aid in the research and education of future engineers. The model is
%   intended to calculate the static deflection and mechanical loading of a
%   BHA (bent motor or push-the-bit RSS) in an inclined wellbore. The model
%   , as it is provided, accounts for:
%
%       - Geometric nonlinearity (large displacement, small strain)
%       - Fully coupled flexibility (lateral-torsional, axial-lateral)
%       - Automatic determination of frictional wellbore contact points
%
%   The primary beneift of this source code is to offer insight into the
%   implemention of:
%       - Nonlinear finite elements for modeling BHA/drill string mechanics
%       - Wellbore contact via a penalty formulation, within a nonlinear
%       finite element model
%
%   Key featrues that have been removed from the model:
%       - 2D/3D wellbores (i.e. wellbore curvature)
%       - Variations in wellbore OD
%       - Buckling algorithms
%       - Directional estimation algoirhtms
%       - Dynamics (linearized or non-linear)
%
%   Implementation of features that have been removed can be added
%   following approaches outlined in the references
%
%   References
%   1)  J.K. Wilson, "Nonlinear Drillstring Modeling with Applications to
%       Induced Vibrations in Unconventional Horizontal Wells".
%       Ph.D. Dissertation, Texas A&M University. 2017
%   2)  J.K. Wilson, "Field Validation of a New Bottomhole-Assembly Model 
%       for Unconventional Shale Plays".
%       SPE-191780-PA. SPE Drilling and Completion. 2019
%   3)  Heisig, 1995, "Postbuckling Analysis of Drillstrings Using the
%       Finite-Element Method".
%       PD-Vol. 65, Drilling Technology, ASME 1995

%% Data Input and Calculation Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    format compact;
    close all;
    clear variables;

    % Run Parameters
    inc                 = 90;       % Wellbore Inclination [°]
    WOB                 = 35;       % Weight on Bit [klbf]
    TOB                 = 5000;     % Torque on Bit [ft-lbf]
    MW                  = 10;       % Mud Density [ppg]
    DB                  = 8.50;     % Wellbore Diameter [in]
    muw                 = 0.0;      % Wellbore Friction Coef. []
    BendAngle           = 1.75;     % Motor Bend Angle [deg]
    STEERf              = 2500;     % RSS Steering Force [lbf]   
    TFO                 = 90;       % Tool Face Orientation [deg]

    % Model Constraints
    eps                 = 1.0e-6;       % Convergence criteria for static iteration
    kw                  = 1.0e8;        % Wall Stiffness
    z                   = 8;            % Number of Intergration segments - Simpson's 1/3 Rule
    sd                  = 1000;         % Soft Spring Stabilization - Spring Divisor
    ks                  = 6000;         % Initial Soft-Spring Stiffness
    src                 = 6;            % Spring Implementation Coefficient 
                                        %(src = 2 for removal of springs on second iteration)
    tc                  = 20;           % Iteration parameter for inner convergence loop
    looplimit           = 300;          % Maximum Number of Static Loops

    % Specifying component index numbers
    global indRSSap;    indRSSap    = 13;       % RSS steering pad
    
%% Importing Data
        [Data,Stab,InputFilePath,InputFileName] = MADSimDataImport();
        
%% Data Organization and Variable Initialization
        noe         = size(Data,1)-1;       % Number of elements in the model
        non         = noe + 1;              % Number of nodes in the model
        nStab       = size(Stab,1);         % Number of stabilizers in the model
    
        Data(1,1)   = 1-MW/65.5;            % Storing Buoyancy Factor
        Data(1,5)   = MW;                   % Storing Mud Weight
        Data(1,6)   = WOB;                  % Storing WOB
        Data(1,7)   = TOB;                  % Storing TOB
        Data(1,10)  = TFO;                  % Storing TFO
        Data(1,11)  = BendAngle;            % Storing Bend Angle
        
        BendEle     = Data(2,11);           % Storing number of elements between bit and bend
        dof         = non*6;                % Degrees of freedom within the model
        TFO         = TFO*pi/180;           % Converting TFO from ° to radians
    
    % Populating array with Model Contstraints
        BCs         = zeros(1,8);
        BCs(1,1)    = eps;
        BCs(1,2)    = kw;
        BCs(1,3)    = z;
        BCs(1,4)    = sd;
        BCs(1,5)    = ks;
        BCs(1,6)    = src;
        BCs(1,7)    = tc;
        BCs(1,8)    = looplimit;

    
%% Bit Offset Calculation for simple Bend Angle Model
    % Dummy displacement vector correlating to the geometry of the motor's bend angle.
    % this value is added to the displacement vector in order to simulate the
    % appropirate clearance between the wellbore and the motor bend
    Uxo = zeros(non,1);
    
    % Determining the bit-to-bend length (Lbb)
    s = 0;
    for i = 1:BendEle
        s = s + Data(i+1,3);
    end
    Lbb = s;
    % Establishing dummy geometry of motor bend (Uxo)
    w = non;
    for i = 1:BendEle
        if i == 1
            s = Lbb;
        else
            s = s - Data(i,3);
        end
        Uxo(w) = s*sin(BendAngle*pi/180);
        w = w - 1;
    end

%% External Force Vector
    % Determining gravitational force components
    [Q] = ExternalLoad(Data,z,inc,0);
        % Adding WOB and TOB to the model
        Q(6*(noe+1)-1) = Q(6*(noe+1)-1) - WOB*1000;
        Q(6*(noe+1)) = Q(6*(noe+1)) - TOB;
    
        % Adding steer force for a push-the-bit RSS
        for i = 1:nStab
            if Stab(i,5) == indRSSap
                nap = non + 1 - Stab(i,1);
                Q((nap-1)*6 + 1) = Q((nap-1)*6 + 1) + STEERf*sin(-TFO);
                Q((nap-1)*6 + 3) = Q((nap-1)*6 + 3) + STEERf*cos(-TFO);
            else
            end
        end


%% Iterative Static Solution
    % Calculates the static delfection of the BHA/Drill string using a
    % Newton-Rhapson iteration
    
    % Starting timing with "tic"
    tic;
    
        [U,out] = MADSimStatic(Data,BCs,Stab,Q,Uxo,DB/12,muw);
        
    % Ending timing with "toc". Output should display elapsed time in
    % Matlab command window
    toc;

    
%% Determining loading from static delfection
    [ExcDataMat] = Loading(U,Uxo,Data,z,kw,Stab,inc,DB/12);

%% Plotting Results
        close all;
        % Create radial clearance vector
        rclear = zeros(non,1);
        for i = 1:non-1
            rclear(i) = (DB - Data(i+1,1))/2;
        end
        
        for i = 1:nStab
            rclear(Stab(i,1)) = (DB - Stab(i,7))/2;
        end
        
        figure(1);
        string = 'Axial/Torsional BHA Deflection';
        set(gcf,'Name',string,'NumberTitle','off');
        set(gcf,'pos',[100 100 1600 800]);
        
        % Axial Deflection
        subplot(4,3,[1;3]);
        plot(ExcDataMat (:,1),ExcDataMat(:,6),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('Axial Deflection (in)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.4f';
        
        % Torsional Displacement
        subplot(4,3,[4;6]);
        plot(ExcDataMat (:,1),ExcDataMat (:,7),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('Torsion (°)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.4f';
        
        % Axial Force
        subplot(4,3,[7;9]);
        plot(ExcDataMat (:,1),ExcDataMat(:,17),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('Axial Force (lbf)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.4f';
        
        % Torque
        subplot(4,3,[10;12]);
        plot(ExcDataMat (:,1),ExcDataMat (:,18),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('Torque (ft-lbf)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.4f';
        
        
        figure(2);
        string = 'High-Side/Right-Side BHA Deflection';
        set(gcf,'Name',string,'NumberTitle','off');
        set(gcf,'pos',[100 100 1600 800]);
        
        % High-Side Deflection
        subplot(4,3,[1;3]);
        plot(ExcDataMat (:,1),ExcDataMat (:,2),'Linewidth',1.5);
        hold on;
        stairs(ExcDataMat (:,1),rclear(:),'-k');
        hold on;
        stairs(ExcDataMat (:,1),-rclear(:),'-k');
        xlabel('Distance from Bit (ft)');
        ylabel('High-Side Deflection (in)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.4f';
        legend('Deflection','Radial Clearance Limit');
        
        % High-Side Sag
        subplot(4,3,[4;6]);
        plot(ExcDataMat (:,1),ExcDataMat (:,3),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('High-Side Sag (°)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.4f';
        
        % Right-Side Deflection
        subplot(4,3,[7;9]);
        plot(ExcDataMat (:,1),ExcDataMat (:,4),'Linewidth',1.5);
        hold on;
        stairs(ExcDataMat (:,1),rclear(:),'-k');
        hold on;
        stairs(ExcDataMat (:,1),-rclear(:),'-k');
        xlabel('Distance from Bit (ft)');
        ylabel('Right-Side Deflection (in)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.4f';
        legend('Deflection','Radial Clearance Limit');
        
        % Right-Side Sag
        subplot(4,3,[10;12]);
        plot(ExcDataMat (:,1),ExcDataMat (:,5),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('Right-Side Sag (°)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.4f';
        
        
        % Plotting Loading************************************************
        figure(3);
        string = 'High-Side/Right-Side BHA Loading';
        set(gcf,'Name',string,'NumberTitle','off');
        set(gcf,'pos',[100 100 1600 900]);
        
        % High-Side Bending
        subplot(6,3,[1;3]);
        plot(ExcDataMat (:,1),ExcDataMat (:,11),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('HS Moment (ft-lbf)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.0f';
        
        % High-Side Contact Force
        subplot(6,3,[4;6]);
        plot(ExcDataMat (:,1),ExcDataMat (:,8),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('HS Contact Force (lbf)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.0f';
        
        % High-Side Shear Force
        subplot(6,3,[7;9]);
        plot(ExcDataMat (:,1),ExcDataMat (:,14),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('HS Shear Force (lbf)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.2f';
        
        % Right-Side Bending
        subplot(6,3,[10;12]);
        plot(ExcDataMat (:,1),ExcDataMat (:,12),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('RS Moment (ft-lbf)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.0f';
        
        % Right-Side Contact
        subplot(6,3,[13;15]);
        plot(ExcDataMat (:,1),ExcDataMat (:,9),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('RS Contact Force (lbf)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.0f';
        
        % Right-Side Shear
        subplot(6,3,[16;18]);
        plot(ExcDataMat (:,1),ExcDataMat (:,15),'Linewidth',1.5);
        xlabel('Distance from Bit (ft)');
        ylabel('RS Shear Force (lbf)');
        xlim([0 max(ExcDataMat (:,1))]);
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%,.2f';

%% Functions
function [Data,Stab,InputFilePath,InputFileName] = MADSimDataImport(varargin)
% MADSimDataImport Imports data from excel file
    % Data          = array containing BHA/ drillstring dimensions and material
    %                 properties
    % Stab          = Array containing stabilizer information
    % InputFilePath = string containing the folder location of the input
    %                 file
    % InputFileName = string containing the input file name

    % Select import folder
    [InputFileName,InputFilePath] = uigetfile({'*.xlsx;*.xls;*.xlsm','Excel Files (*.xlsx;*.xls;*.xlsm)'},'Select BHA File');
    InStr = [InputFilePath InputFileName];
    
    % Reading number of elements (noe), number of stabilizers (nStab), and
    % number of SUrvey points (SurvPoints)
    noe = xlsread(InStr,'MatLabData','C1','basic');
    nStab = xlsread(InStr,'MatLabData','D1','basic');
    
    % Importing BHA Element Data
    string = ['A1:K' num2str(noe+1)];
    Data = xlsread(InStr,'MatLabData',string,'basic');
    
    % Importing Stabilizer Data
    if nStab > 0
        string = ['Q1:W' num2str(nStab)];
        Stab = xlsread(InStr,'MatLabData',string,'basic');
    else
        Stab = 1;
    end
end
function [Q] = ExternalLoad(Data,z,inc,azi)
%ExternalLoad returns the gravitational force vector acting on the
%BHA/drill string assembly
%   Q   = Global gravitational force vector based on constant inclination
%         and azimuth
        inc = inc*pi/180;
        azi = azi*pi/180;
        noe = size(Data,1) - 1;
        BF = Data(1,1);
        Q = zeros(6*(noe+1),1);
        w = noe;
         
        % Determing gravitational parameters
            NormMat = [1,0,0;0,cos(inc),sin(inc);0,-sin(inc),cos(inc)].'*...
                      [cos(pi/2-azi),sin(pi/2-azi),0;-sin(pi/2-azi),cos(pi/2-azi),0;0,0,1].'*...
                      [1,0,0;0,1,0;0,0,1];
            v1 = [NormMat(1,1),NormMat(1,2),NormMat(1,3)];
            v2 = [NormMat(2,1),NormMat(2,2),NormMat(2,3)];
            v3 = [NormMat(3,1),NormMat(3,2),NormMat(3,3)];

            g1 = dot(v1,[0,0,1]);
            g2 = dot(v2,[0,0,1]);
            g3 = dot(v3,[0,0,1]);
            
        for n = 1:noe
            % Zeroing elemental force vector
                Qele = zeros(12,1);
            
            % Extracting elemental data
                E = Data(w+1,5)*144;
                I = Data(w+1,7)/12^4;
                OD = Data(w+1,1)/12;
                nu = Data(w+1,6); 
                G = E/(2*(1+nu));
                ID = Data(w+1,2)/12;
                A = pi/4*(OD^2-ID^2);
                q = Data(w+1,4)*BF;
                L = Data(w+1,3);

            % Integrating elemental force vectors
                h = 1/z;
                eta = 0;
                Ks = 6*(1+nu)*(1+ID/OD)^2/((7+6*nu)*(1+ID/OD)^2+4*(5+3*nu)*(ID/OD)^2); % Shear correction factor
                Lambda = 12*E*I/(Ks*G*A*L^2);
                
                for i = 1:(z+1)
                    % Required shape functions
                        N1 = 1/(1-Lambda)*(2*eta^3-3*eta^2+Lambda*eta+1-Lambda);
                        N2 = L/(1-Lambda)*(eta^3+eta^2*(1/2*Lambda-2)+eta*(1-1/2*Lambda));
                        N3 = 1/(1-Lambda)*(-2*eta^3+3*eta^2-Lambda*eta);
                        N4 = L/(1-Lambda)*(eta^3-eta^2*(1+1/2*Lambda)+eta*Lambda/2);
                        N9 = 1-eta;
                        N10 = eta;

                        h1 = [N1;N2;0;0;0;0;N3;N4;0;0;0;0];
                        h3 = [0;0;N1;N2;0;0;0;0;N3;N4;0;0];
                        h5 = [0;0;0;0;N9;0;0;0;0;0;N10;0];
                    
                    % Integrating using Simpson's 1/3 rule
                        if i == 1 || i == z+1
                            Qele = Qele + h/3*L*...
                                q*(g1*h1+g2*h3+g3*h5);
                        elseif mod(i,2)==0
                            Qele = Qele + 4*h/3*L*...
                                q*(g1*h1+g2*h3+g3*h5);
                        elseif mod(i,2)~=0
                            Qele = Qele + 2*h/3*L*...
                                q*(g1*h1+g2*h3+g3*h5);
                        end

                    eta = eta + h;
                end

                for i = 1:12
                    Q((n-1)*6+i) = Q((n-1)*6+i) + Qele(i);
                end
        w = w-1;     
        end
end
function [U,out] = MADSimStatic(Data,BCs,Stab,Q,Uxo,DB,muw)
%Quasi-Static Calculation of the Drillstring assembly
%   Determines the deflection of the BHA/Drillstring under "snap-shot"
%   loading, including frictional drag. The "Snap-Shot" assumes frictional
%   forces resulting from axial movement and rotation
%   
%   U   = global displacement vector
%   out = Convergence indicator. If = 1, convergance has not been reached

    non         = size(Data,1);
    dof         = 6*non;
    Uin         = zeros(dof,1);
    U           = Uin;
    eps         = BCs(1,1);
    z           = BCs(1,3);
    sd          = BCs(1,4);
    ks          = BCs(1,5);
    src         = BCs(1,6);
    tc          = BCs(1,7);
    looplimit   = BCs(1,8);
    out         = 0;
    
    tt = 1; 
    error = 1;
    counter = 1;
    sscheck = 1;
    
    % Outer convergence loop (counter = tt)
    while error >= eps && counter <= looplimit
        for i = 1:dof
            Uin(i) = U(i);
        end
        convc = 1;
        t = 1;
        
        % Inner convergence loop (counter = t)
        while convc >= eps  && t <tc && counter <= looplimit
            % Calculating global stiffness matrix, internal force vector,
            % and wall contact force vector (K, Fg, Fw)
                K   = KMatrix(BCs,Data,Uin,Stab,Uxo,DB,muw);
                Fg  = FgVector(Data,Uin,z);
                Fw  = FwVector(BCs,Data,Uin,Stab,Uxo,DB,muw);
            
            % Adding soft spring stabilization. Soft-springs are removed
            % from the system after the external loop (tt) reaches the
            % specified limit (src)
            if sd*tt <= ks && tt < src
                for i = 1:non
                    K(6*(i-1)+1,6*(i-1)+1) = K(6*(i-1)+1,6*(i-1)+1) + ks - (sd*tt);
                    K(6*(i-1)+3,6*(i-1)+3) = K(6*(i-1)+3,6*(i-1)+3) + ks - (sd*tt);
                    K(6*(i-1)+5,6*(i-1)+5) = K(6*(i-1)+5,6*(i-1)+5) + ks - (sd*tt);
                end
            else
                sscheck = 0;
                disp('Soft Spring Detached');
            end
            
            % Residual Force Vector
                R = Q - Fg - Fw;
            % delU calculation
                delUin = linsolve(K,R);
            % Adding incremental displacement to total displacement
                Uin = Uin + delUin;
            
            counter = counter + 1;

            %Inner Loop error and convergence check
            for i = 1:non
                if i ==1
                    dispcin = abs(delUin(i))^2;
                    dispccin = abs(Uin(i))^2;
                else
                    dispcin = dispcin + abs(delUin(i))^2;
                    dispccin = dispccin + abs(Uin(i))^2;
                end
            end
            convc = sqrt(dispcin/dispccin);
            fprintf('Internal Error = ');
            disp(convc);
            t = t + 1;
        end
        
        delU = Uin - U;
        U = U + delU;
        
        %Outer Loop error and convergence check
        for i =1:dof
            if i ==1
                dispc = abs(delU(i))^2;
                dispcc = abs(U(i))^2;
            else
                dispc = dispc + abs(delU(i))^2;
                dispcc = dispcc + abs(U(i))^2;
            end
        end

        error = sqrt(dispc/dispcc);
        
        fprintf('External Error = ');
        disp(error);
        
        if error < eps && sscheck == 1
            error = 2*eps;
        else
        end
        
        tt = tt +1;
    end
    
    if convc >= eps
        % Convergence check. If the solution has not converged within the
        % specified tolerance (eps), "out" is set to 1
        out = 1;
    else
    end
end

    function [K] = KMatrix(BCs,Data,U,Stab,Uxo,DB,muw)
    %KMatrix Constructs the stiffness matrix of the drillstring/BHA
    %   K = Stiffness matrix
    
            % Extracting initial parameters and initializing data
                non = size(Data,1);             % Number of nodes
                noe = non - 1;                  % Number of elements
                nStab = size(Stab,1);           % Number of Stabilizers
                TFO = Data(1,10)*pi/180;        % Tool Face [Radians]
                BendAngle = Data(1,11)*pi/180;  % Bend Angle [radians]
                BendEle = Data(2,11);           % Number of elements below the motor's bend
                kw = BCs(1,2);                  % Wall stiffness
                z = BCs(1,3);                   % Integration intervals

                Uele = zeros(12,1);             % Initializing Elemental displacement
                K = zeros(6*(noe+1),6*(noe+1)); % Initializing total stiffness matrix
                Kg = K;                         % Initializing geometric/beam stiffness matrix
                Kw = K;                         % Initializing wall contact stiffness matrix

            % Looping through elements
            n = 1;    
            for w = noe:-1:1
                % Zeroing elemental stiffnes matrices
                    dFgdu = zeros(12,12);
                    dFwdu = zeros(12,12);
                % Extracting elemental data
                    E = Data(w+1,5)*144;
                    L = Data(w+1,3);
                    I = Data(w+1,7)/12^4;
                    OD = Data(w+1,1)/12;
                    ID = Data(w+1,2)/12;
                    J = 2*I;
                    nu = Data(w+1,6); 
                    A = pi/4*(OD^2-ID^2);
                    G = E/(2*(1+nu));
                % Extracting elemental displacement determined in the
                % previous delU calculation
                    for i = 1:12
                        Uele(i) = U((n-1)*6+i);
                    end
                
                % Extracing local displacements for wall contact stiffness
                % calculations
                    x1 = Uele(1) + sin(-TFO)*Uxo(n);
                    x2 = Uele(7) + sin(-TFO)*Uxo(n+1);
                    y1 = Uele(3) + cos(-TFO)*Uxo(n);
                    y2 = Uele(9) + cos(-TFO)*Uxo(n+1);
                    phix1 = Uele(2);
                    phix2 = Uele(8);
                    phiy1 = Uele(4);
                    phiy2 = Uele(10);
                    
                    
                % Applying element transformation for elements below the
                % motor's bend
                    if n >= non-BendEle
                        [Tbend] = TransformBend(BendAngle);
                        Uele = Tbend*Uele;
                    else
                    end

                % Shear correction factor    
                    Ks = 6*(1+nu)*(1+(ID/OD)^2)^2/((7+6*nu)*(1+(ID/OD)^2)^2+4*(5+3*nu)*(ID/OD)^2);
                    Lambda = 12*E*I/(Ks*G*A*L^2);

                % Integrating elemental geometric/beam stiffness
                    h = 1/z;
                    eta = 0;
                    for i = 1:(z+1)
                        % Required shape functions
                            N5 = 6/(L*(1-Lambda))*(eta^2-eta);
                            N6 = 1/(1-Lambda)*(3*eta^2+eta*(Lambda-4)+1-Lambda);
                            N7 = 6/(L*(1-Lambda))*(eta-eta^2); 
                            N8 = 1/(1-Lambda)*(3*eta^2-eta*(Lambda+2));

                            DN1 = 1/(L*(1-Lambda))*(6*eta^2-6*eta+Lambda);
                            DN2 = 1/(1-Lambda)*(3*eta^2+2*eta*(1/2*Lambda-2)+1-1/2*Lambda);
                            DN3 = 1/(L*(1-Lambda))*(-6*eta^2+6*eta-Lambda);
                            DN4 = 1/(1-Lambda)*(3*eta^2-2*eta*(1+1/2*Lambda)+1/2*Lambda);
                            DN5 = 6/(L^2*(1-Lambda))*(2*eta-1);
                            DN6 = 1/(L*(1-Lambda))*(6*eta+Lambda-4);
                            DN7 = 6/(L^2*(1-Lambda))*(1-2*eta); 
                            DN8 = 1/(L*(1-Lambda))*(6*eta-Lambda-2);
                            DN9 = -1/L;
                            DN10 = 1/L;

                            h2 = [N5;N6;0;0;0;0;N7;N8;0;0;0;0];
                            h4 = [0;0;N5;N6;0;0;0;0;N7;N8;0;0];

                            Dh1 = [DN1;DN2;0;0;0;0;DN3;DN4;0;0;0;0];
                            Dh2 = [DN5;DN6;0;0;0;0;DN7;DN8;0;0;0;0];
                            Dh3 = [0;0;DN1;DN2;0;0;0;0;DN3;DN4;0;0];
                            Dh4 = [0;0;DN5;DN6;0;0;0;0;DN7;DN8;0;0];
                            Dh5 = [0;0;0;0;DN9;0;0;0;0;0;DN10;0];
                            Dh6 = [0;0;0;0;0;DN9;0;0;0;0;0;DN10];

                            Dh3til = Dh3;
                            Dh5til = Dh5;
                            Dh6til = Dh6;
                            h4til = h4;
                            h43til = h4 - Dh3;
                            h12til = Dh1 - h2;

                            H1 = (Dh5til.'*Uele) + 1/2*(h2.'*Uele)^2 + 1/2*(h4til.'*Uele)^2 + (Dh1.'*Uele)*(h12til.'*Uele) - (h43til.'*Uele)*(Dh3til.'*Uele);
                            H2 = (Dh6til.'*Uele) - (Dh4.'*Uele)*(h2.'*Uele);

                            DH1 = (Dh5til.') + (h2.'*Uele)*h2.' + (h4til.'*Uele)*h4til.' + (Dh1.'*Uele)*(h12til.') + (h12til.'*Uele)*(Dh1.') - (h43til.'*Uele)*(Dh3til.') - (Dh3til.'*Uele)*(h43til.');
                            DH2 = Dh6til.' - (h2.'*Uele)*Dh4.' - (Dh4.'*Uele)*h2.';

                        % Integrating using Simpsons's 1/3 rule
                            if i == 1 || i == z + 1
                                dFgdu = dFgdu + h/3*L*...
                                    (E*I*(Dh2*(Dh2.') + Dh4*(Dh4.'))...
                                    +G*A*Ks*(h43til*(h43til.') + h12til*(h12til.'))...
                                    +E*A*(H1*(h2*(h2.') + Dh1*(h12til.') + (h12til)*(Dh1.') + h4til*(h4til.') - Dh3til*(h43til.') - (h43til)*(Dh3til.'))...
                                        +(Dh5til + h2*(h2.'*Uele) + Dh1*(h12til.'*Uele) + (h12til)*(Dh1.'*Uele) + h4til*(h4til.'*Uele) - Dh3til*(h43til.'*Uele) - (h43til)*(Dh3til.'*Uele))*DH1)...
                                    +G*J*(H2*(-h2*(Dh4.') - Dh4*(h2.'))...
                                        +(Dh6til - h2*(Dh4.'*Uele) - Dh4*(h2.'*Uele))*DH2));
                            elseif mod(i,2)==0
                                dFgdu = dFgdu + 4*h/3*L*...
                                    (E*I*(Dh2*(Dh2.') + Dh4*(Dh4.'))...
                                    +G*A*Ks*(h43til*(h43til.') + h12til*(h12til.'))...
                                    +E*A*(H1*(h2*(h2.') + Dh1*(h12til.') + (h12til)*(Dh1.') + h4til*(h4til.') - Dh3til*(h43til.') - (h43til)*(Dh3til.'))...
                                        +(Dh5til + h2*(h2.'*Uele) + Dh1*(h12til.'*Uele) + (h12til)*(Dh1.'*Uele) + h4til*(h4til.'*Uele) - Dh3til*(h43til.'*Uele) - (h43til)*(Dh3til.'*Uele))*DH1)...
                                    +G*J*(H2*(-h2*(Dh4.') - Dh4*(h2.'))...
                                        +(Dh6til - h2*(Dh4.'*Uele) - Dh4*(h2.'*Uele))*DH2));
                            elseif mod(i,2)~=0
                                dFgdu = dFgdu + 2*h/3*L*...
                                    (E*I*(Dh2*(Dh2.') + Dh4*(Dh4.'))...
                                    +G*A*Ks*(h43til*(h43til.') + h12til*(h12til.'))...
                                    +E*A*(H1*(h2*(h2.') + Dh1*(h12til.') + (h12til)*(Dh1.') + h4til*(h4til.') - Dh3til*(h43til.') - (h43til)*(Dh3til.'))...
                                        +(Dh5til + h2*(h2.'*Uele) + Dh1*(h12til.'*Uele) + (h12til)*(Dh1.'*Uele) + h4til*(h4til.'*Uele) - Dh3til*(h43til.'*Uele) - (h43til)*(Dh3til.'*Uele))*DH1)...
                                    +G*J*(H2*(-h2*(Dh4.') - Dh4*(h2.'))...
                                        +(Dh6til - h2*(Dh4.'*Uele) - Dh4*(h2.'*Uele))*DH2));
                            end
                        eta = eta + h;
                    end
                    
                % Determining Wall Contact Stiffness
                    nodenumber2 = w;
                    nodenumber1 = w+1;

                    % Determining allowable clearance between element Node 
                    % 1 and Wellbore wall
                        % if node 1 is a stabilizer:
                            rcheck1 = 0;
                            if nStab > 0
                                for i = 1:nStab
                                    if Stab(i,1) == nodenumber1
                                        r1 = (DB - Stab(i,7)/12)/2;     % Clearance between stabilizer and wellbore
                                        OD1 = DB - r1*2;
                                        rcheck1 = 1;
                                    else
                                    end
                                end
                            else
                            end
                        % If Node 1 is NOT a stabilizer
                            if rcheck1 == 0
                                r1 = (DB-OD)/2;     % Clearance between element OD and wellbore
                                OD1 = OD;
                            else
                            end
                    % Determining allowable clearance between element Node 
                    % 2 and Wellbore wall
                        % if node 2 is a stabilizer:
                            rcheck2 = 0;
                            if nStab > 0
                                for i = 1:nStab
                                    if Stab(i,1) == nodenumber2
                                        r2 = (DB - Stab(i,7)/12)/2;     % Clearance between stabilizer and wellbore
                                        OD2 = DB - r2*2;
                                        rcheck2 = 1;
                                    else
                                    end
                                end
                            else
                            end
                        % If Node 2 is NOT a stabilizer
                            if rcheck2 == 0
                                r2 = (DB-OD)/2;     % Clearance between element OD and wellbore
                                OD2 = OD;
                            else
                            end
                    % Clearance between current node displacement and
                    % wellbore wall
                        f1 = r1 - sqrt(x1^2+y1^2);
                        f2 = r2 - sqrt(x2^2+y2^2);

                    % Calculating wall contact stiffness based on nodal
                    % clearance
                        % For the 1st elemental node:
                            if f1 < 0
                                dFwdu(1,1) = kw*((sqrt(x1^2+y1^2))^3-r1*y1^2)/(sqrt(x1^2+y1^2))^3 - muw*kw*r1*x1*y1/(sqrt(x1^2+y1^2))^3;
                                dFwdu(1,3) = kw*x1*y1*r1/(sqrt(x1^2+y1^2))^3 + muw*kw*(r1*x1^2-(sqrt(x1^2+y1^2))^3)/(sqrt(x1^2+y1^2))^3;
                                dFwdu(3,1) = kw*x1*y1*r1/(sqrt(x1^2+y1^2))^3 - muw*kw*(r1*y1^2-(sqrt(x1^2+y1^2))^3)/(sqrt(x1^2+y1^2))^3;
                                dFwdu(3,3) = kw*((sqrt(x1^2+y1^2))^3-r1*x1^2)/(sqrt(x1^2+y1^2))^3 + muw*kw*r1*x1*y1/(sqrt(x1^2+y1^2))^3;
                                dFwdu(4,1) = -muw*kw*OD1*x1*phix1/(2*sqrt(x1^2+y1^2));
                                dFwdu(4,2) = muw*kw*OD1/2*(r1-sqrt(x1^2+y1^2));
                                dFwdu(4,3) = -muw*kw*OD1*y1*phix1/(2*sqrt(x1^2+y1^2));
                                dFwdu(5,1) = muw*kw*(((sqrt(x1^2+y1^2))^3-r1*(x1^2+y1^2))*(phix1*y1-(phiy1))+r1*x1*(phix1*y1-(phiy1)*x1))/(sqrt(x1^2+y1^2))^3;
                                dFwdu(5,2) = muw*kw*(sqrt(x1^2+y1^2)-r1)*(y1-(phiy1)*x1)/sqrt(x1^2+y1^2);
                                dFwdu(5,3) = muw*kw*(((sqrt(x1^2+y1^2))^3-r1*(x1^2+y1^2))*(phix1-(phiy1)*x1)+r1*y1*(phix1*y1-(phiy1)*x1))/(sqrt(x1^2+y1^2))^3;
                                dFwdu(5,4) = -muw*kw*x1*(sqrt(x1^2+y1^2)-r1)/sqrt(x1^2+y1^2);

                                dFwdu(6,1) = muw*kw*OD1*x1/(2*sqrt(x1^2+y1^2));
                                dFwdu(6,3) = muw*kw*OD1*y1/(2*sqrt(x1^2+y1^2));
                            else
                            end
                        
                        % For the 2nd elemental node:
                            if f2 < 0
                                dFwdu(7,7) = kw*((sqrt(x2^2+y2^2))^3-r2*y2^2)/(sqrt(x2^2+y2^2))^3 - muw*kw*r2*x2*y2/(sqrt(x2^2+y2^2))^3;
                                dFwdu(7,9) = kw*x2*y2*r2/(sqrt(x2^2+y2^2))^3 + muw*kw*(r2*x2^2-(sqrt(x2^2+y2^2))^3)/(sqrt(x2^2+y2^2))^3;
                                dFwdu(9,7) = kw*x2*y2*r2/(sqrt(x2^2+y2^2))^3 - muw*kw*(r2*y2^2-(sqrt(x2^2+y2^2))^3)/(sqrt(x2^2+y2^2))^3;
                                dFwdu(9,9) = kw*((sqrt(x2^2+y2^2))^3-r2*x2^2)/(sqrt(x2^2+y2^2))^3 + muw*kw*r2*x2*y2/(sqrt(x2^2+y2^2))^3;
                                dFwdu(10,7) = -muw*kw*OD2*x2*phix2/(2*sqrt(x2^2+y2^2));
                                dFwdu(10,8) = muw*kw*OD2/2*(r2-sqrt(x2^2+y2^2));
                                dFwdu(10,9) = -muw*kw*OD2*y2*phix2/(2*sqrt(x2^2+y2^2));
                                dFwdu(11,7) = muw*kw*(((sqrt(x2^2+y2^2))^3-r2*(x2^2+y2^2))*(phix2*y2-(phiy2))+r2*x2*(phix2*y2-(phiy2)*x2))/(sqrt(x2^2+y2^2))^3;
                                dFwdu(11,8) = muw*kw*(sqrt(x2^2+y2^2)-r2)*(y2-(phiy2)*x2)/sqrt(x2^2+y2^2);
                                dFwdu(11,9) = muw*kw*(((sqrt(x2^2+y2^2))^3-r2*(x2^2+y2^2))*(phix2-(phiy2)*x2)+r2*y2*(phix2*y2-(phiy2)*x2))/(sqrt(x2^2+y2^2))^3;
                                dFwdu(11,10) = -muw*kw*x2*(sqrt(x2^2+y2^2)-r2)/sqrt(x2^2+y2^2);

                                dFwdu(12,7) = muw*kw*OD2*x2/(2*sqrt(x2^2+y2^2));
                                dFwdu(12,9) = muw*kw*OD2*y2/(2*sqrt(x2^2+y2^2));
                            else
                            end

                Keleg = dFgdu;
                Kelew = dFwdu;
                
                % Tranforming elemental stiffness for elements below the
                % motor's bend
                    if n >= non-BendEle
                        Keleg = Tbend.'*Keleg*Tbend;
                    else
                    end

                % Boundary conditions at "top" of BHA/drill string
                    if n == 1
                        Keleg(2,2) = Keleg(2,2) + 1e12;
                        Keleg(4,4) = Keleg(4,4) + 1e12;
                        Keleg(5,5) = Keleg(5,5) + 1e12;
                        Keleg(6,6) = Keleg(6,6) + 1e12;
                    else
                    end

                % Assembling elemental stiffness matrices into global
                % stiffness matrix
                    for i = 1:12
                        for j = 1:12
                            Kg((n-1)*6+i,(n-1)*6+j) = Kg((n-1)*6+i,(n-1)*6+j) + Keleg(i,j);
                        end
                    end 
                    for i = 1:12
                        for j = 1:12
                            Kw((n-1)*6+i,(n-1)*6+j) = Kelew(i,j);
                        end
                    end 
                    
                n = n+1;
            end
        % Total Stiffness Matrix
            K = Kg + Kw;
    end
    function [Fg] = FgVector(Data,U,z)
    %FgVector creates the internal force vector of the BHA/Drill string
    %   Fg = global internal force vector
    
        % Extracing initial parameters and initializing data
            non = size(Data,1);             % Number of nodes
            noe = non - 1;                  % Number of elements
            BendEle = Data(2,11);           % Number of elements below the motor's bend
            BendAngle = Data(1,11)*pi/180;  % Bend angle [radians]
            Uele = zeros(12,1);             % Initializing elemental displacement
            Fg = zeros(6*(noe+1),1);        % Initializing element force vector
            
        % Elemental loop
            n = 1;
            for w = noe:-1:1
                % Zeroing elemental force vector
                    Fge = zeros(12,1);
                    
                % Extracing elemental data
                    E = Data(w+1,5)*144;
                    L = Data(w+1,3);
                    I = Data(w+1,7)/12^4;
                    OD = Data(w+1,1)/12;
                    ID = Data(w+1,2)/12;
                    J = 2*I;
                    nu = Data(w+1,6);
                    A = pi/4*(OD^2-ID^2);
                    G = E/(2*(1+nu));
                    
                % Extracing elemental displacements
                    for i = 1:12
                        Uele(i) = U((n-1)*6+i);
                    end
                % Transforming elemental displacement for elements below
                % the bend
                    if n >= non-BendEle
                        [Tbend] = TransformBend(BendAngle);
                        Uele = Tbend*Uele;
                    else
                    end
                
                % Shear correction factor
                    Ks = 6*(1+nu)*(1+(ID/OD)^2)^2/((7+6*nu)*(1+(ID/OD)^2)^2+4*(5+3*nu)*(ID/OD)^2);
                    Lambda = 12*E*I/(Ks*G*A*L^2);
                
                % Integrating elemental force vector
                    h = 1/z;
                    eta = 0;
                    for i = 1:(z+1)
                        % Required shape functions
                            N5 = 6/(L*(1-Lambda))*(eta^2-eta);
                            N6 = 1/(1-Lambda)*(3*eta^2+eta*(Lambda-4)+1-Lambda);
                            N7 = 6/(L*(1-Lambda))*(eta-eta^2); 
                            N8 = 1/(1-Lambda)*(3*eta^2-eta*(Lambda+2));

                            DN1 = 1/(L*(1-Lambda))*(6*eta^2-6*eta+Lambda);
                            DN2 = 1/(1-Lambda)*(3*eta^2+2*eta*(1/2*Lambda-2)+1-1/2*Lambda);
                            DN3 = 1/(L*(1-Lambda))*(-6*eta^2+6*eta-Lambda);
                            DN4 = 1/(1-Lambda)*(3*eta^2-2*eta*(1+1/2*Lambda)+1/2*Lambda);
                            DN5 = 6/(L^2*(1-Lambda))*(2*eta-1);
                            DN6 = 1/(L*(1-Lambda))*(6*eta+Lambda-4);
                            DN7 = 6/(L^2*(1-Lambda))*(1-2*eta); 
                            DN8 = 1/(L*(1-Lambda))*(6*eta-Lambda-2);
                            DN9 = -1/L;
                            DN10 = 1/L;

                            h2 = [N5;N6;0;0;0;0;N7;N8;0;0;0;0];
                            h4 = [0;0;N5;N6;0;0;0;0;N7;N8;0;0];

                            Dh1 = [DN1;DN2;0;0;0;0;DN3;DN4;0;0;0;0];
                            Dh2 = [DN5;DN6;0;0;0;0;DN7;DN8;0;0;0;0];
                            Dh3 = [0;0;DN1;DN2;0;0;0;0;DN3;DN4;0;0];
                            Dh4 = [0;0;DN5;DN6;0;0;0;0;DN7;DN8;0;0];
                            Dh5 = [0;0;0;0;DN9;0;0;0;0;0;DN10;0];
                            Dh6 = [0;0;0;0;0;DN9;0;0;0;0;0;DN10];

                            Dh3til = Dh3;
                            Dh5til = Dh5;
                            Dh6til = Dh6;
                            h4til = h4;

                            h43til = h4 - Dh3;
                            h12til = Dh1 - h2;

                            H1 = (Dh5til.'*Uele) + 1/2*(h2.'*Uele)^2 + 1/2*(h4til.'*Uele)^2 + (Dh1.'*Uele)*(h12til.'*Uele) - (h43til.'*Uele)*(Dh3til.'*Uele);
                            H2 = (Dh6til.'*Uele) - (Dh4.'*Uele)*(h2.'*Uele);
                            
                        % Integrating using Simpson's 1/3 rule
                            if i == 1 || i == z + 1
                                Fge = Fge + h/3*L*...
                                        (E*I*(Dh2*(Dh2.'*Uele) + Dh4*(Dh4.'*Uele))...
                                        +G*A*Ks*(h43til*(h43til.'*Uele) + h12til*(h12til.'*Uele))...
                                        +E*A*H1*(Dh5til + h2*(h2.'*Uele) + Dh1*(h12til.'*Uele) + (h12til)*(Dh1.'*Uele) + h4til*(h4til.'*Uele) - Dh3til*(h43til.'*Uele) - (h43til)*(Dh3til.'*Uele))...
                                        +G*J*H2*(Dh6til - h2*(Dh4.'*Uele) - Dh4*(h2.'*Uele)));
                            elseif mod(i,2)==0
                                Fge = Fge + 4*h/3*L*...
                                        (E*I*(Dh2*(Dh2.'*Uele) + Dh4*(Dh4.'*Uele))...
                                        +G*A*Ks*(h43til*(h43til.'*Uele) + h12til*(h12til.'*Uele))...
                                        +E*A*H1*(Dh5til + h2*(h2.'*Uele) + Dh1*(h12til.'*Uele) + (h12til)*(Dh1.'*Uele) + h4til*(h4til.'*Uele) - Dh3til*(h43til.'*Uele) - (h43til)*(Dh3til.'*Uele))...
                                        +G*J*H2*(Dh6til - h2*(Dh4.'*Uele) - Dh4*(h2.'*Uele)));
                            elseif mod(i,2)~=0
                                Fge = Fge + 2*h/3*L*...
                                        (E*I*(Dh2*(Dh2.'*Uele) + Dh4*(Dh4.'*Uele))...
                                        +G*A*Ks*(h43til*(h43til.'*Uele) + h12til*(h12til.'*Uele))...
                                        +E*A*H1*(Dh5til + h2*(h2.'*Uele) + Dh1*(h12til.'*Uele) + (h12til)*(Dh1.'*Uele) + h4til*(h4til.'*Uele) - Dh3til*(h43til.'*Uele) - (h43til)*(Dh3til.'*Uele))...
                                        +G*J*H2*(Dh6til - h2*(Dh4.'*Uele) - Dh4*(h2.'*Uele)));
                            end
                        eta = eta + h;
                    end
                    
                    % Transforming elemental force vector for elements
                    % below a bend
                        if n >= non-BendEle
                            Fge = Tbend.'*Fge;
                        else
                        end

                    % Adding elemental force vecotr to global force vector
                        for i = 1:12
                            Fg((n-1)*6+i) = Fg((n-1)*6+i) + Fge(i);
                        end
                n = n+1;
            end
    end
    function [Fw] = FwVector(BCs,Data,U,Stab,Uxo,DB,muw)
    %FwVector Calculates wall contact force vector
    %   Fw = global wall contact force vector
    
        % Extracting intial parameters and initializing data
            noe = Data(1,3);
            nStab = Data(1,4);
            TFO = Data(1,10)*pi/180;
            kw = BCs(1,2);
            Uele = zeros(12,1);
            Fw = zeros(6*(noe+1),1);
            
        % Elemental Loop
        n = 1;
        for w = noe:-1:1
            % Zeroing elemental wall contact force vector
                Fwe = zeros(12,1);
            % Elemental OD
                OD = Data(w+1,1)/12;
            %Extracting elemental displacements
                for i = 1:12
                    Uele(i) = U((n-1)*6+i);
                end
                x1 = Uele(1) + sin(-TFO)*Uxo(n);
                x2 = Uele(7) + sin(-TFO)*Uxo(n+1);
                y1 = Uele(3) + cos(-TFO)*Uxo(n);
                y2 = Uele(9) + cos(-TFO)*Uxo(n+1);
                dx1 = Uele(2);
                dx2 = Uele(8);
                dy1 = Uele(4);
                dy2 = Uele(10);
            
            % Calculating contact forces
                nodenumber2 = w;
                nodenumber1 = w+1;

                % Checking allowable clearance at element node 1 (r1)
                % If node has a stabilizer:
                    rcheck1 = 0;
                    if nStab > 0
                        for i = 1:nStab
                            if Stab(i,1) == nodenumber1
                                r1 = (DB - Stab(i,7)/12)/2;
                                OD1 = DB - r1*2;
                                rcheck1 = 1;
                            else
                            end
                        end
                    else
                    end
                % If node does not have a stabilizer
                    if rcheck1 == 0
                        r1 = (DB-OD)/2;
                        OD1 = OD;
                    else
                    end
                % Checking allowable clearance at element node 1 (r1)
                % If node has a stabilizer:
                    rcheck2 = 0;
                    if nStab > 0
                        for i = 1:nStab
                            if Stab(i,1) == nodenumber2
                                r2 = (DB - Stab(i,7)/12)/2;
                                OD2 = DB - r2*2;
                                rcheck2 = 1;
                            else
                            end
                        end
                    else
                    end
                % If node does not have a stabilizer
                    if rcheck2 == 0
                        r2 = (DB-OD)/2;
                        OD2 = OD;
                    else
                    end
                % Clearance between nodes and wellbore wall
                    f1 = r1 - sqrt(x1^2+y1^2);
                    f2 = r2 - sqrt(x2^2+y2^2);

                % Contact Forces at Node 1
                    if f1 < 0
                        Fn1 = kw*(sqrt(x1^2+y1^2)-r1)/(sqrt(x1^2+y1^2))*[x1;0;y1;0;0;0];
                        Ft1 = muw*kw*(sqrt(x1^2+y1^2)-r1)/(sqrt(x1^2+y1^2))*[-y1;0;x1;0;dx1*y1-(dy1)*x1;0];
                        Fm1 = muw*kw*OD1/2*(sqrt(x1^2+y1^2)-r1)*[0;0;0;-dx1;0;1];
                        Fwe1 = Fn1 + Ft1 + Fm1;
                    else
                        Fwe1 = [0;0;0;0;0;0];
                    end
                % Contact forces at Node 2
                    if f2 < 0
                        Fn2 = kw*(sqrt(x2^2+y2^2)-r2)/(sqrt(x2^2+y2^2))*[x2;0;y2;0;0;0];
                        Ft2 = muw*kw/(sqrt(x2^2+y2^2))*(sqrt(x2^2+y2^2)-r2)*[-y2;0;x2;0;dx2*y2-(dy2)*x2;0];
                        Fm2 = muw*kw*OD2/2*(sqrt(x2^2+y2^2)-r2)*[0;0;0;-dx2;0;1];
                        Fwe2 = Fn2 + Ft2 + Fm2;
                    else
                        Fwe2 = [0;0;0;0;0;0];
                    end
                % Total elemental contact force vector
                    for i=1:6
                        Fwe(i) = Fwe1(i);
                        Fwe(i+6) = Fwe2(i);
                    end
                % Adding nodal contact forces to global contact force
                % vector
                    for i = 1:12
                        Fw((n-1)*6+i) = Fwe(i);
                    end
            n = n + 1;
        end
    end
    function [Tbend] = TransformBend(BendAngle)
        %TransformBend Determines transformation matrix for elements below
        %the motor's bend
        %   Tbend = Elemental Transformation Matrix
            TD1 = [1,0,0,0,0,0,0,0,0,0,0,0;...
                   0,1,0,0,0,0,0,0,0,0,0,0;...
                   0,0,cos(BendAngle),0,sin(BendAngle),0,0,0,0,0,0,0;...
                   0,0,0,cos(BendAngle),0,0,0,0,0,0,0,0;...
                   0,0,-sin(BendAngle),0,cos(BendAngle),0,0,0,0,0,0,0;...
                   0,0,0,0,0,1,0,0,0,0,0,0;...
                   0,0,0,0,0,0,1,0,0,0,0,0;...
                   0,0,0,0,0,0,0,1,0,0,0,0;...
                   0,0,0,0,0,0,0,0,cos(BendAngle),0,sin(BendAngle),0;...
                   0,0,0,0,0,0,0,0,0,cos(BendAngle),0,0;...
                   0,0,0,0,0,0,0,0,-sin(BendAngle),0,cos(BendAngle),0;...
                   0,0,0,0,0,0,0,0,0,0,0,1];
             Tbend = TD1.';
    end

function [ExcDataMat] = Loading(U,Uxo,Data,z,kw,Stab,inc,DB)
    % Loading calculates the internal loads of the BHA/drillstring based on
    % the final calculated deflection
    %   ExcDataMat = Array of loading results
    %       Column 1    = Distance from Bit (ft)
    %       Column 2    = HS def(in)
    %       Column 3    = HS rot/Sag (deg)
    %       Column 4    = RS def(in)
    %       Column 5    = RS rot/Sag (deg)
    %       Column 6    = Axial def(in)
    %       Column 7    = Torsion (deg)
    %       Column 8    = High-Side Contact Force (lbf)
    %       Column 9    = Right-Side Contact Force (lbf)
    %       Column 10   = Total Contact Force (lbf)
    %       Column 11   = HS Bending Moment (ft-lbf)
    %       Column 12   = RS Bending Moment (ft-lbf)
    %       Column 13   = Total Bending Moment (ft-lbf)
    %       Column 14   = HS Shear Load (lbf)
    %       Column 15   = RS Shear Load (lbf)
    %       Column 16   = Total Shear Load (lbf)
    %       Column 17   = Tension (lbf)
    %       Column 18   = Torque (ft-lbf)
    
        % Global variable for determining location of applied steer force
            global indRSSap;
        % Extracing initial parameters
            non         = size(Data,1);
            TFO         = Data(1,10)*pi/180;    % Tool Face Orientation of mud motor/RSS active pad
            BendAngle   = Data(1,11)*pi/180;    % Bend Angle of Mud Motor
            BendEle     = Data(2,11);           % Number of elements between bit and bend
            
        % Initializing Data
            ExcDataMat = zeros(non,18);
            Fwalltot = zeros(non,1);
            
        % Determining if a push-the-bit RSS is present
            ind = find(Stab(:,5) == indRSSap,1);
            if isempty(ind) == true
                RSS = 0;
            else
                RSS = 1;
            end
    
        % Creating position vectors, for plotting
            pos = zeros(non,1);
            s = 0;
            n = 1;
            for w = 2:non
               pos(n,1) = s;
               n = n+1;
               s = s + Data(w,3);
            end
            pos(n,1) = s;
        % Calling functions to determine internal loading and contact
        % forces
            [BeMo,BeMoHs,BeMoRs,ShLo,Tens,Torq,ShLoHS,ShLoRS] = Loads3D(U,Data,z,inc,0);
            [FContact] = FwContact(Data,U,Uxo,Stab,kw,DB,TFO);
            
        % Determining High-Side and Right-Side contact forces
            Fchs = zeros(non,1);
            Fcrs = zeros(non,1);
            w = non;
            for n = 1:non
                Fchs(n) = FContact(w,1)*cos(FContact(w,2));
                Fcrs(n) = -FContact(w,1)*sin(FContact(w,2));
                w = w - 1;
            end
        % Rearranging wall contact force vecotr
            w = non;
            for i = 1:non
                Fwalltot(i) = FContact(w);
                w = w - 1;
            end
            
        % Final results
            ShLoHS(1) = -ShLoHS(1);
            ShLoRS(1) = -ShLoRS(1);
            w = non;
            for i = 1:non
               ExcDataMat (i,1) = pos(i);
               ExcDataMat (i,2) = U(6*(w-1)+3)*12 + cos(TFO)*Uxo(w)*12;
               ExcDataMat (i,4) = -U(6*(w-1)+1)*12 + sin(TFO)*Uxo(w)*12;
               if RSS ~= 1 && w >= non-BendEle
                   ExcDataMat (i,3) = U(6*(w-1)+4)*180/pi + cos(TFO)*BendAngle*180/pi;
                   ExcDataMat (i,5) = -(U(6*(w-1)+2)*180/pi - sin(TFO)*BendAngle*180/pi);
               else
                   ExcDataMat (i,3) = U(6*(w-1)+4)*180/pi;
                   ExcDataMat (i,5) = -U(6*(w-1)+2)*180/pi;
               end
               ExcDataMat (i,6) = U(6*(w-1)+5)*12;
               ExcDataMat (i,7) = U(6*(w-1)+6)*180/pi - U(6)*180/pi;
               ExcDataMat (i,8) = Fchs(i);
               ExcDataMat (i,9) = Fcrs(i);
               ExcDataMat (i,10) = FContact(w);
               ExcDataMat (i,11) = BeMoHs(i);
               ExcDataMat (i,12) = BeMoRs(i);
               ExcDataMat (i,13) = BeMo(i);
               ExcDataMat (i,14) = ShLoHS(i);
               ExcDataMat (i,15) = ShLoRS(i);
               ExcDataMat (i,16) = ShLo(i);
               ExcDataMat (i,17) = Tens(i);
               ExcDataMat (i,18) = -Torq(i);
               w = w - 1;
            end
end

    function [BeMo,BeMoHs,BeMoRs,ShLo,Tens,Torq,ShLohs,ShLors] = Loads3D(U,Data,z,inc,azi)
        % Loads3D Calculates the internal loading in the BHA based on the
        % calculated deflection
        %   BeMo    = Total bending moment at each node [ft-lbf]
        %   BeMoHS  = High-Side bending moment at each node [ft-lbf]
        %   BeMoRS  = Right-Side bending moment at each node [ft-lbf]
        %   ShLo    = Total shear load at each node [lbf]
        %   Tens    = Tension at each node [lbf]
        %   Torq    = Torque at each node [ft-lbf]
        %   ShLohs  = High-Side Shear Load at each node [lbf]
        %   ShLors  = Right-Side shear load at each node [lbf]
        
        % Converting inclination and azimuth to radians
            inc = inc*pi/180;
            azi = azi*pi/180;
        % Extracing initial parameters
            BF = Data(1,1);
            noe = Data(1,3);
            non = noe + 1;
            BendEle = Data(2,11);
            BendAngle = Data(1,11)*pi/180;
        % Initializing Data    
            BeMo = zeros(non,1);
            BeMoHs = zeros(non,1);
            BeMoRs = zeros(non,1);
            BeMoment = zeros(non,1);
            ShearLoad = zeros(non,1);
            Shearhs = zeros(non,1);
            Shearls = zeros(non,1);
            Bendinghs = zeros(non,1);
            Bendingls = zeros(non,1);
            ShLo = zeros(non,1);
            ShLohs = zeros(non,1);
            ShLors = zeros(non,1);
            Tension = zeros(non,1);
            Tens = zeros(non,1);
            Torque = zeros(non,1);
            Torq = zeros(non,1);
            Uele = zeros(12,1);
            
        % Calculate gravitational components
            NormMat = [1,0,0;0,cos(inc),sin(inc);0,-sin(inc),cos(inc)].'*...
                      [cos(pi/2-azi),sin(pi/2-azi),0;-sin(pi/2-azi),cos(pi/2-azi),0;0,0,1].'*...
                      [1,0,0;0,1,0;0,0,1];

            v1 = [NormMat(1,1),NormMat(1,2),NormMat(1,3)];
            v2 = [NormMat(2,1),NormMat(2,2),NormMat(2,3)];
            v3 = [NormMat(3,1),NormMat(3,2),NormMat(3,3)];

            g1 = dot(v1,[0,0,1]);
            g2 = dot(v2,[0,0,1]);
            g3 = dot(v3,[0,0,1]);
            
        % Calculating elemental loading
            w = noe;
            for n = 1:noe
                % Calculate external load vector for element
                    % Zeroing elemental external load vecotr
                        Qele = zeros(12,1);
                    % Extract elemental data
                        E = Data(w+1,5)*144;
                        I = Data(w+1,7)/12^4;
                        OD = Data(w+1,1)/12;
                        ID = Data(w+1,2)/12;
                        nu = Data(w+1,6); 
                        A = pi/4*(OD^2-ID^2);
                        G = E/(2*(1+nu));
                        q = Data(w+1,4)*BF;
                        L = Data(w+1,3);
                        J = 2*I;
                        
                    % Shear correction factor
                        Ks = 6*(1+nu)*(1+(ID/OD)^2)^2/((7+6*nu)*(1+(ID/OD)^2)^2+4*(5+3*nu)*(ID/OD)^2);
                        Lambda = 12*E*I/(Ks*G*A*L^2);
                    % Integrate external load vector for element
                        h = 1/z;
                        eta = 0;
                        for i = 1:(z+1)
                            N1 = 1/(1-Lambda)*(2*eta^3-3*eta^2+Lambda*eta+1-Lambda);
                            N2 = L/(1-Lambda)*(eta^3+eta^2*(1/2*Lambda-2)+eta*(1-1/2*Lambda));
                            N3 = 1/(1-Lambda)*(-2*eta^3+3*eta^2-Lambda*eta);
                            N4 = L/(1-Lambda)*(eta^3-eta^2*(1+1/2*Lambda)+eta*Lambda/2);
                            N9 = 1-eta;
                            N10 = eta;

                            h1 = [N1;N2;0;0;0;0;N3;N4;0;0;0;0];
                            h3 = [0;0;N1;N2;0;0;0;0;N3;N4;0;0];
                            h5 = [0;0;0;0;N9;0;0;0;0;0;N10;0];

                            if i == 1
                                Qele = Qele + h/3*L*...
                                    q*(g1*h1+g2*h3+g3*h5);
                            elseif i == z+1
                                Qele = Qele + h/3*L*...
                                    q*(g1*h1+g2*h3+g3*h5);
                            elseif mod(i,2)==0
                                Qele = Qele + 4*h/3*L*...
                                    q*(g1*h1+g2*h3+g3*h5);
                            elseif mod(i,2)~=0
                                Qele = Qele + 2*h/3*L*...
                                    q*(g1*h1+g2*h3+g3*h5);
                            end
                            eta = eta + h;
                        end

                % Calculate Internal Force Vector
                    % Zeroing elemental internal force vector
                        Fge = zeros(12,1);
                    % Extracing elemental displacement from global
                    % displacement vector
                        for i = 1:12
                            Uele(i) = U((n-1)*6+i);
                        end
                    % Transformation of element below the motor's bend
                        if n >= non-BendEle
                            [Tbend] = TransformBend(BendAngle);
                            Uele = Tbend*Uele;
                        else
                        end

                h = 1/z;
                eta = 0;

                Ks = 6*(1+nu)*(1+(ID/OD)^2)^2/((7+6*nu)*(1+(ID/OD)^2)^2+4*(5+3*nu)*(ID/OD)^2);
                Lambda = 12*E*I/(Ks*G*A*L^2);
                for i = 1:(z+1)
                    % Required shape functions
                        N5 = 6/(L*(1-Lambda))*(eta^2-eta);
                        N6 = 1/(1-Lambda)*(3*eta^2+eta*(Lambda-4)+1-Lambda);
                        N7 = 6/(L*(1-Lambda))*(eta-eta^2); 
                        N8 = 1/(1-Lambda)*(3*eta^2-eta*(Lambda+2));

                        DN1 = 1/(L*(1-Lambda))*(6*eta^2-6*eta+Lambda);
                        DN2 = 1/(1-Lambda)*(3*eta^2+2*eta*(1/2*Lambda-2)+1-1/2*Lambda);
                        DN3 = 1/(L*(1-Lambda))*(-6*eta^2+6*eta-Lambda);
                        DN4 = 1/(1-Lambda)*(3*eta^2-2*eta*(1+1/2*Lambda)+1/2*Lambda);
                        DN5 = 6/(L^2*(1-Lambda))*(2*eta-1);
                        DN6 = 1/(L*(1-Lambda))*(6*eta+Lambda-4);
                        DN7 = 6/(L^2*(1-Lambda))*(1-2*eta); 
                        DN8 = 1/(L*(1-Lambda))*(6*eta-Lambda-2);
                        DN9 = -1/L;
                        DN10 = 1/L;

                        h2 = [N5;N6;0;0;0;0;N7;N8;0;0;0;0];
                        h4 = [0;0;N5;N6;0;0;0;0;N7;N8;0;0];

                        Dh1 = [DN1;DN2;0;0;0;0;DN3;DN4;0;0;0;0];
                        Dh2 = [DN5;DN6;0;0;0;0;DN7;DN8;0;0;0;0];
                        Dh3 = [0;0;DN1;DN2;0;0;0;0;DN3;DN4;0;0];
                        Dh4 = [0;0;DN5;DN6;0;0;0;0;DN7;DN8;0;0];
                        Dh5 = [0;0;0;0;DN9;0;0;0;0;0;DN10;0];
                        Dh6 = [0;0;0;0;0;DN9;0;0;0;0;0;DN10];

                        Dh3til = Dh3;
                        Dh5til = Dh5;
                        Dh6til = Dh6;
                        h4til = h4;
                        h43til = h4 - Dh3;
                        h12til = Dh1 - h2;

                        H1 = (Dh5til.'*Uele) + 1/2*(h2.'*Uele)^2 + 1/2*(h4til.'*Uele)^2 + (Dh1.'*Uele)*(h12til.'*Uele) - (h43til.'*Uele)*(Dh3til.'*Uele);
                        H2 = (Dh6til.'*Uele) - (Dh4.'*Uele)*(h2.'*Uele);
                    % Integration using Simpson's 1/3 rule
                        if i == 1 || i == z+1
                            Fge = Fge + h/3*L*...
                                    (E*I*(Dh2*(Dh2.'*Uele) + Dh4*(Dh4.'*Uele))...
                                    +G*A*Ks*(h43til*(h43til.'*Uele) + h12til*(h12til.'*Uele))...
                                    +E*A*H1*(Dh5til + h2*(h2.'*Uele) + Dh1*(h12til.'*Uele) + (h12til)*(Dh1.'*Uele) + h4til*(h4til.'*Uele) - Dh3til*(h43til.'*Uele) - (h43til)*(Dh3til.'*Uele))...
                                    +G*J*H2*(Dh6til - h2*(Dh4.'*Uele) - Dh4*(h2.'*Uele)));
                        elseif mod(i,2)==0
                            Fge = Fge + 4*h/3*L*...
                                    (E*I*(Dh2*(Dh2.'*Uele) + Dh4*(Dh4.'*Uele))...
                                    +G*A*Ks*(h43til*(h43til.'*Uele) + h12til*(h12til.'*Uele))...
                                    +E*A*H1*(Dh5til + h2*(h2.'*Uele) + Dh1*(h12til.'*Uele) + (h12til)*(Dh1.'*Uele) + h4til*(h4til.'*Uele) - Dh3til*(h43til.'*Uele) - (h43til)*(Dh3til.'*Uele))...
                                    +G*J*H2*(Dh6til - h2*(Dh4.'*Uele) - Dh4*(h2.'*Uele)));
                        elseif mod(i,2)~=0
                            Fge = Fge + 2*h/3*L*...
                                    (E*I*(Dh2*(Dh2.'*Uele) + Dh4*(Dh4.'*Uele))...
                                    +G*A*Ks*(h43til*(h43til.'*Uele) + h12til*(h12til.'*Uele))...
                                    +E*A*H1*(Dh5til + h2*(h2.'*Uele) + Dh1*(h12til.'*Uele) + (h12til)*(Dh1.'*Uele) + h4til*(h4til.'*Uele) - Dh3til*(h43til.'*Uele) - (h43til)*(Dh3til.'*Uele))...
                                    +G*J*H2*(Dh6til - h2*(Dh4.'*Uele) - Dh4*(h2.'*Uele)));
                        end
                        eta = eta + h;
                end

                % Resultant elemental force vector
                    DummyMat = (Fge-Qele);
                % Transformation of elemental forces below the motor's bend
                    if n >= non-BendEle
                        DummyMat = Tbend.'*(DummyMat);
                    else
                    end   
                % Extracting Loading
                    Mb1 = sqrt(DummyMat(2)^2+DummyMat(4)^2);
                    Mb2 = sqrt(DummyMat(8)^2+DummyMat(10)^2);
                    Shear1 = sqrt(DummyMat(1)^2+DummyMat(3)^2);
                    Shear2 = sqrt(DummyMat(7)^2+DummyMat(9)^2);
                    Shearx1 = DummyMat(1);
                    Sheary1 = DummyMat(3);
                    Shearx2 = DummyMat(7);
                    Sheary2 = DummyMat(9);
                    Bendx1 = DummyMat(2);
                    Bendy1 = DummyMat(4);
                    Bendx2 = DummyMat(8);
                    Bendy2 = DummyMat(10);
                    Shearhs(n) = Sheary1;
                    Shearhs(n+1) = Sheary2;
                    Shearls(n) = Shearx1;
                    Shearls(n+1) = Shearx2;
                    Bendingls(n) = Bendx1;
                    Bendingls(n+1) = Bendx2;
                    Bendinghs(n) = Bendy1;
                    Bendinghs(n+1) = Bendy2;
                    BeMoment(n) = Mb1;
                    BeMoment(n+1) = Mb2;
                    ShearLoad(n) = Shear1;
                    ShearLoad(n+1) = Shear2;
                    Tension(n) = -DummyMat(5);
                    Tension(n+1) = DummyMat(11);
                    Torque(n) = DummyMat(6);
                    Torque(n+1) = -DummyMat(12);

                w = w-1;
            end
            
            % Final Results
            w = non;
            for i = 1:non
                BeMo(i) = BeMoment(w);
                ShLo(i) = ShearLoad(w);
                ShLohs(i) = Shearhs(w);
                ShLors(i) = -Shearls(w);
                BeMoHs(i) = Bendinghs(w);
                BeMoRs(i) = -Bendingls(w);
                Tens(i) = Tension(w);
                Torq(i) = Torque(w);
                w = w - 1;
            end
    end
    function [FContact] = FwContact(Data,U,Uxo,Stab,kw,DB,TFO)
        %FwContact Calculates the total resultant contact force
        %   FContact = Vector of normal contact forces along BHA/drillstring

            global indRSSap;    % Global variable for identifying steer force location
        
        
            % Extracting initial parameters
                noe = size(Data,1)-1;
                nStab = size(Stab,1);

            % Initial variables
                Uele = zeros(12,1);
                FContact = zeros(noe+1,2);
            % Elemental loop    
                n = 1;
                for w = noe:-1:1
                    % Extracting Element OD
                        OD = Data(w+1,1)/12;
                    % Extracting Elemental Displacement
                        for i = 1:12
                            Uele(i) = U((n-1)*6+i);
                        end
                        x1 = Uele(1) + sin(-TFO)*Uxo(n);
                        x2 = Uele(7) + sin(-TFO)*Uxo(n+1);
                        y1 = Uele(3) + cos(-TFO)*Uxo(n);
                        y2 = Uele(9) + cos(-TFO)*Uxo(n+1);

                    nodenumber2 = w;
                    nodenumber1 = w+1;
                    % Determine allowable clearance at element node 1 (r1)
                        rcheck1 = 0;
                        if nStab > 0
                            for i = 1:nStab
                                if Stab(i,1) == nodenumber1
                                        r1 = (DB - Stab(i,7)/12)/2;
                                        rcheck1 = 1;
                                else
                                end
                            end
                        else
                        end

                        if rcheck1 == 0
                            r1 = (DB-OD)/2;
                        else
                        end

                    % Determine allowable clearance at element node 2 (r2)
                        rcheck2 = 0;
                        if nStab > 0
                            for i = 1:nStab
                                if Stab(i,1) == nodenumber2
                                        r2 = (DB - Stab(i,7)/12)/2;
                                        rcheck2 = 1;
                                else
                                end
                            end
                        else
                        end

                        if rcheck2 == 0
                            r2 = (DB-OD)/2;
                        else
                        end

                    % Determine clearance at node 1 and node 2 (f1 & f2)
                        f1 = r1 - sqrt(x1^2+y1^2);
                        f2 = r2 - sqrt(x2^2+y2^2);
                    % Calculate normal contact force and contact direction at each node
                        if f1 < 0
                            Fn1 = kw*(sqrt(x1^2+y1^2)-r1);
                            theta1 = atan2(x1,y1)+pi;
                        else
                            Fn1 = 0;
                            theta1 = 0;
                        end

                        if f2 < 0
                            Fn2 = kw*(sqrt(x2^2+y2^2)-r2);
                            theta2 = atan2(x2,y2)+pi;
                        else
                            Fn2 = 0;
                            theta2 = 0;
                        end

                        FContact((n-1)+1,1) = Fn1;
                        FContact((n-1)+1,2) = theta1;
                        FContact((n-1)+2,1) = Fn2;
                        FContact((n-1)+2,2) = theta2;

                    n = n + 1;
                end
            
            % Adding applied steer force for RSS
                for i = 1:nStab
                    if Stab(i,5) == indRSSap
                        nPad = Stab(i,1);
                        Fsteer = Stab(i,2);
                        Fdir = Stab(i,4);
                        FContact(noe+1-nPad+1,1) = Fsteer;
                        FContact(noe+1-nPad+1,2) = Fdir;
                    else
                    end
                end
    end