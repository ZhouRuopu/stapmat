%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses                                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/GetStress.m                                      *
%*                                                                 *
%* - Programmed by:                                                *
%*     Rope Zhou                                                   *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function BeamStress(NG, T)

% Get global data
global cdata;
global sdata;

IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; NUMNP = cdata.NUMNP;
E = sdata.E; AREA = sdata.AREA; nu=sdata.nu;
N1 = sdata.N1; N2 = sdata.N2; N3 = sdata.N3;
STRAIN = sdata.STRAIN; STRESS = sdata.STRESS;

UOrign = sdata.DIS(:, T); IJ = sdata.IJ;

U = zeros(6, NUMNP);
STRAIN = zeros(NUMNP, 6);
STRESS = zeros(NUMNP, 6);
countIJ = zeros(NUMNP,1);

for i = 1:NUMNP
    U(:,i) = UOrign(6*i-5:6*i); 
end

fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT        STRESS11_LN       STRESS12_LN       STRESS13_LN       STRESS12_RN       STRESS12_RN       STRESS12_RN\n' ...
    '       NUMBER\n'], NG);

for N = 1:NUME
    MTYPE = MATP(N);
    I = IJ(N,1); J = IJ(N,2);

    E = E(MTYPE); nu = nu(MTYPE);
   
%   compute the length of truss element
    DX = XYZ(1, N) - XYZ(4, N);
    DY = XYZ(2, N) - XYZ(5, N);
    DZ = XYZ(3, N) - XYZ(6, N);
    D = [DX, DY, DZ];
    XL2 = DX*DX + DY*DY + DZ*DZ;
    XL = sqrt(XL2);

%   TRANSFORMATION MATRIX
    LAMBDA(:,1) = D / XL;  
    LAMBDA(1,2) = N1(N);
    LAMBDA(2,2) = N2(N);
    LAMBDA(3,2) = N3(N);

    LAMBDA(:,3) = CROSS_PRODUCT(LAMBDA(:,1),LAMBDA(:,2));
    LAMBDA(:,3) = LAMBDA(:,3)/sqrt(sum(LAMBDA(:,3).*LAMBDA(:,3)));

    LAMBDA(:,2) = CROSS_PRODUCT(LAMBDA(:,3),LAMBDA(:,1));
    LAMBDA(:,2) = LAMBDA(:,2)/sqrt(sum(LAMBDA(:,2).*LAMBDA(:,2)));

% Compute Strain
    du0dx = [1-1/XL, 1/XL]* [U(1,I), U(1,J)]';
    d2w0dx1 = [6/XL2 + 12/XL2, 4/XL-6/XL, 6/XL2 - 12/XL2, -6/XL+2/XL]* [U(3,I), U(5,I), U(3,J), U(5,J)]';
    d2v0dx1 = [6/XL2 + 12/XL2, 4/XL-6/XL, 6/XL2 - 12/XL2, -6/XL+2/XL]* [U(2,I), U(6,I), U(2,J), U(6,J)]';
    epsx1 = du0dx - d2w0dx1 * (XYZ(3, N)+U(3,I)) - d2v0dx1 * (XYZ(2, N)+U(2,I));

    d2w0dx2 = [6/XL2, 4/XL, 6/XL2, 2/XL]* [U(3,I), U(5,I), U(3,J), U(5,J)]';
    d2v0dx2 = [6/XL2, 4/XL, 6/XL2, 2/XL]* [U(2,I), U(6,I), U(2,J), U(6,J)]';
    epsx2 = du0dx - d2w0dx2 * (XYZ(6, N)+U(3,J)) - d2v0dx2 * (XYZ(5, N)+U(2,J));

    dthetadx = [1-1/XL, 1/XL]* [U(4,I), U(4,J)]';
    gammaxy1 = - dthetadx*(XYZ(3, N)+U(3,I));
    gammaxy2 = - dthetadx*(XYZ(6, N)+U(3,J));
    gammaxz1 = dthetadx*(XYZ(2, N)+U(2,I));
    gammaxz2 = dthetadx*(XYZ(5, N)+U(2,J));

    % 坐标变换
    ST1 = LAMBDA*[epsx1, gammaxy1,  gammaxz1; gammaxy1, 0, 0; gammaxz1, 0, 0];
    ST2 = LAMBDA*[epsx2, gammaxy2,  gammaxz2; gammaxy2, 0, 0; gammaxz2, 0, 0];
    epsx1 = ST1(1,1);
    gammaxy1 = ST1(1,2);
    gammaxz1 = ST1(1,3);
    epsx2 = ST2(1,1);
    gammaxy2 = ST2(1,2);
    gammaxz2 = ST2(1,3);

    % 应变组装
    strain_local_I = [epsx1, 0, 0, gammaxy1, 0, gammaxz1];
    strain_local_J = [epsx2, 0, 0, gammaxy2, 0, gammaxz2];

    STRAIN(I,:) = STRAIN(I,:) + strain_local_I;
    STRAIN(J,:) = STRAIN(J,:) + strain_local_J;
    countIJ(I,:) =  countIJ(I,:)+1;
    countIJ(J,:) =  countIJ(J,:)+1;

    D0 = E*(1-nu)/((1+nu)*(1-2*nu));
    d1 = nu/(1-nu); d2 = (1-2*nu)/(2*(1-nu));
    D = D0.*[1, d1, d1, 0, 0, 0;
             d1, 1,  d1, 0, 0, 0;
             d1, d1, 1, 0, 0, 0;
             0, 0, 0, d2, 0, 0;
             0, 0, 0, 0, d2, 0;
             0, 0, 0, 0, 0, d2;];

    stress_local_I = (D*strain_local_I')';
    stress_local_J = (D*strain_local_J')';

    STRESS(I,:) = STRESS(I,:) + stress_local_I;
    STRESS(J,:) = STRESS(J,:) + stress_local_J;

    
    fprintf(IOUT, ' %10d :  %13.5e    %13.5e    %13.5e          %13.5e    %13.5e    %13.5e \n', N, epsx1, gammaxy1, gammaxz1, epsx2, gammaxy2, gammaxz2);
end

% 应力应变磨平操作
for i = 1:NUMNP
    STRAIN(i,:) = STRAIN(i,:)./countIJ(i,:);
    STRESS(i,:) = STRESS(i,:)./countIJ(i,:);
end

sdata.STRAIN = STRAIN;
sdata.STRESS = STRESS;

end

function V3 = CROSS_PRODUCT(V1, V2)
    V3(1)	= V1(2) * V2(3) - V1(3) * V2(2);
    V3(2)	= V1(3) * V2(1) - V1(1) * V2(3);
    V3(3)	= V1(1) * V2(2) - V1(2) * V2(1);
end