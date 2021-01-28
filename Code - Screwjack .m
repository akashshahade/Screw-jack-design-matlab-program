clc;clear;
fprintf('**** Design of Screwjack ****\n\n');
%----------------------------------------
%----Design of screw----%
W = input('Enter Load to be lifted (KN):');
Sst = input('Enter Screw material Tensile Strength (Mpa):');
Sss = input('Enter Screw material Shear Strength (Mpa):');
Pb = input('Enter Bearing pressure between nut and Screw (Mpa):');
Nf = input('Enter Factor of Safety:');
mu = input('Enter Coefficient of friction between screw and nut:');
mu1 = input('Enter coefficient of collar friction:');
Snt = input('Enter Nut Tensile strength (Mpa):'); %Nut tensile strength
Snc = input('Enter Nut compressive strength (Mpa):'); %Nut compressive srength
Sns = input('Enter Nut shear strength (Mpa):'); %Nut shear strength
max_lift = input('Enter lift height (mm):');

%----Design of screw----%

fprintf('####----------Design of SCREW----------####\n\n');
% Safe values

Sigma_S = Sst/Nf;
Tau_S = Sss/Nf;

%Find Core diameter
p = 10;
dc = sqrt((4*W*1000)/(pi*Sigma_S)) ;
do = round(1.6*dc) ;
do1 = do + 2;
dc1 = do1 - p;
dm = ((do1 + dc1)/2);
disp('Screw Major Diameter (do):');
disp(do1);
disp('Screw Minor Diameter (dc):');
disp(dc1);


q = (p/(pi*dm));
alpha = atan(q);
phi = atan(mu);
alpha_1 = alpha*(180/pi);
phi_1 = phi*(180/pi);
disp('Screw Helix angle (degree):');
disp(alpha_1);
disp('Friction angle (degree):');
disp(phi_1);

% torque required to rotate the screw
e = (alpha + phi);
T1 = (W*1000*tan(e)*(dm/2)) ;
disp('Torque required to rotate the screw T1 (Nmm):');
disp(T1);

%torsional shear stress induced due to torque T1 is
tau = ((16*T1)/(pi*(dc1^3))) ;
disp('Torsional Shear stress induced due to torque T1 (Mpa):');
disp(tau);

%direct compressive stress indced due to axial load
sigmaC = (W*1000/((pi/4)*(dc1^2)));
disp('direct compressive stress indced due to axial load (Mpa)');
disp(sigmaC);
fprintf('\nChecking for permissible limits...\n\n');
if tau < Tau_S && sigmaC < Sigma_S
    fprintf('The induced shear and compressive stresses are less than permissible limits.\nSo Screw design is SAFE.\n\n');
else
    fprintf('Screw Design is not Safe\n\n');
end

%Check for maximum principal stresses

% 01 Maximum shear stress theory
z = 4*(tau^2);
x = (sigmaC^2);
v = ((z+x)^0.5);
tau_max = 0.5*v ;
disp('Tau Max acc to Max shear stress theory (Mpa):');
disp(tau_max);

% 02 Maximum normal Stress theory
b = 0.5*sigmaC ;
sigma_c_max = b + tau_max;
disp('Sigma max acc to MAx normal Stress theory (Mpa):');
disp(sigma_c_max);

fprintf('\n\nChecking for maximum principal stresses...\n\n')
if (tau_max < Tau_S) && (sigma_c_max < Sigma_S)
    fprintf('Screw Design is SAFE.\n\n');
else
    fprintf('Screw Design is not Safe.\n');
end

%----Design of Nut----%

fprintf('####----------Design of NUT----------####\n');

Sigma_t_nut = Snt/Nf;
Sigma_c_nut = Snc/Nf;
Tau_nut_safe = Sns/Nf;
fprintf('\nSafe Tensile stress for Nut(Mpa):');
disp(Sigma_t_nut);
disp('Safe Compressive stress for Nut(Mpa):');
disp(Sigma_c_nut);
disp('Safe Shear stress for Nut(Mpa):');
disp(Tau_nut_safe);

%Find Height of nut, considering bearing pressure on nut.
r = (((do1^2)-(dc1^2))*(pi/4))*(Pb);
n = ceil(W*1000/(r));
fprintf('Number of threads in contact with screw:\n');
disp(n);
H = n*p;
disp('Height of nut(mm):');
disp(H);
t = p/2 ;
disp('Thickness of width of screw(mm):');
disp(t);

disp('Transverse shear stress for screw (Mpa):');
Tau_screw = (W*1000)/(pi*dc1*n*t);
disp(Tau_screw);
disp('Transverse Shear stress for nut(Mpa):');
Tau_nut = (W*1000)/(pi*do1*n*t);
disp(Tau_nut);

if (Tau_screw < Tau_S) && (Tau_nut < Tau_nut_safe)
    fprintf('Design of screw and nut is SAFE.\n');
else
    fprintf('Design of screw and nut is NOT Safe.\n');
end


fprintf('\n####----------Nut collar design----------####\n');
bb = ((4*W*1000)/(pi*Sigma_t_nut));
D1 = ceil(((bb) + (do1^2))^0.5) ;
fprintf('\nCollar inner diameter D1 (mm):');
disp(D1);

vv = ((4*W*1000)/(pi*Sigma_c_nut));
D2 = ceil((vv + (D1^2))^0.5) ;
disp('Collar outer diameter D2 (mm):');
disp(D2);

t1 = ceil((W*1000)/(pi*D1*Tau_nut_safe));
disp('Thickness of Nut Collar t1 (mm):');
disp(t1);

fprintf('\n\n####----------Screw head dimensions----------####\n');
D3 = 1.75*do1 ;
fprintf('\nHead diameter D3 (mm):');
disp(D3);
D4 = D3/4 ;
disp('Cup is fitted with a pin of diamater D4 (mm):');
disp(D4);

fprintf('\n\n####----------Cup Dimensions----------####\n');
Dcup = 1.6*do1 ;
dcup = 0.8*do1;
fprintf('\nOuter diameter of cup Dcup (mm):');
disp(Dcup);
disp('Inner diameter of cup dcup (mm):');
disp(dcup);

fprintf('\n\n####----------Handle Dimensions----------####\n');
fprintf('\nTorque required (T2) to overcome friction at the top of Screw:\n');
Rmean = (((D3/2) + (D4/2))/2);
T2 = (mu1)*(W*1000)*Rmean ;
disp(T2);

fprintf('\nHandle is subjected to a torque of:\n');
T_handle = T1 + T2 ;
disp(T_handle);

fprintf('\nAssuming that a person can apply a force of 300 N, length of handle is gives as (mtr):\n');
L_handle = (T_handle)/300 ;
L_handle_mtr = L_handle / 1000 ;
disp(L_handle_mtr);

disp('Diameter of handle DH (mm):');
Dh = ceil(((32*T_handle)/(pi*30))^(1/3)) ;   %Sigma b of handle element
disp(Dh);

fprintf('\n\n####----------Body Dimensions----------####\n');
fprintf('\nDiameter of body at top D5 (mm):');
D5 = 1.5*D2;
disp(D5);
disp('Thickness of body t3 (mm):');
t3 = 0.25*do1;
disp(t3);
disp('Inside diam of body at bottom D6 (mm):');
D6 = ceil(2.25*D2) ;
disp(D6);
disp('Outside diam of body at bottom D7 (mm):');
D7 = ceil(1.75*D6);
disp(D7);
disp('Thickness of base Tb (mm):');
Tb = 2*t1;
disp(Tb);
disp('Height of body Hb (mm):');
Hb = ceil(max_lift + H + (1.75*H)) ;
disp(Hb);
