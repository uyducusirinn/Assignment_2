%_Author : MELÝH ALTAY

% STUDENT ID : 21732824             X = 8 Y= 2 Z=4
format long g

%ANKARA CARTESIAN COORDINATES
x_ank = 4121948.47632;
y_ank = 2652187.89675;
z_ank = 4069022.78398;
%ISTANBUL CARTESIAN COORDINATES
x_ist = 4208830.22229;
y_ist = 2334850.33116;
z_ist = 4171267.31457
%MERS CARTESIAN COORDINATES
x_mers = 4239149.40946;
y_mers = 2886967.96208;
z_mers = 3778877.05919;
%tTUBI CARTESIAN COORDINATES
x_tubi =  4211317.42954;
y_tubi =  2377865.32801;
z_tubi =  4144663.66099;

R = 6371000; % unit: meters  


%Converting cartesian coordinates to latitude and longitude for Ankara 

latitude = asind(z_ank/R)
longitude = atan2d(y_ank,x_ank)

%deltax, deltay, deltaz (ist-ank)
delta_x_ist = (x_ist - x_ank);   
delta_y_ist = (y_ist - y_ank);
delta_z_ist = (z_ist - z_ank);

%?x, ?y, ?z (mers-ank)
delta_x_mers = (x_mers - x_ank); 
delta_y_mers = (y_mers - y_ank);   % Geodesy Solved  Problems --- Page 2  (I used this formula)
delta_z_mers = (z_mers - z_ank);

%?x, ?y, ?z (tubi-ank)
delta_x_tubi = (x_tubi - x_ank); 
delta_y_tubi = (y_tubi - y_ank);
delta_z_tubi = (z_tubi - z_ank);

% Let's find local geodetic coordinates

%ISTANBUL LOCAL GEODETIC COORDINATES

x_ist_local = -sind(latitude)*((cosd(longitude)*delta_x_ist + sind(longitude)*delta_y_ist)) + cosd(latitude)*delta_z_ist
y_ist_local = -sind(longitude)*delta_x_ist + cosd(longitude)*delta_y_ist
z_ist_local = cosd(latitude)*((cosd(longitude)*delta_x_ist + sind(longitude)*delta_y_ist)) + sind(latitude)*delta_z_ist

%--------------------------------------------------------------------------------------------------------------------------

%MERS LOCAL GEODETIC COORDINATES
x_mers_local = -sind(latitude)*(cosd(longitude)*delta_x_mers + sind(longitude)*delta_y_mers) + cosd(latitude)*delta_z_mers
y_mers_local = -sind(longitude)*delta_x_mers + cosd(longitude)*delta_y_mers
z_mers_local = cosd(latitude)*(cosd(longitude)*delta_x_mers + sind(longitude)*delta_y_mers) + sind(latitude)*delta_z_mers
%--------------------------------------------------------------------------------------------------------------------------

%TUBI LOCAL GEODETIC COORDINATES
x_tubi_local = -sind(latitude)*(cosd(longitude)*delta_x_tubi + sind(longitude)*delta_y_tubi) + cosd(latitude)*delta_z_tubi
y_tubi_local = -sind(longitude)*delta_x_tubi + cosd(longitude)*delta_y_tubi
z_tubi_local = cosd(latitude)*(cosd(longitude)*delta_x_tubi + sind(longitude)*delta_y_tubi) + sind(latitude)*delta_z_tubi


%The chord length (space distance) between origin and target points must be the
%same in global and local coordinate systems:
control = sqrt((delta_x_ist)^2 + (delta_y_ist)^2 + (delta_z_ist)^2)
control2 = sqrt((x_ist_local)^2 + (y_ist_local)^2 + (z_ist_local)^2)

% control ====== control2  (TRUE)


%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------

%%% SECOND QUESTION

%Let's find geodetic polar coordinates


format long g 
%Angle between target and local geodetic north x* direction (azimuth)

a_ist = atand(y_ist_local / x_ist_local)

%distance to target (chord),

k_ist = sqrt(x_ist_local^2 + y_ist_local^2 + z_ist_local^2)

% vertical angle to the target (z? geodetic zenith),

zen_ist = atand((sqrt(x_ist_local^2 + y_ist_local^2)) / z_ist_local)




a_mers = atand(y_mers_local / x_mers_local)

%distance to target (chord),

k_mers = sqrt((x_mers_local)^2 + (y_mers_local)^2 + (z_mers_local)^2)

% vertical angle to the target (z? geodetic zenith),

zen_mers = atand(((sqrt(x_mers_local^2 + y_mers_local^2)) / z_mers_local))





a_tubi = atand(y_tubi_local / x_tubi_local)

%distance to target (chord),

k_tubi = sqrt((x_tubi_local)^2 + (y_tubi_local)^2 + (z_tubi_local)^2)

% vertical angle to the target (z? geodetic zenith),

zen_tubi = atand(((sqrt(x_tubi_local^2 + y_tubi_local^2)) / z_tubi_local))

%----------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------

%THIRD QUESTION


%If this local system is rotated around the vertical axis (zenith) by ? angle

e2 = 0.006739496775 %second eccentricity
c = 6399593.62586  %polar curvature radius

% ISTANBUL 

format long g 
x_rotate_ist = x_ist_local*cosd(a_ist) + y_ist_local*sind(a_ist)
y_rotate_ist = -(x_ist_local)*sind(a_ist) + y_ist_local*cosd(a_ist)   %rotated coordýnates
z_rotate_ist = z_ist_local





nsqr_ist = (e2 * cosd(latitude)*cosd(latitude))
t = tand(latitude)
N_ist = c / sqrt(1 + e2*cosd(latitude)*cosd(latitude))  % N= c / V
Rn_ist = N_ist/ (1 + (nsqr_ist * cosd(latitude)*cosd(latitude)))
%FROM ANKARA TO ISTANBUL NORMAL SECTION LENGTH
Sn_ist = x_rotate_ist*((1 + 1/6*(x_rotate_ist / Rn_ist)^2) - 3/8*((x_rotate_ist /Rn_ist)^2)*(x_rotate_ist / N_ist)* nsqr_ist * t*cos(a_ist))



%MERS
format long g
x_rotate_mers = x_mers_local*cosd(a_mers) + y_mers_local*sind(a_mers)
y_rotate_mers = -(x_mers_local)*sind(a_mers) + y_mers_local*cosd(a_mers)
z_rotate_mers = z_mers_local

nsqr_mers = (e2 * cosd(latitude)*cosd(latitude))
N_mers = c / sqrt(1 + e2*cosd(latitude)*cosd(latitude))
Rn_mers = N_mers/ (1 + (nsqr_mers * cosd(latitude)*cosd(latitude)))
%FROM ANKARA TO MERS NORMAL SECTION LENGTH
Sn_mers = (x_rotate_mers*((1 + 1/6*(x_rotate_mers / Rn_mers)^2) - 3/8*((x_rotate_ist /Rn_mers)^2)*(x_rotate_mers / N_mers)* nsqr_mers * t*cos(a_mers)))*-1


%TUBI

format long g
x_rotate_tubi = x_tubi_local*cosd(a_tubi) + y_tubi_local*sind(a_tubi)
y_rotate_tubi = -(x_tubi_local)*sind(a_tubi) + y_tubi_local*cosd(a_tubi)
z_rotate_tubi = z_tubi_local

nsqr_tubi = (e2 * cosd(latitude)*cosd(latitude))
N_tubi = c / sqrt(1 + e2*cosd(latitude)*cosd(latitude))
Rn_tubi = N_tubi/ (1 + (nsqr_tubi * cosd(latitude)*cosd(latitude)))

%FROM ANKARA TO TUBI NORMAL SECTION LENGTH
Sn_tubi = x_rotate_tubi*((1 + 1/6*(x_rotate_tubi / Rn_tubi)^2) - 3/8*((x_rotate_tubi /Rn_tubi)^2)*(x_rotate_tubi / N_tubi)* nsqr_tubi * t*cos(a_tubi))

%***********************************************************************************************************************************

%we need ýst-tubý normal section length so ý write new code , define new
%lat long . Tubý is the reference point now.

latitude1 = asind(z_tubi/R)
longitude1 = atan2d(y_tubi,x_tubi)


delta_x_ist_tubi = (x_ist - x_tubi);   
delta_y_ist_tubi = (y_ist - y_tubi);
delta_z_ist_tubi = (z_ist - z_tubi);

x_ist_tubi_local = -sind(latitude1)*((cosd(longitude1)*delta_x_ist_tubi + sind(longitude1)*delta_y_ist_tubi)) + cosd(latitude1)*delta_z_ist_tubi
y_ist_tubi_local = -sind(longitude1)*delta_x_ist_tubi + cosd(longitude1)*delta_y_ist_tubi
z_ist_tubi_local = cosd(latitude1)*((cosd(longitude1)*delta_x_ist_tubi + sind(longitude1)*delta_y_ist_tubi)) + sind(latitude1)*delta_z_ist_tubi



format long g 
%Angle between target and local geodetic north x* direction (azimuth)

a_ist_tubi = atand(y_ist_tubi_local / x_ist_tubi_local)

%distance to target (chord),

k_ist_tubi = sqrt(x_ist_tubi_local^2 + y_ist_tubi_local^2 + z_ist_tubi_local^2)

% vertical angle to the target (z? geodetic zenith),

zen_ist_tubi = atand((sqrt(x_ist_tubi_local^2 + y_ist_tubi_local^2)) / z_ist_tubi_local)



format long g 
x_rotate_ist_tubi = x_ist_tubi_local*cosd(a_ist_tubi) + y_ist_tubi_local*sind(a_ist_tubi)
y_rotate_ist_tubi = -(x_ist_tubi_local)*sind(a_ist_tubi) + y_ist_tubi_local*cosd(a_ist_tubi)
z_rotate_ist_tubi = z_ist_tubi_local


nsqr_ist_tubi = (e2 * cosd(latitude1)*cosd(latitude1))
t1 = tand(latitude1)
N_ist_tubi = c / sqrt(1 + e2*cosd(latitude1)*cosd(latitude1))  % N= c / V
Rn_ist_tubi = N_ist_tubi/ (1 + (nsqr_ist_tubi * cosd(latitude1)*cosd(latitude1)))
Sn_ist_tubi = x_rotate_ist_tubi*((1 + 1/6*(x_rotate_ist_tubi / Rn_ist_tubi)^2) - 3/8*((x_rotate_ist_tubi /Rn_ist_tubi)^2)*(x_rotate_ist_tubi / N_ist_tubi)* nsqr_ist_tubi * t1*cos(a_ist_tubi))







% I define delta values


delta_x_mers_tubi = (x_mers - x_tubi);   
delta_y_mers_tubi = (y_mers - y_tubi);
delta_z_mers_tubi = (z_mers - z_tubi);

x_mers_tubi_local = -sind(latitude1)*((cosd(longitude1)*delta_x_mers_tubi + sind(longitude1)*delta_y_mers_tubi)) + cosd(latitude1)*delta_z_mers_tubi
y_mers_tubi_local = -sind(longitude1)*delta_x_mers_tubi + cosd(longitude1)*delta_y_mers_tubi
z_mers_tubi_local = cosd(latitude1)*((cosd(longitude1)*delta_x_mers_tubi + sind(longitude1)*delta_y_mers_tubi)) + sind(latitude1)*delta_z_mers_tubi

format long g 
%Angle between target and local geodetic north x* direction (azimuth)

a_mers_tubi = atand(y_mers_tubi_local / x_mers_tubi_local)

%distance to target (chord),

k_mers_tubi = sqrt(x_mers_tubi_local^2 + y_mers_tubi_local^2 + z_mers_tubi_local^2)

% vertical angle to the target (z? geodetic zenith),

zen_mers_tubi = atand((sqrt(x_mers_tubi_local^2 + y_mers_tubi_local^2)) / z_mers_tubi_local)



format long g 
x_rotate_mers_tubi = x_mers_tubi_local*cosd(a_mers_tubi) + y_mers_tubi_local*sind(a_mers_tubi)
y_rotate_mers_tubi = -(x_mers_tubi_local)*sind(a_mers_tubi) + y_mers_tubi_local*cosd(a_mers_tubi)
z_rotate_mers_tubi = z_mers_tubi_local


nsqr_mers_tubi = (e2 * cosd(latitude1)*cosd(latitude1))
t1 = tand(latitude1)
N_mers_tubi = c / sqrt(1 + e2*cosd(latitude1)*cosd(latitude1))  % N= c / V
Rn_mers_tubi = N_mers_tubi/ (1 + (nsqr_mers_tubi * cosd(latitude1)*cosd(latitude1)))
Sn_mers_tubi = (x_rotate_mers_tubi*((1 + 1/6*(x_rotate_mers_tubi / Rn_mers_tubi)^2) - 3/8*((x_rotate_mers_tubi /Rn_mers_tubi)^2)*(x_rotate_mers_tubi / N_mers_tubi)* nsqr_mers_tubi * t1*cos(a_mers_tubi)))*-1



U1 = (Sn_ist + Sn_tubi + Sn_ist_tubi) /2  % ISTA-ANKR-TUBI
U2 = (-Sn_mers_tubi + Sn_mers + Sn_tubi) /2 % TUBI-ANKR-MERS


%alfa = 2*atand(sqrt((sin((U1-Sn_ist) /R) * sin((U1-Sn_tubi)/R))) / ((sin(U1/R) * sin((U1-Sn_ist_tubi) / R) )))


% BURADAKI TUM ACILARI HESAP MAKINESINDE HESAPLADIM. BU KISIMDA
% PARANTEZLERI YANLIS KOYDUGUM ICIN BOYLE YAPTIM.
% 8.DERS NOTU SAYFA 9 FORMULLER KULLANDIM.


%ELLIPSOIDAL ANGLES ISTA-ANKR-TUBI
alfa1 = 22.1783635322119
beta1 = 154.147123149056
teta1 = 3.67451885466995

%ELLIPSOIDAL ANGLES TUBI-ANKR-MERS
alfa2 = 28.2123360076082
beta2=  130.667572285129
teta2 = 21.1220782059163





%%% LAST QUESTION


%Compute ellipsoidal (spherical) excesses for both triangles

%Solution this problem :  Difference between sum of the spherical triangle
%angles and 180 degrees is called spherical excess.


Difference_1 = (alfa1 + beta1 + teta1) - 180
Difference_2 = (alfa2 + beta2 + teta2) - 180









