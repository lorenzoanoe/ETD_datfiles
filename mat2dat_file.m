sys = 8 ;           % Earth-Moon system
% sys = 1 ;           % Sun-Jupiter system

grid_scale = 400 ;  % lower definition, faster generation
% % % grid_scale = 800 ;  % higher definition, slower generation --> only few datasets are available for the Earth-Moon system and Gamma = 0.80 : 0.04 : 1.36

% Systems list
% 1 - Sun-Jupiter           mu = 9.54e-04
% 2 - Jupiter-Ganymede      mu = 7.80e-05   ABC in collision
% 3 - Jupiter-Europa        mu = 2.53e-05   ABC in collision
% 4 - Sun-Neptune           mu = 5.15e-05
% 5 - Sun-Mars              mu = 3.22e-07
% 6 - Sun-Earth             mu = 3.00e-06
% 7 - Sun-Psyche            mu = 1.15e-11
% 8 - Earth-Moon            mu = 1.22e-02
% 9 - Sun-Pluto             mu = 6.55e-09
% 10 - Pluto-Charon         mu = 1.04e-01   ABC in collision
% e.g. with sys = 8 the Earth-Moon system is chosen from the list

G = 6.674e-20 ;             % [km^3 / (kg*s^2)] universal gravitational constant
m_1vec(1) = 1.9885e+30 ;        m_2vec(1) = 1.89819e27 ;    R_refvec(1) = 778.4e6 ;       R_P_kmvec(1) = 71492 ;
m_1vec(2) = 1.89819E+27 ;       m_2vec(2) = 1.48E+23 ;      R_refvec(2) = 1070400 ;       R_P_kmvec(2) = 5262 ;
m_1vec(3) = 1.89819E+27 ;       m_2vec(3) = 4.80E+22 ;      R_refvec(3) = 670900 ;        R_P_kmvec(3) = 1561 ;
m_1vec(4) = 1.9885e+30 ;        m_2vec(4) = 1.0243E+26 ;    R_refvec(4) = 4498252900 ;    R_P_kmvec(4) = 49528 ;
m_1vec(5) = 1.9885e+30 ;        m_2vec(5) = 6.417E+23 ;     R_refvec(5) = 227900000 ;     R_P_kmvec(5) = 3402.4 ;
m_1vec(6) = 1.9885e+30 ;        m_2vec(6) = 398600/G ;      R_refvec(6) = 149.6e6 ;       R_P_kmvec(6) = 6378.0 ;
m_1vec(7) = 1.9885e+30 ;        m_2vec(7) = 2.287e19 ;      R_refvec(7) = 4.36921e11 ;    R_P_kmvec(7) = 278 ;
m_1vec(8) = 5.9724E+24 ;        m_2vec(8) = 7.348E+22 ;     R_refvec(8) = 384748 ;        R_P_kmvec(8) = 1737.4 ;
m_1vec(9) = 1.9885e+30 ;        m_2vec(9) = 1.303E+22 ;     R_refvec(9) = 5.87e9 ;        R_P_kmvec(9) = 2376.6 ;
m_1vec(10) = 1.303E+22 ;        m_2vec(10) = 1.52E+21 ;     R_refvec(10) = 19591.4 ;      R_P_kmvec(10) = 606 ;

mu_1 = m_1vec(sys) * G ;             % [km^3 / s^2] gravitational constant of m1
mu_2 = m_2vec(sys) * G ;             % [km^3 / s^2] gravitational constant of m2
R_ref = R_refvec(sys) ;
R_P = R_P_kmvec(sys) / R_ref ;       % [km] physical radius of m2

mu = mu_2 / (mu_1+mu_2) ;               % mass ratio of primaries
T = 2*pi*sqrt(R_ref^3/(mu_1+mu_2)) ;    % [s] orbital period
n = 1 / sqrt(R_ref^3/(mu_1+mu_2)) ;     % [rad/s] mean motion
r_Hill = (mu/3)^(1/3) ;                 % dimensionless radius of Hill's sphere

LP = lagrangePoints(mu) ; % position of the Lagrangian points in the baricentric synodic frame
x_CJL = LP(:,1) ;
y_CJL = LP(:,2) ;
CJ_L = (x_CJL).^2 + y_CJL.^2 + 2*(1-mu)./sqrt((x_CJL+mu).^2+y_CJL.^2) + 2*mu./sqrt((x_CJL-(1-mu)).^2+y_CJL.^2) ; % Jacobi constant on the Lagrangian points (velocity=0)

% Position of the Lagrangian points in the synodic frame centered in M1
L1 = LP(1) + mu ;
L2 = LP(2) + mu ;
L3 = LP(3) + mu ;
x_L4 = cosd(60) ;   % x_L4 = LP(4,1) + mu ;
y_L4 = sind(60) ;   % y_L4 = LP(4,2) ;
x_L5 = cosd(-60) ;
y_L5 = sind(-60) ;

% Creation of the grid
A = 3.5*r_Hill ;    % width of the grid generation
B = 4.5*r_Hill ;    % height of the grid generation
step_NA = round(r_Hill/grid_scale, 1, 'significant') ;
vett_x20 = [ 0 : -step_NA : -A, step_NA : step_NA : A*1.2 ] ;
vett_x20 = sort(vett_x20) ;
vett_y20 = [ 0 : -step_NA : -B, step_NA : step_NA : B ] ;
vett_y20 = sort(vett_y20) ;
l_x = length(vett_x20) ;
l_y = length(vett_y20) ;

% Upload Matlab structure for the capture set
addpath(append(pwd, "\datasets"));

output_dir = fullfile(pwd, "datasets_dat") ;
if ~exist(output_dir, 'dir')
    mkdir(output_dir) ;
end

for GAMMA = 0:0.02:1.36      % three-body energy parameter
    
    if grid_scale ~= 400
        mat_filename = append("strsys", num2str(sys), "Gamma", num2str(GAMMA*100), "V16_step", num2str(grid_scale), ".mat") ;
    else
        mat_filename = append("strsys", num2str(sys), "Gamma", num2str(GAMMA*100), "V16.mat") ;
    end

    load( fullfile(pwd, "datasets", mat_filename) ) ;

    [~, base_name, ~] = fileparts(mat_filename) ;
    dat_filename = fullfile(output_dir, append(base_name, ".dat")) ;
    saveStructAsDat(str, dat_filename) ;
    
end

function saveStructAsDat(str, out_filename)
    fields = fieldnames(str) ;
    fields(5) = [] ; % keep this to delete "Contour" polyshape variable
    fields = [fields(2:end); fields(1)];
    fid = fopen(out_filename, 'w') ;
    if fid == -1
        error("Unable to open %s for writing.", out_filename) ;
    end

    for i = 1:numel(fields)
        field_value = str.(fields{i}) ;
        fprintf(fid, '%s\n', formatStructField(field_value)) ;
    end

    fclose(fid) ;
end

function text = formatStructField(value)
    if isempty(value)
        text = '' ;
        return
    end

    if ischar(value) || isstring(value)
        text = strjoin(string(value), '\t') ;
        return
    end

    if iscell(value)
        text = joinCellValues(value) ;
        return
    end

    if isnumeric(value) || islogical(value)
        text = joinNumericValues(value) ;
        return
    end

    if isstruct(value)
        text = '' ;
        return
    end

    text = strtrim(evalc('disp(value)')) ;
end

function s = joinNumericValues(v)
    if isvector(v)
        s = sprintf('%g\t', v) ;
        s = s(1:end-1) ;
    else
        s = [] ;
        for row = 1:size(v, 1)
            row_str = sprintf('%g\t', v(row, :)) ;
            row_str = row_str(1:end-1) ;
            s = [s, row_str, '\n'] ;
        end
        s = s(1:end-1) ;
    end
end

function s = joinCellValues(v)
    s = [] ;
    for k = 1:numel(v)
        item = v{k} ;
        if isnumeric(item) || islogical(item)
            item_txt = joinNumericValues(item) ;
        elseif ischar(item) || isstring(item)
            item_txt = strjoin(string(item), '\t') ;
        else
            item_txt = strtrim(evalc('disp(item)')) ;
        end
        s = [s, item_txt, '\t'] ;
    end
    if ~isempty(s)
        s = s(1:end-1) ;
    end
end