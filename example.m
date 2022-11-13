% plot results for paper
% experimental validation
%
% images
% widths and FEHMs of wires
%
% influence of the F-number
%
% author: Martin F. Schiffner
% date: 2022-11-13
% modified: 2022-11-13

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
element_width = 279.8e-6;
element_pitch = 304.8e-6;

f_bounds = [ 2.25, 6.75 ] * 1e6;

c_avg = 1538;

f_s = 20e6;

data_RF = randn( 3000, 128 );

theta_incident = deg2rad( [ 80, 90, 110 ] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% 1.) usual orthogonal grid
%--------------------------------------------------------------------------
% a) independent parameters
grid_orthogonal_usual.delta_x = element_pitch / 4;
grid_orthogonal_usual.delta_z = grid_orthogonal_usual.delta_x;
grid_orthogonal_usual.N_x = 512;
grid_orthogonal_usual.N_z = 512;

% b) dependent parameters
grid_orthogonal_usual.M_x = ( grid_orthogonal_usual.N_x - 1 ) / 2;
grid_orthogonal_usual.positions_x = (-grid_orthogonal_usual.M_x:grid_orthogonal_usual.M_x) * grid_orthogonal_usual.delta_x;
grid_orthogonal_usual.positions_z = (1:grid_orthogonal_usual.N_z) * grid_orthogonal_usual.delta_z;
[ grid_orthogonal_usual.X, grid_orthogonal_usual.Z ] = meshgrid( grid_orthogonal_usual.positions_x, grid_orthogonal_usual.positions_z );

grid_orthogonal_usual.positions = [ grid_orthogonal_usual.X( : ), grid_orthogonal_usual.Z( : ) ];

%--------------------------------------------------------------------------
% 2.) optimal orthogonal grid
%--------------------------------------------------------------------------
% a) independent parameters
grid_orthogonal_optimal.delta_x = element_pitch / 4;
grid_orthogonal_optimal.delta_z = grid_orthogonal_usual.delta_x;
grid_orthogonal_optimal.N_x = 512;
grid_orthogonal_optimal.N_z = 512;

% b) dependent parameters
grid_orthogonal_optimal.M_x = ( grid_orthogonal_optimal.N_x - 1 ) / 2;
grid_orthogonal_optimal.positions_x = (-grid_orthogonal_optimal.M_x:grid_orthogonal_optimal.M_x) * grid_orthogonal_optimal.delta_x;
grid_orthogonal_optimal.positions_z = (1:grid_orthogonal_optimal.N_z) * grid_orthogonal_optimal.delta_z;
[ grid_orthogonal_optimal.X, grid_orthogonal_optimal.Z ] = meshgrid( grid_orthogonal_optimal.positions_x, grid_orthogonal_optimal.positions_z );

grid_orthogonal_optimal.positions = [ grid_orthogonal_optimal.X( : ), grid_orthogonal_optimal.Z( : ) ];

%--------------------------------------------------------------------------
% 3.) proposed rhombic grid
%--------------------------------------------------------------------------
grid_rhombic.positions = [];

%--------------------------------------------------------------------------
% 4.) combine grids in cell array
%--------------------------------------------------------------------------
grids = { grid_orthogonal_usual, grid_orthogonal_optimal, grid_rhombic };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beamform RF data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
index_t0 = 16;
F_number = 1;

% specify cell array for B-mode images
images = cell( size( grids ) );

% iterate grids
for index_grid = 1:numel( grids )

    % call delay-and-sum (DAS) algorithm
    images{ index_grid } = das_pw( grids{ index_grid }.positions, data_RF, f_s, theta_incident, element_width, element_pitch, c_avg, f_bounds, index_t0, window, F_number );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate results to orthogonal grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

