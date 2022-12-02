% compute results for paper [1]
%
% images
% relative RMSEs and mean structural SSIMs
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%   [1] M. F. Schiffner, "Rhombic grids reduce the number of voxels in fast pulse-echo ultrasound imaging,"
%       2022 IEEE Int. Ultrasonics Symp. (IUS), Venice, Italy, Oct. 2022, pp. 1–4.
%       DOI: https://doi.org/10.1109/IUS54386.2022.9958278
%       arXiv: https://arxiv.org/abs/2210.04818
%       YouTube: https://www.youtube.com/watch?v=T6dkazW5ZuM
%
% -------------------------------------------------------------------------
% ABOUT:
% -------------------------------------------------------------------------
%   author: Martin F. Schiffner
%   date: 2022-11-13
%   modified: 2022-12-02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load RF data acquired from tissue phantom
load( 'data_RF.mat' );

% bandwidth
f_bounds = [ 2.25, 6.75 ] * 1e6;

% F-number for dynamic receive aperture
F_number = 1;

% field of view (FOV)
FOV_x = [ - 64, 64 ] * element_pitch;
FOV_z = ( 64 / 4 + [ 0, 128 ] ) * element_pitch;

%--------------------------------------------------------------------------
% interpolation
%--------------------------------------------------------------------------
% a) numbers of zeros per axis to pad to avoid overlap
N_pad = [ 511, 511 ] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% 1.) predict spectral passband to determine grid spacings
%--------------------------------------------------------------------------
% bounds on the wavenumber
k_lb = 2 * pi * f_bounds( 1 ) / c_avg;
k_ub = 2 * pi * f_bounds( 2 ) / c_avg;

fprintf( 'k_lb = %.1f rad / m\n', k_lb );
fprintf( 'k_ub = %.1f rad / m\n', k_ub );

% bounds on the receive angle phi [ Sect. II-A1 ]
phi_lb = - atan( 1 / ( 2 * F_number ) );
phi_ub = - phi_lb;

fprintf( 'phi_ub = %.1f°\n', rad2deg( phi_ub ) );

% bounds on the angular spatial frequencies [ Eq. (3) ]
k_hat_x_lb = k_ub * ( sin( theta_incident( 1 ) ) + sin( phi_lb ) );
k_hat_x_ub = k_ub * ( sin( theta_incident( end ) ) + sin( phi_ub ) );
k_hat_z_lb = k_lb * ( min( cos( theta_incident ) ) + cos( phi_lb ) );
k_hat_z_ub = k_ub * ( max( cos( theta_incident ) ) + 1 );

fprintf( 'k_hat_x_lb: %.1f\n', k_hat_x_lb );
fprintf( 'k_hat_x_ub: %.1f\n', k_hat_x_ub );
fprintf( 'k_hat_z_lb: %.1f\n', k_hat_z_lb );
fprintf( 'k_hat_z_ub: %.1f\n', k_hat_z_ub );

%--------------------------------------------------------------------------
% 1.) orthogonal grid (usual spacings)
%--------------------------------------------------------------------------
% a) independent parameters
grid_orthogonal_usual.delta_x = element_pitch / 4;
grid_orthogonal_usual.delta_z = grid_orthogonal_usual.delta_x;

% b) dependent parameters
grid_orthogonal_usual.N_x = floor( diff( FOV_x ) / grid_orthogonal_usual.delta_x );
grid_orthogonal_usual.N_z = floor( diff( FOV_z ) / grid_orthogonal_usual.delta_z );
grid_orthogonal_usual.volume = grid_orthogonal_usual.delta_x * grid_orthogonal_usual.delta_z;
grid_orthogonal_usual.offset_x = FOV_x( 1 ) + ( diff( FOV_x ) - ( grid_orthogonal_usual.N_x - 1 ) * grid_orthogonal_usual.delta_x ) / 2;
grid_orthogonal_usual.offset_z = FOV_z( 1 ) + ( diff( FOV_z ) - ( grid_orthogonal_usual.N_z - 1 ) * grid_orthogonal_usual.delta_z ) / 2;
grid_orthogonal_usual.positions_x = grid_orthogonal_usual.offset_x + (0:( grid_orthogonal_usual.N_x - 1 )) * grid_orthogonal_usual.delta_x;
grid_orthogonal_usual.positions_z = grid_orthogonal_usual.offset_z + (0:( grid_orthogonal_usual.N_z - 1 )) * grid_orthogonal_usual.delta_z;
[ grid_orthogonal_usual.X, grid_orthogonal_usual.Z ] = meshgrid( grid_orthogonal_usual.positions_x, grid_orthogonal_usual.positions_z );

grid_orthogonal_usual.positions = [ grid_orthogonal_usual.X( : ), grid_orthogonal_usual.Z( : ) ];

%--------------------------------------------------------------------------
% 2.) orthogonal grid (optimal spacings)
%--------------------------------------------------------------------------
% a) independent parameters
grid_orthogonal_optimal.delta_x = 2 * pi / ( k_hat_x_ub - k_hat_x_lb );
grid_orthogonal_optimal.delta_z = 2 * pi / ( k_hat_z_ub - k_hat_z_lb );

% b) dependent parameters
grid_orthogonal_optimal.N_x = floor( diff( FOV_x ) / grid_orthogonal_optimal.delta_x );
grid_orthogonal_optimal.N_z = floor( diff( FOV_z ) / grid_orthogonal_optimal.delta_z );
grid_orthogonal_optimal.volume = grid_orthogonal_optimal.delta_x * grid_orthogonal_optimal.delta_z;
grid_orthogonal_optimal.offset_x = FOV_x( 1 ) + ( diff( FOV_x ) - ( grid_orthogonal_optimal.N_x - 1 ) * grid_orthogonal_optimal.delta_x ) / 2;
grid_orthogonal_optimal.offset_z = FOV_z( 1 ) + ( diff( FOV_z ) - ( grid_orthogonal_optimal.N_z - 1 ) * grid_orthogonal_optimal.delta_z ) / 2;
grid_orthogonal_optimal.positions_x = grid_orthogonal_optimal.offset_x + (0:( grid_orthogonal_optimal.N_x - 1 )) * grid_orthogonal_optimal.delta_x;
grid_orthogonal_optimal.positions_z = grid_orthogonal_optimal.offset_z + (0:( grid_orthogonal_optimal.N_z - 1 )) * grid_orthogonal_optimal.delta_z;
[ grid_orthogonal_optimal.X, grid_orthogonal_optimal.Z ] = meshgrid( grid_orthogonal_optimal.positions_x, grid_orthogonal_optimal.positions_z );

grid_orthogonal_optimal.positions = [ grid_orthogonal_optimal.X( : ), grid_orthogonal_optimal.Z( : ) ];

%--------------------------------------------------------------------------
% 3.) proposed 120° rhombic grid
%--------------------------------------------------------------------------
% angles
angle_1_deg = 30;
angle_2_deg = 90;

% lengths
length_1 = k_hat_z_ub - k_hat_z_lb;
length_2 = k_hat_z_ub - k_hat_z_lb;

fprintf( 'length_1: %.10f / m (%.1f / m)\n', length_1, length_1 );
fprintf( 'length_2: %.10f / m (%.1f / m)\n', length_2, length_2 );

% spectral basis
dir_1 = [ cos( deg2rad( angle_1_deg ) ), sin( deg2rad( angle_1_deg ) ) ];
dir_2 = [ cos( deg2rad( angle_2_deg ) ), sin( deg2rad( angle_2_deg ) ) ];

% vectors
vec_1 = length_1 * dir_1;
vec_2 = length_2 * dir_2;

% sampling basis
A = [ vec_1; vec_2 ];
B = 2 * pi * eye( 2 );
U = A \ B;

U_norm = vecnorm( U );
grid_rhombic.volume = abs( det( U ) );

grid_rhombic.N_1 = floor( diff( FOV_x ) / U( 1, 1 ) );
grid_rhombic.N_2 = floor( diff( FOV_z ) / U( 2, 2 ) );

grid_rhombic.offset_1 = FOV_x( 1 ) + ( diff( FOV_x ) - ( grid_rhombic.N_1 - 1 ) * U( 1, 1 ) ) / 2;
grid_rhombic.offset_2 = FOV_z( 1 ) + ( diff( FOV_z ) - ( grid_rhombic.N_2 - 1 ) * U( 2, 2 ) ) / 2;

%
N_samples_rhombic_rows = zeros( grid_rhombic.N_2, 1 );
grid_rhombic.positions = zeros( 1, 2 );
index = 0;
for l_2 = 0:( grid_rhombic.N_2 - 1 )

    l_1_lb = ceil( ( FOV_x( 1 ) - grid_rhombic.offset_1 - l_2 * U( 1, 2 ) ) / U( 1, 1 ) );
    l_1_ub = floor( ( FOV_x( 2 ) - grid_rhombic.offset_1 - l_2 * U( 1, 2 ) ) / U( 1, 1 ) );

    for l_1 = l_1_lb:l_1_ub
        index = index + 1;
        grid_rhombic.positions( index, : ) = [ grid_rhombic.offset_1, grid_rhombic.offset_2 ] + l_1 * U( :, 1 ).' + l_2 * U( :, 2 ).';
    end
    N_samples_rhombic_rows( l_2 + 1 ) = l_1_ub - l_1_lb + 1;
end

%--------------------------------------------------------------------------
% 4.) combine grids in cell array
%--------------------------------------------------------------------------
grids = { grid_orthogonal_usual, grid_orthogonal_optimal, grid_rhombic };
str_grids = { 'Orthogonal grid (usual spacings)', 'Orthogonal grid (optimal spacings)', '120° rhombic grid' };

%--------------------------------------------------------------------------
% 5.) display grids
%--------------------------------------------------------------------------
% create new figure
figure( 1 );

% iterate grids
for index_grid = 1:numel( grids )

    % create new subplot
    subplot( 1, numel( grids ), index_grid );
    scatter( grids{ index_grid }.positions( :, 1 ), grids{ index_grid }.positions( :, 2 ) );

    % draw field of view
    line( FOV_x, FOV_z( 1 ) * ones( 1, 2 ), 'Color', 'r' );
    line( FOV_x, FOV_z( 2 ) * ones( 1, 2 ), 'Color', 'r' );
    line( FOV_x( 1 ) * ones( 1, 2 ), FOV_z, 'Color', 'r');
    line( FOV_x( 2 ) * ones( 1, 2 ), FOV_z, 'Color', 'r' );

    % set limits
    xlim( [ FOV_x( 1 ) - 5e-3, FOV_x( 2 ) + 5e-3 ] );
    ylim( [ FOV_z( 1 ) - 5e-3, FOV_z( 2 ) + 5e-3 ] );
    title( str_grids{ index_grid } );

end % for index_grid = 1:numel( grids )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beamform RF data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
index_t0 = 16;

% specify cell array for B-mode images
images_single = cell( numel( grids ), numel( theta_incident ) );
times_elapsed = zeros( numel( grids ), numel( theta_incident ) );

images_cpwc = cell( numel( grids ), 1 );

% iterate grids
for index_grid = 1:numel( grids )

    % iterate steering angles
    for index_theta = 1:numel( theta_incident )

        % call delay-and-sum (DAS) algorithm
        [ images_single{ index_grid, index_theta }, times_elapsed( index_grid, index_theta ) ] = das_pw( grids{ index_grid }.positions, data_RF( :, :, index_theta ), f_s, theta_incident( index_theta ), element_width, element_pitch, c_avg, f_bounds, index_t0, [], F_number );

    end % for index_theta = 1:numel( theta_incident )

    % 2.) create (compound) images
    images_cpwc{ index_grid } = sum( cat( 2, images_single{ index_grid, : } ), 2 );

end % for index_grid = 1:numel( grids )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate all B-mode images to usual orthogonal grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% 0.) frequency axes
%--------------------------------------------------------------------------
N_dsft_x = grid_orthogonal_usual.N_x + N_pad( 1 );
N_dsft_z = grid_orthogonal_usual.N_z + N_pad( 2 );

% shift
index_shift_x = ceil( N_dsft_x / 2);

% axes of discrete angular spatial frequencies
axis_k_x = 2 * pi * (( index_shift_x - N_dsft_x ):( index_shift_x - 1 ) ) / ( N_dsft_x * grid_orthogonal_usual.delta_x );
axis_k_z = 2 * pi * (0:( N_dsft_z - 1 ) ) / ( N_dsft_z * grid_orthogonal_usual.delta_z );

% grid of discrete angular spatial frequencies
[ K_x, K_z ] = meshgrid( axis_k_x, axis_k_z );

%--------------------------------------------------------------------------
% 1.) interpolation function
%--------------------------------------------------------------------------
% specify cell array
filter_dft_single = cell( numel( theta_incident ), 1 );

% iterate steering angles
for index_theta = 1:numel( theta_incident )

    % compute wavenumbers
    K = ( K_x.^2 + K_z.^2 ) ./ ( 2 * ( sin( theta_incident( index_theta ) ) * K_x + cos( theta_incident( index_theta ) ) * K_z ) );

    % restrict to relevant bandwidth
    indicator_k = ( K >= k_lb ) & ( K <= k_ub );

    sin_phi = 10 * ones( size( K ) );
    cos_phi = 10 * ones( size( K ) );
    sin_phi( indicator_k ) = K_x( indicator_k ) ./ K( indicator_k ) - sin( theta_incident( index_theta ) );
    cos_phi( indicator_k ) = K_z( indicator_k ) ./ K( indicator_k ) - cos( theta_incident( index_theta ) );

    filter_dft_single{ index_theta } = indicator_k & ( sin_phi >= sin( phi_lb ) ) & ( sin_phi <= sin( phi_ub ) ) & cos_phi >= cos( phi_lb );    

end % for index_theta = 1:numel( theta_incident )

filter_dft_compound = any( cat( 3, filter_dft_single{ : } ), 3 );
filter_compound = ifft2( filter_dft_compound );
filter_compound_dB = 20 * log10( abs( fftshift( filter_compound ) ) / max( abs( filter_compound( : ) ) ) );

figure( 2 );
subplot( 1, 2, 1 );
imagesc( filter_dft_compound );
title( 'Passband filter (k-space)' );
subplot( 1, 2, 2 );
imagesc( filter_compound_dB, [ -70, 0 ] );
title( 'Interpolation kernel (x-space)' );
colormap gray;

%--------------------------------------------------------------------------
% 2.) usual orthogonal grid
%--------------------------------------------------------------------------
% a) samples
samples_usual = grid_orthogonal_usual.volume * reshape( images_cpwc{ 1 }, [ grid_orthogonal_usual.N_z, grid_orthogonal_usual.N_x ] );
samples_usual_dB = 20 * log10( abs( samples_usual ) / max( abs( samples_usual( : ) ) ) );

% b) spectrum
samples_usual_dft = fftshift( fft2( samples_usual ), 2 );
samples_usual_dft_dB = illustration.dB( samples_usual_dft, 20 );

%--------------------------------------------------------------------------
% 3.) optimal orthogonal grid
%--------------------------------------------------------------------------
% a) samples (coarse)
samples_optimal = reshape( images_cpwc{ 2 }, [ grid_orthogonal_optimal.N_z, grid_orthogonal_optimal.N_x ] );
samples_optimal_dB = 20 * log10( abs( samples_optimal ) / max( abs( samples_optimal(:) ) ) );

% b) spectrum
samples_optimal_dft = grid_orthogonal_optimal.volume * dft_positions( images_cpwc{ 2 }, grids{ 2 }.positions, axis_k_x, axis_k_z );

% apply phase shift
samples_optimal_dft = samples_optimal_dft .* exp( 1j * axis_k_x * grid_orthogonal_usual.offset_x ) .* exp( 1j * axis_k_z.' * grid_orthogonal_usual.offset_z );
samples_optimal_dft_dB = illustration.dB( samples_optimal_dft, 20 );

% c) samples (fine)
samples_optimal_interp = ifft2( ifftshift( samples_optimal_dft .* filter_dft_compound, 2 ) );
samples_optimal_interp = samples_optimal_interp( 1:grid_orthogonal_usual.N_z, 1:grid_orthogonal_usual.N_x );
samples_optimal_interp_dB = illustration.dB( samples_optimal_interp, 20 );

%--------------------------------------------------------------------------
% 4.) rhombic grid
%--------------------------------------------------------------------------
% a) spectrum
samples_rhombic_dft = grid_rhombic.volume * dft_positions( images_cpwc{ 3 }, grid_rhombic.positions, axis_k_x, axis_k_z );

% apply phase shift
samples_rhombic_dft = samples_rhombic_dft .* exp( 1j * axis_k_x * grid_orthogonal_usual.offset_x ) .* exp( 1j * axis_k_z.' * grid_orthogonal_usual.offset_z );

% logarithmic compression
samples_rhombic_dft_dB = illustration.dB( samples_rhombic_dft, 20 );

% b) samples (fine)
samples_rhombic_interp = ifft2( ifftshift( samples_rhombic_dft .* filter_dft_compound, 2 ) );
samples_rhombic_interp = samples_rhombic_interp( 1:grid_orthogonal_usual.N_z, 1:grid_orthogonal_usual.N_x );
samples_rhombic_interp_dB = illustration.dB( samples_rhombic_interp, 20 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute errors and show results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% 1.) metrics
%--------------------------------------------------------------------------
error_orthogonal_optimal = samples_usual - samples_optimal_interp;
error_orthogonal_optimal_dB = 20 * log10( abs( error_orthogonal_optimal ) / max( abs( samples_usual( : ) ) ) );

error_rhombic = samples_usual - samples_rhombic_interp;
error_rhombic_dB = 20 * log10( abs( error_rhombic ) / max( abs( samples_usual( : ) ) ) );

rel_RMSE_orthogonal_optimal = norm( error_orthogonal_optimal( : ) ) / norm( samples_usual( : ) );
rel_RMSE_rhombic = norm( error_rhombic( : ) ) / norm( samples_usual( : ) );

fprintf( 'Rel. RMSE (optimal vs. usual): %.10f %% (%.1f %%)\n', rel_RMSE_orthogonal_optimal * 1e2, rel_RMSE_orthogonal_optimal * 1e2 );
fprintf( 'Rel. RMSE (rhombic vs. usual): %.10f %% (%.1f %%)\n', rel_RMSE_rhombic * 1e2, rel_RMSE_rhombic * 1e2 );

% SSIM_orthogonal_usual = ssim( samples_orthogonal_fine_dB, samples_orthogonal_fine_interp_dB );
SSIM_orthogonal_optimal = ssim( samples_optimal_interp_dB, samples_usual_dB );
SSIM_rhombic = ssim( samples_rhombic_interp_dB, samples_usual_dB );

fprintf( 'Mean SSIM Index (optimal vs. usual): %.10f %% (%.1f %%)\n', SSIM_orthogonal_optimal * 1e2, SSIM_orthogonal_optimal * 1e2 );
fprintf( 'Mean SSIM Index (rhombic vs. usual): %.10f %% (%.1f %%)\n', SSIM_rhombic * 1e2, SSIM_rhombic * 1e2 );

c_limits = [-60,0];

%--------------------------------------------------------------------------
% 2.) show results
%--------------------------------------------------------------------------
figure( 1 );
subplot( 4, 3, 2 );
imagesc( grid_orthogonal_optimal.positions_x, grid_orthogonal_optimal.positions_z, samples_optimal_dB, c_limits );
subplot( 4, 3, 4 );
imagesc( axis_k_x, axis_k_z, samples_usual_dft_dB, c_limits );
subplot( 4, 3, 5 );
imagesc( axis_k_x, axis_k_z, samples_optimal_dft_dB, c_limits );
subplot( 4, 3, 6 );
imagesc( axis_k_x, axis_k_z, samples_rhombic_dft_dB, c_limits );
subplot( 4, 3, 7 );
imagesc( grid_orthogonal_usual.positions_x, grid_orthogonal_usual.positions_z, samples_usual_dB, c_limits );
subplot( 4, 3, 8 );
imagesc( grid_orthogonal_usual.positions_x, grid_orthogonal_usual.positions_z, samples_optimal_interp_dB, c_limits );
title( sprintf( 'mean SSIM index: %.2f %%', SSIM_orthogonal_optimal * 1e2 ) );
subplot( 4, 3, 9 );
imagesc( grid_orthogonal_usual.positions_x, grid_orthogonal_usual.positions_z, samples_rhombic_interp_dB, c_limits );
title( sprintf( 'mean SSIM index: %.2f %%', SSIM_rhombic * 1e2 ) );
% subplot( 4, 3, 10 );
% imagesc( positions_x, positions_z, error_orthogonal_fine_dB, c_limits );
subplot( 4, 3, 11 );
imagesc( grid_orthogonal_usual.positions_x, grid_orthogonal_usual.positions_z, error_orthogonal_optimal_dB, c_limits );
title( sprintf( 'rel. RMSE: %.2f %%', rel_RMSE_orthogonal_optimal * 1e2 ) );
subplot( 4, 3, 12 );
imagesc( grid_orthogonal_usual.positions_x, grid_orthogonal_usual.positions_z, error_rhombic_dB, c_limits );
title( sprintf( 'rel. RMSE: %.2f %%', rel_RMSE_rhombic * 1e2 ) );
colormap gray;
