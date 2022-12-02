function [ image, time_elapsed ] = das_pw( positions, data_RF, f_s, theta_incident, element_width, element_pitch, c_avg, f_bounds, index_t0, window, F_number )
% DAS_PW Delay-And-Sum (DAS) Beamforming [ Fourier domain, steered plane wave ]
%
% Computes a B-mode image using
% the delay-and-sum (DAS) algorithm in
% the Fourier domain.
%
% The Fourier domain permits
%  1.) the independent receive focusing of different frequencies,
%  2.) exact corrections of the round-trip times-of-flight independent of the sampling rate, and
%  3.) the usage of frequency-dependent apodization weights.
%
% A uniform linear transducer array is assumed.
%
% This function is a simplified version of
% the function das_pw in [1] that supports
% arbitrary grids [2].
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
%   01.) positions:           voxel positions (m) [ 2d array; 1st dimension: voxel; 2nd dimension: position in 2d space ]
%   02.) data_RF:             RF data (a.u.)      [ 2d array; 1st dimension: time, 2nd dimension: array elements ]
%   03.) f_s:                 sampling rate of the RF data (Hz)
%   04.) theta_incident:      steering angle of the incident plane wave (rad)
%   05.) element_width:       element width of the linear array (m)
%   06.) element_pitch:       element pitch of the linear array (m)
%   07.) c_avg:               speed of sound (m/s)
%
% OPTIONAL
%   08.) f_bounds:            frequency bounds (Hz) [ row vector; 1st component: lower bound, 2nd component: upper bound ]
%   09.) index_t0:            time index of the sample extracted from the focused RF signal (1)
%   10.) window:              window function for receive apodization ( object of class windows.window )
%   11.) F_number:            receive F-number
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%   01.) image:               complex-valued B-mode image
%   02.) F_number_values:     value of the F-number for each frequency
%   03.) signal:              focused RF signal for last image voxel
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%   [1] M. F. Schiffner and G. Schmitz, "Frequency-dependent F-number increases the contrast and the spatial resolution in fast pulse-echo ultrasound imaging,"
%       2021 IEEE Int. Ultrasonics Symp. (IUS), Xi’an, China, Sep. 2021, pp. 1–4.
%       DOI: https://doi.org/10.1109/IUS52206.2021.9593488
%       arXiv: https://arxiv.org/abs/2111.04593
%       YouTube: https://www.youtube.com/watch?v=T6BoYRvQ6rg
%
%   [2] M. F. Schiffner, "Rhombic grids reduce the number of voxels in fast pulse-echo ultrasound imaging,"
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
%   modified: 2022-11-20

% print status
time_start = tic;
str_date_time = sprintf( '%04d-%02d-%02d: %02d:%02d:%02d', fix( clock ) );
fprintf( '\t %s: Delay-and-Sum (DAS) Beamforming [ Fourier domain, steered plane wave ]... ', str_date_time );

%--------------------------------------------------------------------------
% 1.) check arguments
%--------------------------------------------------------------------------
% ensure at least 7 and at most 11 arguments
narginchk( 7, 11 );

% ensure real-valued matrix w/ two columns for positions
if ~( ismatrix( positions ) && ( size( positions, 2 ) == 2 ) && isreal( data_RF ) )
    errorStruct.message = 'positions must be a matrix with two columns!';
    errorStruct.identifier = 'das_pw:NoTwoColumnMatrix';
    error( errorStruct );
end

% ensure real-valued matrix for data_RF
if ~( ismatrix( data_RF ) && isreal( data_RF ) )
    errorStruct.message = 'data_RF must be a matrix!';
    errorStruct.identifier = 'das_pw:NoTwoColumnMatrix';
    error( errorStruct );
end

% ensure positive real-valued scalar f_s
if ~( isscalar( f_s ) && isreal( f_s ) && f_s > 0 )
    errorStruct.message = 'f_s must be positive and scalar!';
    errorStruct.identifier = 'das_pw:NoTwoColumnMatrix';
    error( errorStruct );
end

% ensure positive real-valued scalar theta_incident
if ~( isscalar( theta_incident ) && isreal( theta_incident ) && theta_incident > - pi / 2 && theta_incident < pi / 2 )
    errorStruct.message = 'theta_incident must be scalar and range from - pi / 2 to pi / 2!';
    errorStruct.identifier = 'das_pw:NoTwoColumnMatrix';
    error( errorStruct );
end

% ensure existence of nonempty f_bounds
if nargin < 8 || isempty( f_bounds )
    f_bounds = [ 0, f_s / 2];
end

% ensure existence of nonempty index_t0
if nargin < 9 || isempty( index_t0 )
	index_t0 = 0;
end

% ensure nonnegative integer for index_t0
mustBeNonnegative( index_t0 );
mustBeInteger( index_t0 );
mustBeScalarOrEmpty( index_t0 );

% ensure existence of nonempty window
if nargin < 10 || isempty( window )
    window = windows.tukey( 0.2 );
end

% ensure existence of nonempty F_number
if nargin < 11 || isempty( F_number )
    F_number = 1;
end

% ensure nonnegative integer for index_t0
mustBeNonnegative( F_number );
mustBeScalarOrEmpty( F_number );

%--------------------------------------------------------------------------
% 2.) geometry
%--------------------------------------------------------------------------
% number of array elements
N_elements = size( data_RF, 2 );
M_elements = ( N_elements - 1 ) / 2;

% half-width of the elements
element_width_over_two = element_width / 2;

% centroids of element faces
positions_ctr_x = (-M_elements:M_elements) * element_pitch;
positions_lbs_x = positions_ctr_x - element_width_over_two;
positions_ubs_x = positions_ctr_x + element_width_over_two;

% propagation direction
e_theta_x = sin( theta_incident );
e_theta_z = cos( theta_incident );

% reference position
position_ctr_x_ref = sign( e_theta_x ) * positions_ctr_x( 1 );

% lateral distances [ N_pos, N_elements ]
dist_lateral = positions_ctr_x - positions( :, 1 );
indicator_lower = dist_lateral < 0;

% incident wave travel times
t_in_plus_shift = ( e_theta_x * ( positions( :, 1 ) - position_ctr_x_ref ) + e_theta_z * positions( :, 2 ) ) / c_avg + index_t0 / f_s;

% scattered wave travel times
t_sc = sqrt( dist_lateral.^2 + positions( :, 2 ).^2 ) / c_avg;
tof = t_in_plus_shift + t_sc;

% maximum relative time shift in the electronic focusing ( zero padding in DFT )
positions_z_min = min( positions( :, 2 ) );
N_samples_t_add = ceil( ( sqrt( ( 2 * M_elements * element_pitch )^2 + positions_z_min^2 ) - positions_z_min ) * f_s / c_avg );

%--------------------------------------------------------------------------
% 3.) create frequency axis
%--------------------------------------------------------------------------
% number of samples in time domain
N_samples_t = size( data_RF, 1 );

% ensure odd number of points in DFT
N_points_dft = N_samples_t + N_samples_t_add + 1 - mod( N_samples_t + N_samples_t_add, 2 );

% boundary frequency indices
index_Omega_lb = ceil( f_bounds( 1 ) * N_points_dft / f_s ) + 1;
index_Omega_ub = floor( f_bounds( 2 ) * N_points_dft / f_s ) + 1;
indices_Omega = (index_Omega_lb:index_Omega_ub).';

% frequency axes
axis_f_bp = f_s * ( indices_Omega - 1 ) / N_points_dft;
axis_omega_bp = 2 * pi * axis_f_bp;

%--------------------------------------------------------------------------
% 4.) frequency-dependent F-number
%--------------------------------------------------------------------------
% maximum F-numbers for each axial position
% F_ub = sqrt( axis_f_bp .* positions( :, 2 ).' / ( 2.88 * c_avg ) );

%--------------------------------------------------------------------------
% 5.) electronic focusing
%--------------------------------------------------------------------------
% compute DFT of RF data
data_RF_dft = fft( data_RF, N_points_dft, 1 );
data_RF_dft_analy = 2 * data_RF_dft( indices_Omega, : );

% desired half-widths of the receive subapertures
width_aperture_over_two_desired = positions( :, 2 ).' ./ ( 2 * F_number );

% initialize image w/ zeros
image = complex( zeros( size( positions, 1 ), 1 ) );

% iterate voxel positions
for index_pos = 1:size( positions, 1 )

    % print progress in percent
    if mod( index_pos, 200 ) || index_pos == size( positions, 1 )
        fprintf( '%5.1f %%', ( index_pos - 1 ) / size( positions, 1 ) * 1e2 );
    end

    %----------------------------------------------------------------------
    % a) determine actual receive aperture
    %----------------------------------------------------------------------
    % map desired bounds of the receive aperture to element indices
    index_aperture_lb = max( ceil( M_elements + ( positions( index_pos, 1 ) - width_aperture_over_two_desired( :, index_pos ) ) / element_pitch ) + 1, 1 );
    index_aperture_ub = min( floor( M_elements + ( positions( index_pos, 1 ) + width_aperture_over_two_desired( :, index_pos ) ) / element_pitch ) + 1, N_elements );

    % no aperture
    if index_aperture_ub < 1 || index_aperture_lb > N_elements

        % erase progress in percent
        if mod( index_pos, 200 ) || index_pos == size( positions, 1 )
            fprintf( '\b\b\b\b\b\b\b' );
        end

        % skip calculations
        continue;

    end % if index_aperture_ub < 1 || index_aperture_lb > N_elements

    % actual width of the receive aperture
    width_aperture_lower_over_two = positions( index_pos, 1 ) - positions_lbs_x( index_aperture_lb ).';
    width_aperture_upper_over_two = positions_ubs_x( index_aperture_ub ).' - positions( index_pos, 1 );

    %----------------------------------------------------------------------
    % b) compute phase shifts
    %----------------------------------------------------------------------
    weights = exp( 1j * axis_omega_bp * tof( index_pos, : ) );

    %----------------------------------------------------------------------
    % c) compute apodization weights for current voxel
    %----------------------------------------------------------------------
    % check type of window function
    if isa( window, 'windows.boxcar' )

        %--------------------------------------------------------------
        % i.) simple solution for boxcar window
        %--------------------------------------------------------------
        window_samples = compute_samples( window, dist_lateral( index_pos, : ) ./ width_aperture_over_two_desired( :, index_pos_z ) );

    else

        %--------------------------------------------------------------
        % ii.) complex solution for other windows
        %--------------------------------------------------------------
        % sample window functions
        indicator_lower_act = indicator_lower( index_pos, : );
        window_samples_lower = compute_samples( window, dist_lateral( index_pos, indicator_lower_act ) ./ width_aperture_lower_over_two );
        window_samples_upper = compute_samples( window, dist_lateral( index_pos, ~indicator_lower_act ) ./ width_aperture_upper_over_two );
        window_samples = [ window_samples_lower, window_samples_upper ];

    end % if isa( window, 'windows.boxcar' )

    % normalize window samples
    window_samples = window_samples / sum( window_samples, 2 );

    %------------------------------------------------------------------
    % d) apply weights and focus RF signals
    %------------------------------------------------------------------
    data_RF_dft_analy_focused = sum( window_samples .* weights .* data_RF_dft_analy, 2 );

%         figure( index_pos_z + 3 );
%         temp = zeros( index_Omega_ub, N_elements );
%         temp( indices_Omega, : ) = window_samples .* weights;
%         window_samples_td = ifft( temp, N_points_dft, 1, 'symmetric' );
%         imagesc( illustration.dB( window_samples_td, 20 ), [ -60, 0 ] );
%         temp = zeros( index_Omega_ub, 1 );
%         temp( indices_Omega ) = data_RF_dft_analy_focused;
%         temp = ifft( temp, N_points_dft, 1, 'symmetric' );
%         plot( (1:N_samples_t), temp, (1:N_samples_t), abs( hilbert( temp ) ) );

    %------------------------------------------------------------------
    % d) inverse DFT of focused signal at time index 0
    %------------------------------------------------------------------
    image( index_pos ) = sum( data_RF_dft_analy_focused, 1 );

	% erase progress in percent
    if mod( index_pos, 200 ) || index_pos == size( positions, 1 )
        fprintf( '\b\b\b\b\b\b\b' );
    end

end % for index_pos = 1:size( positions, 1 )

% infer and print elapsed time
time_elapsed = toc( time_start );
fprintf( 'done! (%f s)\n', time_elapsed );

end % function [ image, time_elapsed ] = das_pw( positions, data_RF, f_s, theta_incident, element_width, element_pitch, c_avg, f_bounds, index_t0, window, F_number )
