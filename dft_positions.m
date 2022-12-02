function result = dft_positions( samples, positions, axis_k_x, axis_k_z )
%DFT_POSITIONS Compute discrete Fourier transform for arbitrary grids
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
%   01.) samples:           samples
%   02.) positions:         sample positions (m)
%   03.) axis_k_x:          angular spatial frequencies along the x-axis (rad / m)
%   04.) axis_k_z:          angular spatial frequencies along the z-axis (rad / m)
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%   01.) result:            Fourier transform of impulse-modulated bivariate function
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
%   date: 2022-09-13
%   modified: 2022-12-02

% phase shift
axis_k_x_times_x = positions( :, 1 ) .* axis_k_x;
axis_k_z_times_z = positions( :, 2 ) .* axis_k_z;

% complex exponentials
exp_x = exp( -1j * axis_k_x_times_x );
exp_z = exp( -1j * axis_k_z_times_z );

temp = samples .* exp_x;

% specify cell array for results
result = cell( numel( axis_k_z ), 1 );

% iterate axial frequencies
for index_k_z = 1:numel( axis_k_z )

    % compute sum
    result{ index_k_z } = sum( temp .* exp_z( :, index_k_z ), 1 );

end % for index_k_z = 1:numel( axis_k_z )

% concatenate
result = cat( 1, result{ : } );

end % function result = dft_positions( samples, positions, axis_k_x, axis_k_z )
