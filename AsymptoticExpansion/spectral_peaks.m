function peak_positions = spectral_peaks(mode_number, n_ext, r, n_bead, Modestart)
% Uses the size parameter to estimate peak positions in wavelength (nm)

r = r * 1e-6; % working in microns

%lookup code
m = n_bead / n_ext;
B = pi * r * 2 * 1e9; % in nanometers

% generates modes of the same mode number 
if Modestart == 1;
   
    RSP_TE = wgm_schiller(mode_number, m,  'TE');
    lambda_te = (1./RSP_TE) * B * (n_ext) ;
    
    RSP_TM = wgm_schiller(mode_number, m,  'TM');
    lambda_tm = (1./RSP_TM) * B* (n_ext);
    
end

% if mode numbers are different but consecutive  
if Modestart ==0;
    
    for i=1:length(mode_number)
        RSP_TE(i) = wgm_schiller(mode_number(i)+1, m, 'TE');
        RSP_TM(i) = wgm_schiller(mode_number(i), m, 'TM');
    end
    lambda_te = (1./RSP_TE) * B * (n_ext);
    lambda_tm = (1./RSP_TM) * B* (n_ext);

end

peak_positions = [lambda_te, lambda_tm];
peak_positions = sort(peak_positions, 'ascend');

end

