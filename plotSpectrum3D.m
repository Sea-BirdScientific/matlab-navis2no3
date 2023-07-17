function plotSpectrum3D(dat, cal)

% Create 3D plot: Wavelength (WL) x Pressure (P) x UV_INTENSITY
%
% We're looking for saturated spectra.

w = cal.WL;   % Wavelength
p = dat.P;    % Pressure
pStart =  dat.spectra_pix_range(1); % start of the pixel span of spectrometer channels
pEnd   = dat.spectra_pix_range(2);     % end   of the pixel span of spectrometer channels
R = pStart:pEnd;
wr = w(R);
spec = dat.UV_INTEN;   % UV Intensity (Counts)
disp(max(max(spec)));

surf(wr, p, spec);
ylabel('Pressure  (dbar)'); 
xlabel('Wavelength  (nm)');

az = -18.42;
el =  29.44;
view(az,el);