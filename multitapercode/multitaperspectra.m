function spec = multitaperspectra(frames, tapers, weights, NFFT)

% function spec = multitaperspectra(frames, tapers, weights, NFFT)
%
% Compute multiple window (or multitaper) power spectra of multiple frames.
%
% Inputs: 
%   frames      Matrix of size (num_frames x frame_size), obtained e.g. 
%               using the 'buffer' command or 'enframe' in VoiceBox toolbox.
%               ***IMPORTANT*** the frames should NOT be windowed (such as
%               Hamming what one typically does), just give the 'raw' frames.
%               as input to this function.
%   tapers      Matrix of size (frame_size x num_windows) containing the 
%               tapers (window functions) as column vectors
%   weights     Vector of length num_windows of taper weights
%   NFFT        Number of FFT bins.
%
% Outputs:
%   spec        Power spectra estimates for each frame, one frame per
%               column.

spec = zeros(NFFT, size(frames', 2));
for (taper_nbr = 1:size(tapers,2))
    spec = spec + weights(taper_nbr) * abs(fft((frames').*repmat(tapers(:,taper_nbr),1,size(frames',2)), NFFT)).^2;
end
spec = spec(1:NFFT/2+1, :);