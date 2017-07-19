function [f,t,cl,sc] = npd_tw_010916(freqFour,ij,seg_samp_no)
% function [f,t,cl,sc] = sp2_fn2_R2a(d1,d2,seg_samp_no,samp_rate,flags);
%
% Function with core routines for periodogram or multi-taper spectral estimates,
%  with R2 and directionality analysis.
% see NOTE below about calling this routine.
%
% Copyright (C) 2015, David M. Halliday.
% This file is part of NeuroSpec.
%
%    NeuroSpec is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    NeuroSpec is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with NeuroSpec; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%    NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org
%
%  Inputs are two matrices containing pre processed time series or point process data.
%  Matrices have L columns with T rows.
%   L = Number of segments: seg_tot.
%   T = DFT segment length: seg_size.
% 
% Input arguments
%  d1          Channel 1 data matrix.
%  d2          Channel 2 data matrix.
%  seg_samp_no Vector with number of data points in each segment.
%  samp_rate   Sampling rate (samples/sec)
%  flags       Structure with flags to control processing options.
%              Flags supported:
%                wind: Apply data window (split cosine taper) to raw data.
%                        (0:No window; 1,2,3,...: Type).  1:10%; 2:20%; 3:50%; 4:100% taper.
%                 MT:  Multi-taper option. No of tapers is MT_tap.
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates.
%  t    column matrix with      time domain parameter estimates.
%  cl   single structure with scalar values related to analysis.
%  sc   optional matrix for spectral coefficients.
%
% Output parameters
%  f column 1       frequency in Hz.
%  f column 2       Log  input/d1 spectrum.
%  f column 3       Log output/d2 spectrum.
%  f column 4       Coherence.
%  f column 5       Phase.
%  f column 6 - 12  Frequency domain directionality analysis
%
%  t column 1       Lag in ms.
%  t column 2       Cumulant density.
%  t column 3       Time domain directionality analysis, rho estimate
%
%  cl               Structure with values for confidence limits and other analysis parameters.
%
%  sc column 1      f11.
%  sc column 2      f22.
%  sc column 3      f21 (Complex).
%
% References:
% 1. Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
%     Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
% 2. Bloomfield, P. Fourier Analysis of Time Series: An Introduction.
%    2nd edition. Wiley, New York, 2000.
% 3. Halliday DM (2015) Nonparametric directionality measures for time series and point process data,
%     Journal of Integrative Neuroscience, 14(2), In Press. DOI: 10.1142/S0219635215300127
%
% References 1 and 3 referred to in comments as PBMB and JIN, respectively
%
% function [f,t,cl,sc] = sp2_fn2_R2a(d1,d2,seg_samp_no,samp_rate,flags);
%
% NOTE: This routine is not intended to support analysis of raw data.
% It is intended as a support routine for the 2 channel spectral analysis function:
% sp2a2_R2_mt.m. Refer to this function for further details.

%%%%
fd1 = squeeze(freqFour.fourierspctrm(:,ij(1),2:end))';
fd1 = [fd1;fd1(end:-1:1,:)];
fd2 = squeeze(freqFour.fourierspctrm(:,ij(2),2:end))';
fd2 = [fd2;fd2(end:-1:1,:)];

samp_tot = numel(fd1);
seg_tot = size(fd1,2);
samp_rate = max(freqFour.freq)*2;
seg_size = size(fd1,1);
flags.MT_tap = freqFour.cumtapcnt(1); 
psd_fac = 1/(2*pi*samp_tot);

f11=psd_fac*sum(abs(fd1.*fd1),2);   % Spectrum 1, PBMB (5.2), Mag squared for  input auto spectra.
f22=psd_fac*sum(abs(fd2.*fd2),2);   % Spectrum 2, PBMB (5.2), Mag squared for output auto spectra.
f21=psd_fac*sum(fd2.*conj(fd1),2);  % Cross spectrum (complex valued), PBMB (5.2).

% Call function sp2_fnR2b for directionality analysis
[f_R2,t_R2,cl_R2]=sp2_fnR2b(fd1,fd2,f11,f22,psd_fac,samp_tot);

% Coherence and R2 value, coherence at zero frequency assumed zero
chyx=[0; abs(f21(2:seg_size).*f21(2:seg_size))./(f11(2:seg_size).*f22(2:seg_size))];
R2_ch=(1/seg_size)*sum(chyx);   % Integral from -pi to +pi, estimate for JIN (2.2)

% Construct output spectral matrix f.
f_index = (2:seg_size/2+1)';    % Indexing for output, DC component not output.
% f_index=(1:seg_size)';    % Indexing for output, DC component not output.
deltaf=samp_rate/seg_size;    % Spacing of Fourier frequencies in Hz.
f(:,1)=(f_index-1)*deltaf;    % Column 1 - frequencies in Hz.
f(:,2)=log10(f11(f_index));   % Column 2 - Log spectrum ch 1.
f(:,3)=log10(f22(f_index));   % Column 3 - Log spectrum ch 2.
f(:,4)=chyx(f_index);         % Column 4 - Coherence, PBMB (5.5).
f(:,5)=angle(f21(f_index));   % Column 5 - Phase, PBMB (5.7).

% Append R2 results, 7 columns in matrix f_R2 from function sp2_fnR2b
f=[f f_R2];

% Construct output time domain matrix t.
% Column 1 of t matrix is time in msec. Range (-T/2)*dt to (T/2-1)*dt.
deltat=1000.0/samp_rate;    % dt in msec.
t(:,1)=((1:seg_size)'-seg_size/2-1)*deltat;
% Estimate cumulant density using inverse DFT of cross spectrum.
[t(:,2),q_var]=calc_q21(f11(1:seg_size/2+1),f22(1:seg_size/2+1),f21(1:seg_size/2+1),seg_size,samp_tot);

% Append rho estimate from function sp2_fnR2b
t=[t t_R2];

% Construct optional output spectral matrix sc, DC value is output.
if (nargout>3)
  sc(:,1)=f11(1:seg_size/2+1);  % Column 1 - f11
  sc(:,2)=f22(1:seg_size/2+1);  % Column 2 - f22
  sc(:,3)=f21(1:seg_size/2+1);  % Column 3 - f21 (Complex valued)
end

% Calculate effective number of segments, L', used to approximate variance of estimates.
seg_tot_var=1/sum((seg_samp_no(:,1)/samp_tot).*(seg_samp_no(:,1)/samp_tot));
% Correct for multi-taper option
  seg_tot_var=seg_tot_var*flags.MT_tap;


% Construct cl structure, confidence limits for parameter estimates.
cl=cl_R2;                    % Include values returned by sp2_fnR2b
cl.col_R20=10;               % column in f matrix containing zero lag coherence component
cl.col_R2p=11;               % column in f matrix containing forward (positive lag) coherence component
cl.col_R2n=12;               % column in f matrix containing reverse (negative lag) coherence component
cl.col_rho=3;                % column in t matrix containing rho estimate
cl.R2_ch=R2_ch;              % R2 estimated from coherence between original processes
% cl.type=flags.sp_type;       % Analysis type.
cl.seg_size=seg_size;        % T.
cl.seg_tot=seg_tot;          % L.
cl.seg_tot_var=seg_tot_var;  % Effective no of segments (L').
cl.samp_tot=samp_tot;        % R.
cl.samp_rate=samp_rate;      % Sampling rate.
cl.dt=deltat;                % Delta t.
cl.df=deltaf;                % Delta f.
cl.f_c95=0.8512*sqrt(1/seg_tot_var);  % 95% Confidence limit for spectral estimates, PBMB (6.2)
                                      % Confidence interval in plot is TWICE this value
cl.ch_c95=1-0.05^(1/(seg_tot_var-1)); % 95% Confidence limit for coherence, PBMB (6.6).
cl.q_c95=1.96*sqrt(q_var);            % 95% Confidence limits for cumulant, PBMB (6.11).

% % Extra values associated with multi-taper option
% if flags.MT
%   cl.MT=flags.MT;          % Flag indicating multi-taper (MT) estiamte - 0: No, 1: Yes
%   cl.MT_NW=flags.NW;       % NW for MT estimates
%   cl.MT_tap=flags.MT_tap;  % No of tapers in MT estimate
%   cl.MT_eig=MT_eig;        % Eigenvalues for individual eigen spectra in MT estimate  
% else
%   cl.MT=flags.MT;
%   cl.MT_NW=0;
%   cl.MT_tap=0;
%   cl.MT_eig=0;
% end

return

%------------------------------------------------------------------------------
function [q21,q21_var] = calc_q21(f11,f22,f21,seg_size,samp_tot);
%function [q21,q21_var] = calc_q21(f11,f22,f21,seg_size,samp_tot);
%
% Function to estimate cumulant density from complex cross spectrum.
% Also includes estimate of variance.

% f21 is (1:seg_size/2+1)
% Expand to hermitian vector.
f21_all=[f21; conj(f21(seg_size/2:-1:2))];

% Inverse DFT.
cov=ifft(f21_all);

% Estimate cumulant, shifted by T/2 so that time zero is in centre.
% 2pi/T factor is 2*pi, since ifft routine includes 1/T term.
q21([seg_size/2+1:seg_size,1:seg_size/2],1)=real(cov(1:seg_size))*2*pi; % PBMB (5.9).

% Estimate variance of cumulant density estimate.
var_fac=4*pi*pi/(seg_size*samp_tot);                             % Factor (2pi/T)(2pi/R).
q21_var=var_fac*2*sum(f11(1:seg_size/2+1).*f22(1:seg_size/2+1)); % PBMB (6.10).

return
