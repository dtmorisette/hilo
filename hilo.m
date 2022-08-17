function d_out = hilo(d, varargin)
% HILO  Process high-low C-V data into interface trap density vs. energy
%
% Notes:
%   * Has not yet been tested on p-type data; should work but...
%
% Syntax
%
%   d = HILO('-v')  prints version string
%   d = HILO(d)     processes the given CV data with default parameters
%   d = HILO(d, 'param1', value1, 'param2', value2, ...) 
%                   allows additional parameters to be specified. 
%
%   d is a struct containing all the data required to perform HiLo CV
%   analysis. The minimum set of fields is:
%
%       d.Vg   vector   Gate voltage (V)
%       d.Cqs  vector   Quasi-static capacitance (F)
%       d.Chf  vector   HF capacitance (F)
%       
%   Optional fields:
%
%       d.f    scalar   Frequency for HF capacitance measurement (Hz)
%       d.G    vector   Parallel model measured conductance
%
%   Other fields are ignored. All vector fields must have the same length.
%
%
% Allowable parameters are:
%
%   'Temperature' (default: 21°C)
%       Specifies the temperature at which the measurent was performed. 
%       If specified, any existing value is overwritten. If not specified 
%       and no existing value is present, the default value is used.
%
%   'Substrate' (default: 'SiC')
%       Name of the material used for the substrate. One of 'Si' or 'SiC'
%
%   'Oxide' (default: 'SiO2')
%       Name of the material used for the gate insulator.
%       One of 'SiO2' or 'Al2O3'
%
%   'Area' (default: 1.59e-3 cm^2)
%       The area of the measured device. Default is area of 450 um circle
%       If specified, any existing value is overwritten. If not specified 
%       and no existing value is present, the default value is used.
%
%   'Frequency' (default: 100 kHz)
%       The frequency at which the HF capacitance data was taken. 
%       If specified, any existing value is overwritten. If not specified 
%       and no existing value is present, the default value is used.
%
%   'Cox' 
%       The oxide capacitance. (default: calculated by averaging 5 data 
%       points near the maximum in the quasi-static capacitance)
%       If specified, any existing value is overwritten. If not specified 
%       and no existing value is present, the default value is used.
%
%   'Doping' (default: auto-calculated)
%       The doping level of the substrate. Specifying this option overrides
%       the default calculated from the slope of 1/CHF^2 vs Vg. Positive
%       values indicate p-type, negative values indicate n-type.
%
%   'Offset' (default: 0)
%       The offset applied to the calculated surface potential, specified
%       in units of kT/q. Used to explore the effect of arbitrary shifts 
%       in the surface potential.
%
%   'Range' (default: [])
%       A list of indices of the data that will be retained for analysis.
%       All other data in the vectors supplied in 'd' will be discarded.
%       An empty array [] indicates all data should be used.
%
%   'RsCorrect' (default: true)
%       Use the (optionally) provided conductance data (d.G) to correct
%       both Chf and G data for series resistance. Original data is moved
%       to d.Chf_raw and d.G_raw
%
%   'AdjRange' (default: [])
%       A list of indices of the data to use to calculate gain and offset
%       correction parameters. If the list is empty no offset correciton is
%       performed. The listed data points should be in an range where 
%       where Cqs and Chf should be equal, such as strong depletion. The 
%       adjusted Chf overwrites the existing Chf data.
%
%   'AdjOffsetOnly' (default: true)
%       If False, the corrected d.Chf' = A * d.Chf + B. The coeffients 
%       A and B are found to minimize the error term Cqs-Chf in a 
%       least-squared sense. If 'Range' is used in conjunction with 
%       'AdjRange' the latter is corrected as needed. 
%       If True, A is assumed to be 1, or equivalenty no gain correction
%       is applied. The average offset in the supplied range is removed.
%
%   'Debug' (default: true)
%       The debug mode generates diagnostic plots and prints diagnostic
%       information to the command line.
%
%   'Method' (default: 'Intercept')
%       Selects the method that will be used to calculate the surface
%       potential. Options are:
%
%       'FullFit'   Fits the entire 1/Cs_hf^2 dataset to the exact 
%                   MOS-CV model
%       'Intercept' Fits the linear section of the 1/Cs_hf^2 data to a
%                   linear model and forces the x-intercept to -1 kT/q
%
% Data added to structure by HILO:
%
%   The following are added only if not already present:
%
%       d.T         Measurement temperature (°C)
%       d.area      Device area (cm^2)
%       d.f         Frequency used during HF capacitance measurement (Hz) 
%
%   The following are added from the parameters given on the command line:
%
%       d.offset    The arbitrary offset added to the internally determined
%                   surface potential (in normalized units: kT/q)
%
%   The following are added during the analysis
%
%       d.Rs        Series resistance
%       d.Cox       Oxide capacitance (if not specified on command line)
%       d.invCsq    1/Chf^2 normalized to area (cm^4/F^2)
%       d.Nd        Donor concentration as extracted from 1/Chf^2 vs Vg
%       d.Na        Acceptor concentration as extracted from 1/Chf^2 vs Vg
%       d.w         Depth associated with extracted doping profile (cm)
%       d.doping    Doping profile (- for n-type, + for p-type) (cm^-3)
%       d.Ef_Ei     Bulk Ef-Ei (eV) i.e. doping parameter (-kT/q * Uf)
%       d.Lde       Extrinsic Debye length (cm)
%       d.Cs_fb     Semiconductor capacitance at flatband (F)
%       d.Cs_hf     HF semiconductor capacitance (F)
%       d.Cs_qs     QS semiconductor capacitance (F)
%       d.phi_s     Surface potential (V)
%       d.Ec_Ef     Ec-Ef at surface (trap position relative to Ec) (eV)
%       d.Vfb       Flatband voltage (V)
%       d.Cfb       Flatband capacitance (V)
%       d.Ecb_Ef    Flatband energy (on Ec-Ef scale used for Dit)
%       d.Cs_ideal  Ideal semiconductor depletion capacitance (F)
%       d.Chf_ideal Ideal HF capacitance (F)
%       d.Cfb_ideal Ideal flatband capacitance (F)
%       d.Dit_hilo  Dit by traditional Hi-Lo technique (cm^-2 eV^-1)
%       d.Dit_Cpsi  Dit by C-psi technique (cm^-2 eV^-1)
%       d.Nit_hilo  Nit (Hi-Lo Dit integrated over Ec-Ef) midgap to flatband
%       d.Nit_Cpsi  Nit (C-psi Dit integrated over Ec-Ef) midgap to flatband
%
% References:
%
%   C. N. Berglund, ?Surface states at steam-grown silicon-silicon dioxide 
%       interfaces,? Electron Devices, IEEE Transactions on, vol. 13, 
%       no. 10, pp. 701?705, 1966.
%
%	J. R. Brews, ?Correcting interface-state errors in MOS doping profile 
%       determinations,? Journal of Applied Physics, vol. 44, no. 7, 
%       pp. 3228?3231, 1973.
%
%   H. Yoshioka, T. Nakamura, and T. Kimoto, ?Accurate evaluation of 
%       interface state density in SiC metal-oxide-semiconductor structures
%       using surface potential based on depletion capacitance,?
%       Journal of Applied Physics, vol. 111, no. 1, pp. 014502 1?5, 2012.
%
% Copyright (c) 2013 Dallas T. Morisette (morisett@purdue.edu).
% Released under the terms of the FreeBSD License. 
% See LICENSE file for details.
%

    version = '$VERSION_STRING$';

    ip = inputParser;
    ip.addRequired('d', @(x)isstruct(x) || ischar(x));
    ip.addParamValue('Temperature',    -1, @isscalar);
    ip.addParamValue('Substrate',   'SiC'           );
    ip.addParamValue('Oxide',      'SiO2'           );
    ip.addParamValue('Area',           -1, @isscalar);
    ip.addParamValue('Frequency',      -1, @isscalar);
    ip.addParamValue('Cox',            -1, @isscalar);
    ip.addParamValue('Doping',          0, @isscalar);
    ip.addParamValue('Offset',          0, @isscalar);
    ip.addParamValue('Range',          [], @isvector);
    ip.addParamValue('RsCorrect',    true           );
    ip.addParamValue('AdjRange',       [], @isvector);
    ip.addParamValue('AdjOffsetOnly',true           );
    ip.addParamValue('Debug',        true           );
    ip.addParamValue('Method','Intercept'           );
    
    ip.parse(d, varargin{:});
    
    if ischar(d)
        switch d
            case {'-v', '--version'}
                disp(version);
            otherwise
                disp(['Unknown option: ' d]);
        end
        return;
    end                
    
    debug = ip.Results.Debug;
    method = lower(ip.Results.Method);
    
    VerifyLengths(d);
    nFig = 1;
    
    % If the 'Temperature' parameter is specified, override any existing
    % value in d.T. If not specified, and if d.T does not exist, set 
    % to default of 21°C
    
    if ip.Results.Temperature ~= -1
        d.T = ip.Results.Temperature;
    elseif ~isfield(d, 'T')
        d.T = 21;
    end
    
    % Get fundamental physical constants and material properties
    c  = GetConstants(d.T);
    s  = GetMaterial(ip.Results.Substrate, d.T);
    
    % If the 'Area' parameter is specified, override any existing
    % value in d.area. If not specified, and if d.area does not exist, set 
    % to default of 1.59e-3 cm^2 (area of 450um diameter circle)
    
    if ip.Results.Area ~= -1
        d.area = ip.Results.Area;
    elseif ~isfield(d, 'area');
        d.area = ((450e-4)^2)*pi/4;
    end
    
    % If the 'Frequency' parameter is specified, override any existing
    % value in d.f. If not specified, and if d.f does not exist, set 
    % to default of 100e3 (100 kHz)

    if ip.Results.Frequency ~= -1
        d.f = ip.Results.Frequency;
    elseif ~isfield(d, 'f');
        d.f = 100e3;
    end
    
    d.offset = ip.Results.Offset;
    
    cropRange = ip.Results.Range;
    adjustRange = ip.Results.AdjRange;
    
    % If a range is supplied, crop data to range
    if ~isempty(cropRange)
        n = length(d.Vg);
        mask1 = ismember(1:n,cropRange);
        d = CropRange(d, cropRange);
        if ~isempty(adjustRange)
            mask2 = ismember(1:n,adjustRange);
            mask = mask1 & mask2;
            adjustRange = find(mask(cropRange));
        end
    end    
    
    % If RsCorrect was requested, extract Rs and correct d.Chf and d.G
    % Saves old Chf and G data in d.Chf_raw and d.G_raw, respectively
    % Will only be performed if d.Chf_raw does not exist, since running
    % the correction more than once would corrupt the raw data and would
    % not improve the correction (extracted Rs should be close to zero if
    % the correction has already been performed).
    
    if ip.Results.RsCorrect && ~isfield(d, 'Chf_raw')
        d.Chf_raw = d.Chf;
        d.G_raw   = d.G;

        d.Rs = FindRs(d.Chf, d.G, d.f, 5);

        den = (d.G*d.Rs - 1).^2 + (2*pi*d.f*d.Chf*d.Rs).^2;
        d.Chf = d.Chf./den;
        d.G = (d.G - d.Rs*d.G.^2 - d.Rs*(2*pi*d.f*d.Chf).^2)./den;
        if (debug)
            fprintf('Rs = %.2f ohms\n', d.Rs);
        end
    end
      
    % Perform gain/offset correction
    % Chf' = A * Chf + B
    if ~isempty(adjustRange)
        if ip.Results.AdjOffsetOnly % A = 1
            d.Chf = d.Chf - mean(d.Chf(adjustRange)-d.Cqs(adjustRange));
        else
            pf = [d.Chf(adjustRange) ones(length(d.Chf(adjustRange)),1)]\d.Cqs(adjustRange);
            d.Chf = d.Chf*pf(1)+pf(2);
        end
    end
    
    % If the 'Cox' parameter is specified, override any existing value in
    % d.Cox. If not specified, and if d.Cox does not exist, calculate based 
    % on the average of 5 points in strong accumulation.
    
    if ip.Results.Cox > 0
        d.Cox = ip.Results.Cox;
    elseif ~isfield(d, 'Cox');
        d.Cox = FindCox(d.Cqs, 5);
    end
    if debug
        fprintf('Cox = %.3f pF\n', d.Cox*1e12);
    end
    
    % Extract doping using slope of linear fit of 1/Chf^2 vs Vg plot
    d.invCsq = 1./(d.Chf/d.area).^2;
    [pf, rng] = FindLinearFit(d.Vg, d.invCsq, 5);
    if ip.Results.Doping == 0
        doping = 2/(c.q*c.eps0*s.k*pf(1));
    else
        doping = ip.Results.Doping;
    end
    
    if doping < 0
        d.Nd = -doping;
        d.Na = 0;
        d.type = 'n';
    else
        d.Na = doping;
        d.Nd = 0;
        d.type = 'p';
    end
    
    % Calculate doping profile
    % J. R. Brews, J. of Appl. Phys., vol. 44, no. 7, pp. 3228-3231, 1973.
    d.w = s.k*c.eps0*d.area*(1./d.Chf - 1/d.Cox);
    d.doping = 2/(c.q*s.k*c.eps0)*(1-d.Cqs/d.Cox)./(1-d.Chf/d.Cox) ./ ...
               dydx(d.Vg,d.invCsq,11);
    
    if debug
        f = figure(nFig);
        f.Name = 'InvCsq_vs_VG';
        nFig = nFig + 1;
        
        if d.type == 'n'
            x = linspace(min(d.Vg),-pf(2)/pf(1));
        else
            x = linspace(-pf(2)/pf(1), max(d.Vg));
        end
        plot(d.Vg, d.invCsq, 'ob', ...
             d.Vg(rng), d.invCsq(rng), 'og', ...
             x, polyval(pf,x), 'r');
        xlabel('Gate Voltage (V)');
        ylabel('1/C_{HF}^2 (cm^4/F^2)');
        grid on;
        
        if d.type == 'n'
            str = 'Nd';
        else
            str = 'Na';
        end
        fprintf('%s = %.3e cm^-3\n', str, abs(doping));

        f = figure(nFig);
        f.Name = 'doping_profile';
        nFig = nFig + 1;
        
        semilogy(d.w*1e7, abs(d.doping));
        xlabel('Depth (nm)');
        if d.type == 'n'
            str = 'Donor';
        else
            str = 'Acceptor';
        end
        ylabel(sprintf('%s Concentration (cm^{-3})', str));
        grid on;
    end
    
    % Flatband capacitance and voltage
    Uf = -asinh((d.Nd - d.Na)/s.ni/2);
    d.Ef_Ei = -c.kT*Uf;
    nb = s.ni*exp( d.Ef_Ei/c.kT);
    pb = s.ni*exp(-d.Ef_Ei/c.kT);
    d.Lde = sqrt(s.k*c.eps0*c.kTq/(c.q*(pb+nb)));    
    d.Cs_fb = s.k*c.eps0/d.Lde*d.area;

    % Surface potential
    d.Cs_hf = d.Chf*d.Cox./(d.Cox-d.Chf);
    d.Cs_qs = d.Cqs*d.Cox./(d.Cox-d.Cqs);
    d.phi_s = cumtrapz(d.Vg,1-d.Cqs/d.Cox);
    
    switch method
        case 'fullfit'
            Us = d.phi_s/c.kTq;
            options = optimset('Display', 'off', 'TolX',1e-7);
            dU = fminbnd(@(dU)rms((d.Cs_fb./d.Cs_hf).^2 - 1./(CdNorm(Us-dU,Uf)).^2), ...
                          -5,5, options);
        case 'intercept'
            [pf, ~] = FindLinearFit(d.phi_s/c.kTq, 1./(d.Cs_hf).^2, 5);
            dU = -pf(2)/pf(1) + 1;
        
        otherwise
            dU = 0;
    end
    d.phi_s = d.phi_s - dU*c.kTq + d.offset*c.kTq;
    d.Ec_Ef = s.Ec_Ei - d.Ef_Ei - d.phi_s;      

    if debug
        invCssq = 1./(d.Cs_hf/d.area).^2;
        [pf, rng] = FindLinearFit(d.phi_s, invCssq, 5);
        
        f = figure(nFig);
        f.Name = 'InvCsq_vs_SurfacePotential';
        nFig = nFig + 1;
        
        if d.type == 'n'
            x = linspace(min(d.phi_s), -pf(2)/pf(1));
        else
            x = linspace(-pf(2)/pf(1), max(d.phi_s));
        end
        plot(d.phi_s/c.kTq, invCssq, 'ob', ...
             d.phi_s(rng)/c.kTq, invCssq(rng), 'og', ...
             x/c.kTq, polyval(pf,x), 'r');
        xlabel('Surface Potential (kT/q)');
        ylabel('1/C_{SHF}^2 (cm^4/F^2)');
        grid on;        
    end
    
    try
        % Flatband voltage interpolation
        index = find(abs(diff(sign(d.phi_s)))==2);
        rng = index-5:index+5;    
        d.Vfb = interp1(d.phi_s(rng),d.Vg(rng), 0);

        % Flatband capacitance interpolation
        index = find(abs(diff(sign(d.Vg-d.Vfb)))==2);
        rng = index-5:index+5;    
        d.Cfb = interp1(d.Vg(rng),d.Chf(rng),d.Vfb);

        if debug
            fprintf('Vfb = %.2f V\n', d.Vfb);
        end
    catch err
    end
    
    % Flatband location on energy scale
    d.Ecb_Ef = s.Ec_Ei-d.Ef_Ei;
    if debug
        fprintf('Ec-Ef (flatband) = %.3f eV\n', d.Ecb_Ef);
    end

    % Ideal capacitances (Cs, Chf, Cfb)
    d.Cs_ideal = d.Cs_fb*CdNorm(d.phi_s/c.kTq,Uf);
    d.Chf_ideal = d.Cox*d.Cs_ideal./(d.Cox+d.Cs_ideal);
    d.Cfb_ideal = d.Cox*d.Cs_fb/(d.Cox + d.Cs_fb);
    
    % Calculate Dit by Hi-Lo and C-psi
    d.Dit_hilo = (d.Cs_qs - d.Cs_hf   )/d.area/c.q;
    d.Dit_Cpsi = (d.Cs_qs - d.Cs_ideal)/d.area/c.q;

    % Calculate Nit. Ignores data from flatband to majority band edge.
    % Data in the near band edge region is likely unreliable, but as a 
    % result the calculated Nit is underestimated. 
    
    if (d.type == 'n')
        rng = find(d.Ec_Ef >= d.Ecb_Ef & d.Ec_Ef < s.Eg/2);
    else
        rng = find(d.Ec_Ef <= s.Eg - d.Ecb_Ef & d.Ec_Ef > s.Eg/2);
    end
    d.Nit_hilo = trapz(d.Ec_Ef(rng),d.Dit_hilo(rng));
    d.Nit_Cpsi = trapz(d.Ec_Ef(rng),d.Dit_Cpsi(rng));
    
    if debug
        f = figure(nFig);
        f.Name = 'Dit_vs_Energy';
        nFig = nFig + 1;
        
        semilogy(d.Ec_Ef, [d.Dit_hilo d.Dit_Cpsi],'.');
        xlabel('Trap Energy (E_c - E_t, eV)');
        ylabel('Interface Trap Density (cm^{-2} eV^{-1})');
        ylim([1e10 10^(ceil(log10(max([d.Dit_hilo; d.Dit_Cpsi]))))]);
        grid on;
        hold on;
        semilogy(d.Ecb_Ef*[1 1],ylim,'r--');
        hold off;
        
        f = figure(nFig);
        f.Name = 'Cs_vs_SurfacePotential';
        nFig = nFig + 1;
        
        semilogy(d.phi_s/c.kTq, [d.Cs_qs d.Cs_hf d.Cs_ideal]);
        xlabel('Surface Potential (kTq)')
        ylabel('Semiconductor Capacitance (F)');
        if d.type == 'n'
            xlim([-s.Eg*0.4 d.Ecb_Ef]/c.kT);
        else
            xlim([d.Ecb_Ef -s.Eg*0.4]/c.kT);
        end
        grid on;
        
        f = figure(nFig);
        f.Name = 'C_vs_VG';
        nFig = nFig + 1;
        
        plot(d.Vg, [d.Cqs d.Chf d.Chf_ideal]*1e12, ...
             d.Vg, d.Cox*ones(size(d.Vg))*1e12,'k--');
        if isfield(d,'Vfb')
             hold on;
             plot(d.Vfb, d.Cfb*1e12, '.g');
             hold off;
        end
        xlabel('Gate Voltage (V)');
        ylabel('Capacitance (pF)');    
        grid on;
        
        fprintf('Nit (HiLo) = %.2e cm^-2\n', d.Nit_hilo);
        fprintf('Nit (Cpsi) = %.2e cm^-2\n', d.Nit_Cpsi);
    end
    
    d_out = d;
end

function rng = FindAccumulation(C,n)
% FINDACCUMULATION(C,n)  Finds n points around the maximum of a smoothed 
%   version of the given capacitance data. Selected points are biased 
%   towards the nearest end of the data set. Used to locate the area of
%   strong accumulation to extract Rs and Cox

    % Idiot check - make sure n is not greater than length of data set
    if length(C) < n
        n = length(C);
    end
    
    % Smooth given data and locate maximum
    C = smooth(C);
    i = find(C == max(C));
    i = i(1);
    
    % Select data points to average
    m = length(C);
    if i < m/2          % max in low end of data set
        if i <= n
            rng = 1:n;
        else
            rng = (i-n+1):i;
        end
    else                        % max in high end of data set
        if i > m-n
            rng = (m-n+1):m;
        else
            rng = i:(i+n-1);
        end
    end
end

function Cox = FindCox(C,n)
    rng = FindAccumulation(C,n);
    Cox = mean(C(rng))*1.001; % Add 0.1% to eliminate NaN errors
end

function Rs = FindRs(C,G,f,n)
    rng = FindAccumulation(C,n);
    w2 = (2*pi*f)^2;
    Rs = G./(w2*C.^2 + G.^2);
    Rs = mean(Rs(rng));
end

function d = CropRange(d,rng)
% CROPRANGE(d,rng)  Crop all vectors in the struct d to the given range rng
%
    fields = fieldnames(d);
    for i = 1:length(fields)
        field = fields{i};
        tmp = d.(field);
        if isvector(tmp) && ~isscalar(tmp) && isnumeric(tmp)
            d.(field) = tmp(rng);
        end
    end
end

function VerifyLengths(d)
% VERIFYLENGTHS(d)  Ensure that all numeric vectors are the same length
% in the given structure d
%
    n = 0;
    fields = fieldnames(d);
    for i = 1:length(fields);
        tmp = d.(fields{i});
        if isvector(tmp) && ~isscalar(tmp) && isnumeric(tmp)
            if n == 0
                n = length(tmp);
            else
                assert(n == length(tmp), 'All vectors must be the same length');
            end
        end
    end
end

function C = CdNorm(Us,Uf)
% CDNORM  Calculate the normalized depletion capacitance from 
% the surface potential Us and doping parameter Uf.
% The normalization constant is the flatband semiconductor capacitance
%
% Uf = (Ei(bulk) - Ef)/kT = { ln(ni/ND) n-type or ln(NA/ni) p-type }
% Us = (Ei(bulk) - Ei(surface))/kT
%
    if Uf < 0   % n-type
        C =  sign(Us).*(exp( Us)-1)./sqrt(2*(exp( Us)-Us-1));
    else        % p-type
        C = -sign(Us).*(exp(-Us)-1)./sqrt(2*(exp(-Us)+Us-1));
    end
end
