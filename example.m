% Example HiLo analysis session

% Import data from CSV file. Headers are automatically parsed and used as field names
% dAU is a MATLAB structure array with fields corresponding to CSV file headers

data = ImportCSV('data.csv');

% translateAuburnData is a helper script that transfers the fields actually needed by 
% the hilo.m script into the field names that it expects. 

d0 = TranslateData(data, 'crop', 1);

% Check for depletion and accumlation offset between HF and QS capacitance. If 
% gain/offset correction is warranted select appropriate range (strong depletion and 
% strong accumulation)

plot(d0.Vg, d0.Cqs - d0.Chf)
rng = find((d0.Vg < -2.0) | (d0.Vg > 3.0));

% Perform analysis with specified area, with series resistance and gain/offset 
% correction enabled

d = hilo(d0, 'Area', 2.827e-3, 'RsCorrect', true, 'AdjRange', rng, 'Debug', true);

% Output:
%   Cox = 150.941 pF
%   Rs = 111.45 ohms
%   Nd = 4.816e+15 cm^-3
%   Vfb = 0.74 V
%   Ec-Ef (flatband) = 0.203 eV
%   Nit (HiLo) = 5.92e+10 cm^-2
%   Nit (Cpsi) = 1.89e+11 cm^-2

% Save results to CSV file

SaveAsCSV(d,['processed/' d.fname]);
