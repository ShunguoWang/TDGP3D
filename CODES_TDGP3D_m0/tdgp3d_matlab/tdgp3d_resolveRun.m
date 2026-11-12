function [S, F] = tdgp3d_resolveRun(runDir, o)
% TDGP3D_RESOLVERUN  (helper) Settings and file names of a 3D TDGP run folder.
%
% runDir holds tdgp.settings, MT3D_V6.dat (the forward mesh), the ModEM data file
% and the output files <outputBase>_ensemble.txt and <outputBase>_domain.txt.
% o is the option struct of the calling function; its fields settingsFile,
% ensembleFile, outDir and prefix override the defaults when not empty.
%
% Output S: the settings (tdgp3d_readSettings), with the domain the run used
% (<outputBase>_domain.txt: xmin ... zmax, lambda_increase, logrho_bg) overlaid.
% Output F: file names
%   F.settings  <runDir>/tdgp.settings
%   F.domain    <runDir>/<outputBase>_domain.txt   ('' if not there)
%   F.ensemble  <runDir>/<outputBase>_ensemble.txt
%   F.mesh      <runDir>/MT3D_V6.dat
%   F.data      <runDir>/<dataFile>                ('' if not there: no sites on the plots)
%   F.outDir    runDir
%   F.prefix    outputBase

if isempty(o.settingsFile)
    F.settings = fullfile(runDir, 'tdgp.settings');
else
    F.settings = o.settingsFile;
end
S = tdgp3d_readSettings(F.settings);

F.domain = fullfile(runDir, [S.outputbase '_domain.txt']);
if exist(F.domain, 'file') == 2
    S = tdgp3d_readSettings(F.domain, S);
else
    F.domain = '';
end

if isempty(o.ensembleFile)
    F.ensemble = fullfile(runDir, [S.outputbase '_ensemble.txt']);
else
    F.ensemble = o.ensembleFile;
end
F.mesh = fullfile(runDir, 'MT3D_V6.dat');
F.data = S.datafile;
if exist(F.data, 'file') ~= 2
    F.data = fullfile(runDir, S.datafile);
end
if exist(F.data, 'file') ~= 2
    F.data = '';
end
if isempty(o.outDir)
    F.outDir = runDir;
else
    F.outDir = o.outDir;
end
if isempty(o.prefix)
    F.prefix = S.outputbase;
else
    F.prefix = o.prefix;
end

if exist(F.ensemble, 'file') ~= 2
    error('tdgp3d_resolveRun:noEnsemble', 'Ensemble file not found: %s', F.ensemble);
end
if exist(F.mesh, 'file') ~= 2
    error('tdgp3d_resolveRun:noMesh', 'Forward mesh file not found: %s (the MT3D_V6.dat the run used)', F.mesh);
end

end
