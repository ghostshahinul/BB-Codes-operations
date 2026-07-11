function setup_BB_project()
%SETUP_BB_PROJECT Add BB-code project folders to the MATLAB path.

    projectRoot = fileparts(mfilename('E:\PhD Life\LDI project\BB-code-git\BB-Codes-operations'));

    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'scripts'));
    addpath(fullfile(projectRoot, 'tests'));
    addpath(fullfile(projectRoot, 'data'));

    fprintf('BB Codes project paths added.\n');
    fprintf('Project root: %s\n', projectRoot);
end