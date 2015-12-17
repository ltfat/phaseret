function repr_table2
% The EBU SQAM database is available freely from 
% https://tech.ebu.ch/publications/sqamcd
% Just download and unpack the archive in the following directory 
relDatabasepath = 'Databases';
% such that all the files reside in the following subdirectory
subdirs = {'SQAM'};
% The resulting files are flac files.
% The following Bash command converts all flac files in the
% current directory to wav files provided "flac - Command-line FLAC
% encoder/decoder" is installed
% ls *.flac | xargs flac -df


basepath = fileparts(which(mfilename));
databasePath = [basepath,filesep,relDatabasepath];


exportdir = [basepath,filesep,'texexport'];

if ~exist(exportdir,'dir')
    mkdir(exportdir);
end

a = 256;
M = 2048;
maxsamples = 10*44100;
fprintf('------------------ GAUSS WINDOW -----------------------\n');
comparetospsi(databasePath,subdirs,'exportdir',exportdir,'gauss','expname','sqam','M',M,'a',a,'storewavs','maxsamples',maxsamples);
fprintf('------------------ TRUNC. GAUSS WINDOW -----------------------\n');
comparetospsi(databasePath,subdirs,'exportdir',exportdir,'truncgauss','expname','sqam','M',M,'a',a,'storewavs','maxsamples',maxsamples);
fprintf('------------------ HANN WINDOW -----------------------\n');
comparetospsi(databasePath,subdirs,'exportdir',exportdir,'hann','expname','sqam','M',M,'a',a,'storewavs','maxsamples',maxsamples);
fprintf('------------------ HAMMING WINDOW -----------------------\n');
comparetospsi(databasePath,subdirs,'exportdir',exportdir,'hamming','expname','sqam','M',M,'a',a,'storewavs','maxsamples',maxsamples);