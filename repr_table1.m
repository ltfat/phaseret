function repr_table1
% The database is available from:
% http://data.cstr.ed.ac.uk/mocha/
% http://data.cstr.ed.ac.uk/mocha/unchecked/
% Just download and unpack the tarballs in the following directory: 
relDatabasepath = 'Databases';

basepath = fileparts(which(mfilename));
databasePath = [basepath,filesep,relDatabasepath];
subdirs = {'fsew0','msak0'};

exportdir = [basepath,filesep,'texexport'];

if ~exist(exportdir,'dir')
    mkdir(exportdir);
end

a = 128;
M = 1024;
fprintf('------------------ GAUSS WINDOW -----------------------\n');
comparetospsi(databasePath,subdirs,'exportdir',exportdir,'gauss','expname','mocha','storewavs','M',M,'a',a);
fprintf('------------------ TRUNC. GAUSS WINDOW -----------------------\n');
comparetospsi(databasePath,subdirs,'exportdir',exportdir,'truncgauss','expname','mocha','storewavs','M',M,'a',a);
fprintf('------------------ HANN WINDOW -----------------------\n');
comparetospsi(databasePath,subdirs,'exportdir',exportdir,'hann','expname','mocha','storewavs','M',M,'a',a);
fprintf('------------------ HAMMING WINDOW -----------------------\n');
comparetospsi(databasePath,subdirs,'exportdir',exportdir,'hamming','expname','mocha','storewavs','M',M,'a',a);
