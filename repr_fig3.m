function repr_fig3(winslice,wsslice)
% The speech database is available from:
% http://data.cstr.ed.ac.uk/mocha/
% http://data.cstr.ed.ac.uk/mocha/unchecked/
% and the EBU SQAM database is available freely from 
% EBU SQAM
% https://tech.ebu.ch/publications/sqamcd
% Just download and unpack the tarballs in the following directory: 
relDatabasepath = 'Databases';
% In case of the EBU SQAM database, the resulting files are flac.
% The following Bash command converts all flac files in the
% current directory to wav files provided flac - Command-line FLAC encoder/decoder
% is installed
% ls *.flac | xargs flac -df

basepath = fileparts(which(mfilename));
databasePath = [basepath,filesep,relDatabasepath];
subdirs = {'fsew0','msak0'};

exportdir = [basepath,filesep,'texexport'];

if ~exist(exportdir,'dir')
    mkdir(exportdir);
end
maxwavno = [];


win = {'gauss','truncgauss','hann','hamming'};
ws = {'nowarmstart','wsheapint'};

if nargin>0
    win = win(winslice);
end
if nargin>1
    ws = ws(wsslice);
end

a = 128;
M = 1024;
maxit = 200;

for wsId = 1:numel(ws)
    for ii=1:numel(win)
        w = win{ii};
        s = ws{wsId};
        fprintf('------------------ %s WINDOW -----------------------\n',upper(w));
        comparetoall(databasePath,subdirs,'exportdir',exportdir,w,...
                     'storewavs','maxwavno',maxwavno,'maxit',maxit,s,...
                     'expname','mocha','M',M,'a',a);   
    end
end



