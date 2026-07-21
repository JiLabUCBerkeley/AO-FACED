%set overwrite to 1 if you should overwrite no matter what.
function saveFigs_mod(saveDir,overwrite, filestr, namestr, reverseOrder)


if nargin < 5 reverseOrder = 0; end

if nargin <2 overwrite = 0; end
if nargin < 3
    if isunix filestr = '.fig';
    elseif ispc filestr = '.emf';
    else error('unsupported computer'); end
end

if isempty(strfind(filestr,'.')) error('wrong filestring'); end

%if ~overwrite overwrite = (menu('warning!!! overwrite files?','yes','no')==1); end

fs = get(0,'Children');
if reverseOrder fs = flipud(fs); end

if nargin < 1
    saveDir = '/autofs/space/flanders_003/users/vjsriniv/OCT_data/newdata/';
end

i = 0;
for k=1:1:length(fs)
    while (1)
        i = i + 1;

        istr = num2str(i);
        zplength = 4-length(istr);
        if zplength < 0 error('index too high'); end
        istr = [repmat('0',1,zplength),istr];     

        filename = [namestr '_' istr filestr]
        if ~exist([saveDir filename]) || overwrite
          if strcmp(filestr,'.pdf') figure(fs(k)); orient landscape; end            
          saveas(fs(k),[saveDir filename]); break;
        end     

     end
end