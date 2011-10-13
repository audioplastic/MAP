function createMLFFiles(experimentName, filePath, withSP)

d = dir([filePath, '*.map']);

if withSP == 1
    fid = fopen([experimentName '_sp.mlf'],'w');
else
    fid = fopen([experimentName '.mlf'],'w');
end;

fprintf(fid,'#!MLF!#\n');

for I = 1:size(d,1)
    fprintf(fid,['"*/' d(I).name(1:end-3) 'lab"\n']);
    fprintf(fid,'sil\n');
    labels = d(I).name(5:end-5);
    for J = 1:length(labels);
        switch labels(J)
            case 'O'
                S = 'oh';
            case '1'
                S = 'one';
            case '2'
                S = 'two';
            case '3'
                S = 'three';
            case '4'
                S = 'four';
            case '5'
                S = 'five';
            case '6'
                S = 'six';
            case '7'
                S = 'seven';
            case '8'
                S = 'eight';
            case '9'
                S = 'nine';
            case 'Z'
                S = 'zero';
        end;
        fprintf(fid,S);
        if (J < length(labels)) && (withSP == 1)
            fprintf(fid,'\nsp');
        end;
        fprintf(fid,'\n');
    end;
    fprintf(fid,'sil\n.\n');
end;

fclose(fid);
