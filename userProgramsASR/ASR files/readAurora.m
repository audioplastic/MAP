function x=readAurora(n)
fid=fopen(n,'r','b'); % aurora files are 16 bit big-endian
if fid
  [x,count]=fread(fid,inf,'int16');
  fclose(fid);
else
  warndlg(['Cannot load Aurora binary signal file ' n])
end  