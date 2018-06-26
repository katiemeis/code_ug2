% Read lines into a cell array and make a hash from it
fid = fopen('gene_hash.txt');
C = textscan(fid,'%s %s');

key={C{2}{:}};
value=str2num({C{1}{:}});
D=containers.Map(key,value);

% Read the adjacency file line by line and load the value to a zero matrix based on the hash
M=zeros(length(D.keys));

fid2=fopen('GSE27715_1e5_0_cyto_named.adj');
line=fgetl(fid2);

while ischar(line)
    linecell=textscan(line,'%s');
    M(str2num(D(linecell{1}{1})),str2num(D(linecell{1}{2})))=str2num(linecell{1}{3});
    line = fgetl(fid2);
end
dlmwrite('t.txt',M);
fclose(fid);
fclose(fid2);


%{Test the cell array and the hash
celldisp(C)
D.size
D.keys
D.values
D('Ascl5')

if D.isKey('Ascl5')
    disp(D('Ascl5'))
else
    disp('The key dose not exist')
end
%}
