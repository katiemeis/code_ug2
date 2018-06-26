function [hashtable, revhash, M] = readAdjacencyMatrix(hash_file, adj_file, edgeweights)
% [hashtable, M] = readAdjacencyMatrix(hash_file, adj_file)
%
% Reads in a CYTOSCAPE format
% adjacency matrix as a numeric matrix in MATLAB. Uses the hash_file to
% build a hash table to index the matrix. Returns the hashtable and the
% matrix M.

% Check whether correct number of arguments have been specified
error(nargchk(2, 3, nargin));

% give default value to edgeweights 
if nargin == 2
    edgeweights = 0;
end

% check correct value of edgeweights
if (edgeweights ~= 0) && (edgeweights ~= 1)
    error('Incorrect value %d for edgeweights\n', edgeweights);
end    

% Read lines into a cell array and make a hash from it
fid = fopen(hash_file);
C = textscan(fid,'%s %s');

key={C{2}{:}};
value={C{1}{:}};
hashtable=containers.Map(key,value);

%Make a reverse hash
revkey=values(hashtable);
revvalue=keys(hashtable);
revhash=containers.Map(revkey,revvalue);

% Read the adjacency file line by line and load the value to a zero matrix based on the hash
M=zeros(length(hashtable.keys));

fid2=fopen(adj_file);
line=fgetl(fid2);

while ischar(line)

    linecell=textscan(line,'%s');

    if edgeweights
        M(str2num(hashtable(linecell{1}{1})),str2num(hashtable(linecell{1}{2})))=str2num(linecell{1}{3});
    else    
        M(str2num(hashtable(linecell{1}{1})),str2num(hashtable(linecell{1}{2})))=1;
    end

    line = fgetl(fid2);
end
fclose(fid);
fclose(fid2);
end
