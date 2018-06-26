%Graphs everything on its own axis. Produces 4 separate figures


%read lines into cell array and map name to probe ID
fid = fopen('GeneListPlotting.txt', 'r');
probes = textscan(fid, '%s%s%s%s%s%s%s', 'Delimiter', '\t');

key = {probes{1,3}{:,1}};
value = {probes{1,1}{:,1}};
name2probe = containers.Map(key, value); 

%read adjacency file
fid2 = fopen('GSE49991-matrix-matlab.txt', 'r');
fs = ['%q', repmat('%f', 1, 177)];
data = textscan(fid2, fs, 'Delimiter', ',');

%map probe name to row number/index
key = data{1,1};
index = [1:length(data{1,1})];
value = index(:);
probe2row = containers.Map(key, value);

%find the row we want
probeID1 = name2probe('Zfp7');
rownum1 = probe2row(probeID1);

probeID2 = name2probe('Aqp7');
rownum2 = probe2row(probeID2);

dmatrix = {data(1,2:end)};
C = dmatrix{1,1};
datamatrix = cell2mat(C);
plotvect1 = datamatrix(rownum1, :);
plotvect2 = datamatrix(rownum2, :);

%average triplicates
newplotvect1 = reshape(plotvect1, 3, []);
avgpoints1 = sum(newplotvect1, 1)./size(newplotvect1, 1);

newplotvect2 = reshape(plotvect2, 3, []);
avgpoints2 = sum(newplotvect2, 1)./size(newplotvect2, 1);

%plot erythroid conditions gene1
figure(1)
x = [0:2:24 27:3:48 52:4:72 96 120 168];
y = avgpoints1(1, 1:30);
plot(x,y)
xlabel('Time')
ylabel('Expression')
title('Erythroid for gene x')

%plot neutrophil conditions gene1
figure(2)
x = [2:2:24 27:3:48 52:4:72 96 120 168];
y = avgpoints1(1, 31:end);
plot(x,y)
xlabel('Time')
ylabel('Expression')
title('Neutrophil for gene x')

%plot erythroid conditions gene2
figure(3)
x = [0:2:24 27:3:48 52:4:72 96 120 168];
y = avgpoints2(1, 1:30);
plot(x,y)
xlabel('Time')
ylabel('Expression')
title('Erythroid for gene y')

%plot neutrophil conditions gene2
figure(4)
x = [2:2:24 27:3:48 52:4:72 96 120 168];
y = avgpoints2(1, 31:end);
plot(x,y)
xlabel('Time')
ylabel('Expression')
title('Neutrophil for gene y')

