function [subM,bg]=subsetmatrix(nodelist,M,revhash)
%This function subset a matrix and plot the submatrix
%M is the original matrix 
%nodelist is a list of elements from matrix M, i.e. all genes involved in the shortest path
%revhash is a reverse hash of gene names in which keys are index and values are gene names

subM=zeros(length(nodelist));
for i=[1:length(nodelist)]
    for j=[1:length(nodelist)]
        subM(i,j)=M(nodelist(i),nodelist(j));
    end
end

%Make Node ID cell array for graph visualization 
genename_nodelist={};
for i=1:length(nodelist)
    genename_nodelist{i}=revhash(num2str(nodelist(i)));
end
%Plot subM
bg = biograph(subM,genename_nodelist,'ShowArrows','off','ShowWeights','off');
%use view(bg) to show the graph)
end

