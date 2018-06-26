function nodeslist=nodesinshortestpath(genesofinterest,sparsematrix)
%find shortest paths between all gene pairs and put all nodes involved in the shortest paths into a vector
%input is a list of gene indcies from a gene hash table and a sparse matrix in which 0s represent no significant interactions and 1s represent significant interactions between two genes
%output is a vector of all nodes involved in the shortest paths betweeneach gene pair
path=[];
for i=genesofinterest
    for j=genesofinterest
        [d,p]=graphshortestpath(sparsematrix,i,j);
        path=[path p];
    end
end
nodeslist=unique(path);
end

