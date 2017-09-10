function [T,dc] = SpacePCA(SpecA,SpecB)

PH = Polyhedron(SpecA, SpecB);
Vertices = PH.V; %complexiy of vertex enumeration    which algo is in MPT
nVertices = size(Vertices,1);
dc = mean(Vertices);    
Vertices = Vertices-ones(nVertices,1)*dc;
Covariance = cov(Vertices);
[~,~,T] = svd(Covariance);
T = T'; 

    