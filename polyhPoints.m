function PH = polyhPoints(n)

if nargin == 0
    n = 7;
end

[lat, lon] = ginput(n);
PH = Polyhedron('V',[lat lon]);
plot(PH)
    