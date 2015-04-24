// Gmsh project created on Thu Apr 23 16:02:20 2015
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, -0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {1, 1, 0, 1.0};
Point(5) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 5};
Line(4) = {5, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Physical Surface(7) = {6};
Physical Point(8) = {1};
Physical Point(9) = {2};
Physical Point(10) = {3};
Physical Point(11) = {5};
