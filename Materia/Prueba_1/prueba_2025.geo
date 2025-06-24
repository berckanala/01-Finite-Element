// Gmsh project created on Thu May 08 12:23:32 2025
SetFactory("OpenCASCADE");

//+
Point(1)={0.0, 0.0, 0.0, 1.0};
//+
Point(2)={280.0, 60.0, 0.0, 1.0};
//+
Point(3)={280.0, 160.0, 0.0, 1.0};
//+
Point(4)={61.0, 160.0, 0.0, 1.0};
//+
Point(5)={60.5, 140.0, 0.0, 1.0};
//+
Point(6)={60.0, 160.0, 0.0, 1.0};
//+
Point(7)={0.0, 160.0, 0.0, 1.0};
//+

Line(1)={1,2};
//+
Line(2)={2,3};
//+
Line(3)={3,4};
//+
Line(4)={4,5};
//+
Line(5)={5,6};
//+
Line(6)={6,7};
//+
Line(7)={7,1};

//+
Coherence;
//+
Curve Loop(10) = { 1, 2, 3, 4, 5, 6,7};
//+
Coherence;
//+
Plane Surface(20) = {10};
//+
Physical Surface("Concrete",201)={20};
//+
Physical Curve("BC_1", 123) = {7};
//+
Physical Curve("BC_F", 125) = {6,3};
