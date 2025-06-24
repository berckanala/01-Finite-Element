SetFactory("OpenCASCADE");

n_malla=5;

// =================
// Definición de puntos
// =================
Point(1)  = {0,0,0,1};
Point(2)  = {1,0,0,1};
Point(3)  = {1,1,0,1};
Point(4)  = {0,1,0,1};


// =================
// Definición de líneas
// =================
Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,1};



Curve Loop(105) = {1, 2, 3, 4};
Plane Surface(105) = {105};


Circle(5) = {1, 1, 0, 0.707, 0, 2*Pi};
//+
Curve Loop(106) = {5};
//+
Surface(106) = {106};
//+
BooleanDifference{ Surface{105}; Delete; }{ Surface{106}; Delete; }
//+


Transfinite Curve {1,5,4} = n_malla;
Transfinite Curve {3} = 4;
Transfinite Curve {2} = 2;
Transfinite Surface {105} = {2, 5, 4, 1};
Recombine Surface {105};  // Opcional si quieres elementos cuadriláteros

Physical Line("BC_R1") = {1,2,3,4,5};
//+
Physical Surface("Omega",107) = {105};
