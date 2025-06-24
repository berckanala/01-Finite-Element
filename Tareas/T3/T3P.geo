SetFactory("OpenCASCADE");

// Parámetros opcionales (no usados directamente en Gmsh)
n_fino = 20;
n_medio = 30;
n_grueso = 20;

// =================
// Definición de puntos
// =================
Point(1)  = {0,200,0,1};
Point(2)  = {200,0,0,1};
Point(3)  = {600,0,0,1};
Point(4)  = {1000,0,0,1};
Point(5)  = {1000 ,1000,0,1};
Point(6)  = {600 ,1000,0,1};
Point(8)  = {0,1000,0,1};
Point(10) = {-600 ,1000,0,1};
Point(11) = {-1000 ,1000,0,1};
Point(12) = {-1000,0,0,1};
Point(13) = {-600,0,0,1};
Point(14) = {-200,0,0,1};


Point(17) = {-1000,200,0,1};
Point(18) = {1000,200,0,1};

Point(19) = {-600,200,0,1};
Point(21) = {600,200,0,1};

// =================
// Definición de líneas
// =================
Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,18};
Line(5)  = {18,5};
Line(6)  = {5,6};
Line(7)  = {6,8};
Line(9)  = {8,10};
Line(11) = {10,11};
Line(12) = {11,17};
Line(13) = {17,12};
Line(14) = {12,13};
Line(15) = {13,14};
Line(16) = {14,1};

Line(17) = {10,19};
Line(18) = {19,13};
Line(23) = {6,21};
Line(24) = {21,3};
Line(25) = {17,19};
Line(26) = {19,1};
Line(28) = {1,21};
Line(30) = {21,18};
Line(31) = {8,1};

// ============================
// Curve Loops (cerrados correctamente)
// ============================

Curve Loop(101) = {31, 9, 17, -26};
Plane Surface(101) = {101};

Curve Loop(1002) = {-18, 15, 16, -26};
Plane Surface(102) = {1002};

Curve Loop(103) = {11, 12, 25, -17};
Plane Surface(103) = {103};

Curve Loop(104) = {13, 14, -18, -25};
Plane Surface(104) = {104};

Curve Loop(105) = {1, 2, 24, -28};
Plane Surface(105) = {105};

Curve Loop(106) = {-24, 3, 4, 30};
Plane Surface(106) = {106};

Curve Loop(107) = {23, 30, 5, 6};
Plane Surface(107) = {107};

Curve Loop(108) = {31, -28, -23, 7};
Plane Surface(108) = {108};

// ============================
// Transfinite curves (opcional)
// ============================

Transfinite Curve {1,2,3,4,5,6,7,9,11,12,13,14,15,16,
                  17,18,23,24,25,26,28,30,31} = 5;

// Transfinite surface por bloque
Transfinite Surface {101,102,103,104,105,106,107,108};
Recombine Surface {101,102,103,104,105,106,107,108};
// =============================
// PHYSICAL GROUPS
// =============================

// Todas las superficies
Physical Surface("Steel",201)={101,102,103,104,105,106,107,108};
//+
Physical Curve("BC_1", 123) = {12,13,5,4};
//+
Physical Point("BC_R1", 202) = {8};
Physical Point("BC_R2", 203) = {1};