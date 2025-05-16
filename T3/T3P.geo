SetFactory("OpenCASCADE");

// =====================
// Parámetros opcionales
// =====================
n_fino = 60;
n_medio = 30;
n_grueso = 10;

// =================
// Definición de puntos
// =================
Point(1)  = {0,200,0,1};
Point(2)  = {200,0,0,1};
Point(3)  = {600,0,0,1};
Point(4)  = {1000,0,0,1};
Point(5)  = {1000 ,1000,0,1};
Point(6)  = {600 ,1000,0,1};
Point(7)  = {200 ,1000,0,1};
Point(8)  = {0,1000,0,1};
Point(9)  = {-200 ,1000,0,1};
Point(10) = {-600 ,1000,0,1};
Point(11) = {-1000 ,1000,0,1};
Point(12) = {-1000,0,0,1};
Point(13) = {-600,0,0,1};
Point(14) = {-200,0,0,1};

Point(15) = {0,800,0,1};
Point(16) = {0,700,0,1};

Point(17) = {-1000,200,0,1};
Point(18) = {1000,200,0,1};

Point(19) = {-600,200,0,1};
Point(20) = {-200,200,0,1};
Point(21) = {600,200,0,1};
Point(22) = {200,200,0,1};

// =============
// Definición de líneas
// =============
Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,18};
Line(5)  = {18,5};
Line(6)  = {5,6};
Line(7)  = {6,7};
Line(8)  = {7,8};
Line(9)  = {8,9};
Line(10) = {9,10};
Line(11) = {10,11};
Line(12) = {11,17};
Line(13) = {17,12};
Line(14) = {12,13};
Line(15) = {13,14};
Line(16) = {14,1};

Line(17) = {10,19};
Line(18) = {19,13};
Line(19) = {9,20};
Line(20) = {20,14};
Line(21) = {7,22};
Line(22) = {22,2};
Line(23) = {6,21};
Line(24) = {21,3};
Line(25) = {17,19};
Line(26) = {19,20};
Line(27) = {20,1};
Line(28) = {1,22};
Line(29) = {22,21};
Line(30) = {21,18};
Line(31) = {8,1};

// ============================
// Mallado transfinito de curvas
// ============================
Transfinite Curve {1, 22, 28, 21, 8, 31, 27, 16, 20, 9, 19} = n_fino;
Transfinite Curve {29, 23, 7, 24, 2, 26, 10, 17, 18, 15} = n_medio;
Transfinite Curve {25, 11, 12, 13, 14, 30, 6, 5, 3, 4} = n_grueso;

// ===================
// Superficie principal
// ===================
Curve Loop(10) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
Plane Surface(20) = {10};

// ===================
// Mallado estructurado
// ===================


// ===================
// Physical groups (opcional)
// ===================
Physical Surface("Steel", 201) = {20};
Physical Curve("BC_1", 123) = {12,13,5,4};
Physical Point("Load", 202) = {15};
Physical Point("Sensor", 203) = {16};

// ================
// Opcional: Quad9
// ================
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
