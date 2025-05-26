SetFactory("OpenCASCADE");

// Parámetros opcionales (no usados directamente en Gmsh)
n_fino = 15;
n_medio = 10;
n_latt=40;
n_latb=10;


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

Point(18) = {1000,200,0,1};
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

Line(23) = {6,21};
Line(24) = {21,3};
Line(28) = {1,21};
Line(30) = {21,18};
Line(31) = {8,1};

// ============================
// Curve Loops (cerrados correctamente)
// ============================

Curve Loop(105) = {1, 2, 24, -28};
Plane Surface(105) = {105};

Curve Loop(106) = {-24, 3, 4, 30};
Plane Surface(106) = {106};

Curve Loop(107) = {23, 30, 5, 6};
Plane Surface(107) = {107};

Curve Loop(108) = {31, -28, -23, 7};
Plane Surface(108) = {108};

// ============================
// Transfinite curves
// ============================

// Asigna nodos consistentes en todas las curvas transfinita
Transfinite Curve {2,28,7} = n_fino;  // Refinado para las curvas exteriores

Transfinite Curve {3,30,6} = n_medio;  // Refinado para la parte inferior
Transfinite Curve {31,23,5} = n_latt;  // Refinado para la parte lateral
Transfinite Curve {1,24,4} = n_latb;  // Refinado para la parte lateral


// ============================
// Transfinite surfaces y recombinación
// ============================

Transfinite Surface {105, 106, 107, 108};
Recombine Surface {105, 106, 107, 108};

// =============================
// PHYSICAL GROUPS
// =============================

// Grupo físico de todas las superficies de acero
Physical Surface("Steel1", 201) = {105, 108};
Physical Surface("Steel2", 202) = {107,106};

// Grupo físico de curvas de borde para condiciones (por ejemplo, carga o apoyo)
Physical Curve("BC_1", 123) = {4, 5};  // Ajustado: solo curvas válidas existentes

// Apoyos puntuales
Physical Curve("BC_R1", 203) = {31};  // Nodo superior izquierdo
