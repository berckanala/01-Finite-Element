SetFactory("OpenCASCADE");

// Parámetros de discretización
n_fino  = 10;
n_medio = 10;
r       = 1.5;  // Razón de progresión geométrica

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
Line(10)  = {1,2};  // inferior
Line(20)  = {2,3};  // derecha
Line(30)  = {3,4};  // superior
Line(40)  = {4,1};  // izquierda

// ============================
// Curve Loops y Superficie
// ============================
Curve Loop(105) = {40,30,20,10};
Plane Surface(105) = {105};

// Aplicar transfinite con progresión geométrica
Transfinite Curve {40} = n_fino Using Progression 1/r;
Transfinite Curve {20} = n_fino Using Progression r;  // Inversa para que coincida bien con la malla

Transfinite Curve {30} = n_fino Using Progression 1/r;
Transfinite Curve {10} = n_fino Using Progression r;  // Inversa para que coincida bien con la malla

// Superficie estructurada
Transfinite Surface {105};

Physical Surface("Omega") = {105};

Physical Line("1") = {10};
Physical Line("2") = {20};
Physical Line("3") = {30};
Physical Line("4") = {40};
