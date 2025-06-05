SetFactory("OpenCASCADE");

// Parámetros de discretización
n_fino  = 10;
n_medio = 10;
r       = 1.3;  // Razón de progresión geométrica

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
Line(1)  = {1,2};  // inferior
Line(2)  = {2,3};  // derecha
Line(3)  = {3,4};  // superior
Line(4)  = {4,1};  // izquierda

// ============================
// Curve Loops y Superficie
// ============================
Curve Loop(105) = {4,3,2,1};
Plane Surface(105) = {105};

// Aplicar transfinite con progresión geométrica
Transfinite Curve {4} = n_fino Using Progression 1/r;
Transfinite Curve {2} = n_fino Using Progression r;  // Inversa para que coincida bien con la malla

Transfinite Curve {3} = n_fino Using Progression 1/r;
Transfinite Curve {1} = n_fino Using Progression r;  // Inversa para que coincida bien con la malla

// Superficie estructurada
Transfinite Surface {105};

