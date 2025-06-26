SetFactory("OpenCASCADE");

//--------------------------------------
// PARÁMETROS
xApoyoIzq = 50;       // Coordenada X del apoyo inferior izquierdo
xApoyoDer = 170;      // Coordenada X del apoyo inferior derecho
yApoyoInf = 0;        // Coordenada Y base para los apoyos inferiores
Ancho_apoyos = 25;    // Altura de extrusión de los apoyos
AAA = Ancho_apoyos/2; // Mitad del ancho

//--------------------------------------
// PUNTOS DEL CONTORNO PRINCIPAL
Point(1) = { 0, 0, 0, 1.0};
Point(2) = {220, 0, 0, 1.0};
Point(3) = {220, 4, 0, 1.0};
Point(4) = { 0, 4, 0, 1.0};
Point(5) = { 0, 26, 0, 1.0};
Point(6) = {220, 26, 0, 1.0};
Point(7) = {220, 30, 0, 1.0};
Point(8) = { 0, 30, 0, 1.0};

//--------------------------------------
// LÍNEAS DEL CONTORNO
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};
Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 5};

//--------------------------------------
// APOYOS PARAMÉTRICOS
Point(20) = {110, 30, 0, 1.0};                // Apoyo superior
Point(30) = {xApoyoIzq, yApoyoInf, 0, 1.0};   // Apoyo inferior izquierdo
Point(40) = {xApoyoDer, yApoyoInf, 0, 1.0};   // Apoyo inferior derecho

// Apoyo superior circular
Point(21) = {100, 40, 0, 1.0};
Point(22) = {120, 40, 0, 1.0};
Point(23) = {110, 40, 0, 1.0};
Circle(401) = {21, 23, 20};
Circle(402) = {20, 23, 22};
Line(403) = {21, 22};
Line Loop(501) = {401, 402, -403};
Plane Surface(601) = {501};

// Apoyo inferior izquierdo circular
Point(31) = {xApoyoIzq - 10, yApoyoInf - 10, 0, 1.0};
Point(32) = {xApoyoIzq + 10, yApoyoInf - 10, 0, 1.0};
Point(33) = {xApoyoIzq, yApoyoInf - 10, 0, 1.0};
Circle(411) = {31, 33, 30};
Circle(412) = {30, 33, 32};
Line(413) = {31, 32};
Line Loop(511) = {411, 412, -413};
Plane Surface(611) = {511};

// Apoyo inferior derecho circular
Point(41) = {xApoyoDer - 10, yApoyoInf - 10, 0, 1.0};
Point(42) = {xApoyoDer + 10, yApoyoInf - 10, 0, 1.0};
Point(43) = {xApoyoDer, yApoyoInf - 10, 0, 1.0};
Circle(421) = {41, 43, 40};
Circle(422) = {40, 43, 42};
Line(423) = {41, 42};
Line Loop(521) = {421, 422, -423};
Plane Surface(621) = {521};

//--------------------------------------
// SUPERFICIES PRINCIPALES
Curve Loop(701) = {1, 2, 3, 4};
Plane Surface(7010) = {701};
Curve Loop(702) = {3, 5, 6, 7};
Plane Surface(7020) = {702};
Curve Loop(703) = {6, 8, 9, 10};
Plane Surface(7030) = {703};

//--------------------------------------
// Physical Surfaces
Physical Surface("Apoyo superior", 7032) = {601};
Physical Surface("Apoyo inferior izq", 7033) = {611};
Physical Surface("Apoyo inferior der", 7034) = {621};
Physical Surface("Ala superior", 7035) = {7030};
Physical Surface("Ala inferior", 7036) = {7010};
Physical Surface("Alma", 7037) = {7020};

//--------------------------------------
// EXTRUSIONES SIMPLIFICADAS — UNA SOLA POR VOLUMEN
Translate {0, 0, -10} { Surface{7010}; }
Extrude {0, 0, 20} { Surface{7010}; }

Translate {0, 0, -10} { Surface{7030}; }
Extrude {0, 0, 20} { Surface{7030}; }

Translate {0, 0, -2} { Surface{7020}; }
Extrude {0, 0, 4} { Surface{7020}; }

Translate {0, 0, -AAA} { Surface{601}; }
Extrude {0, 0, Ancho_apoyos} { Surface{601}; }

Translate {0, 0, -AAA} { Surface{611}; }
Extrude {0, 0, Ancho_apoyos} { Surface{611}; }

Translate {0, 0, -AAA} { Surface{621}; }
Extrude {0, 0, Ancho_apoyos} { Surface{621}; }//+
