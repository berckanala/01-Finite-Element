SetFactory("OpenCASCADE");


//--------------------------------------
// PARÁMETROS
r = 9;                // Radio de los círculos del alma
xApoyoIzq = 20;       // Coordenada X del apoyo inferior izquierdo
xApoyoDer = 200;      // Coordenada X del apoyo inferior derecho
yApoyoInf = 0;        // Coordenada Y base para los apoyos inferiores
Ancho_apoyos = 10;    // Altura de extrusión de los apoyos
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
// CÍRCULOS DEL ALMA — CENTROS EN Y = 15
For i In {0:3}
  x = 45 + i * 45;
  id = 11 + i;
  pid = 101 + i*10;
  cid = 201 + i*10;
  lid = 301 + i*10;

  Point(id) = {x, 15, 0, 1.0};
  Point(pid+0) = {x + r, 15, 0, 1.0};
  Point(pid+1) = {x, 15 + r, 0, 1.0};
  Point(pid+2) = {x - r, 15, 0, 1.0};
  Point(pid+3) = {x, 15 - r, 0, 1.0};

  Circle(cid+0) = {pid+0, id, pid+1};
  Circle(cid+1) = {pid+1, id, pid+2};
  Circle(cid+2) = {pid+2, id, pid+3};
  Circle(cid+3) = {pid+3, id, pid+0};
  Line Loop(lid) = {cid+0, cid+1, cid+2, cid+3};
EndFor

//--------------------------------------
// APOYOS PARAMÉTRICOS
Point(20) = {110, 30, 0, 1.0};                // Apoyo superior
Point(30) = {xApoyoIzq, yApoyoInf, 0, 1.0};   // Apoyo inferior izquierdo
Point(40) = {xApoyoDer, yApoyoInf, 0, 1.0};   // Apoyo inferior derecho

// Apoyo superior
Point(21) = {110, 40, 0, 1.0};
Point(22) = {100, 40, 0, 1.0};
Point(23) = {120, 40, 0, 1.0};
Circle(401) = {22, 21, 20};
Circle(402) = {20, 21, 23};
Line(403) = {22, 21};
Line(404) = {21, 23};
Line(405) = {22, 23};
Line Loop(501) = {401, 402, -404, -403};
Plane Surface(601) = {501};

// Apoyo inferior izquierdo
Point(31) = {xApoyoIzq, yApoyoInf - 10, 0, 1.0};
Point(32) = {xApoyoIzq - 10, yApoyoInf - 10, 0, 1.0};
Point(33) = {xApoyoIzq + 10, yApoyoInf - 10, 0, 1.0};
Circle(411) = {32, 31, 30};
Circle(412) = {30, 31, 33};
Line(413) = {32, 31};
Line(414) = {31, 33};
Line(415) = {32, 33};
Line Loop(511) = {411, 412, -414, -413};
Plane Surface(611) = {511};

// Apoyo inferior derecho
Point(41) = {xApoyoDer, yApoyoInf - 10, 0, 1.0};
Point(42) = {xApoyoDer - 10, yApoyoInf - 10, 0, 1.0};
Point(43) = {xApoyoDer + 10, yApoyoInf - 10, 0, 1.0};
Circle(421) = {42, 41, 40};
Circle(422) = {40, 41, 43};
Line(423) = {42, 41};
Line(424) = {41, 43};
Line(425) = {42, 43};
Line Loop(521) = {421, 422, -424, -423};
Plane Surface(621) = {521};

//--------------------------------------
// SUPERFICIES PRINCIPALES
Curve Loop(701) = {1, 2, 3, 4};
Plane Surface(7010) = {701};

Curve Loop(702) = {3, 5, 6, 7};
Plane Surface(7020) = {702, 301, 311, 321, 331};

Curve Loop(703) = {6, 8, 9, 10};
Plane Surface(7030) = {703};

//--------------------------------------
// EXTRUSIONES SIMPLIFICADAS — UNA POR VOLUMEN
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
Extrude {0, 0, Ancho_apoyos} { Surface{621}; }

//--------------------------------------
// UNIFICAR LOS VOLUMENES DE LA VIGA
Compound Volume{1,2,3};

//--------------------------------------
// DEFINIR VOLUMENES FISICOS
// Ojo: después de Compound, la viga es solo un volumen (1)
Physical Volume("Viga", 7079) = {1};
Physical Volume("BC_1", 7080) = {4};
Physical Volume("BC_R1", 7081) = {5,6};
