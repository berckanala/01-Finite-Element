// The essential boundary length corresponds to the entire surface where the nut can make contact. 
// Since the nut is approximately 10 mm long on each side, the effective contact length is around 17 mm on both sides. 
// However, the nut can be positioned anywhere along the full length of the inner surface of the wrench head.

SetFactory("OpenCASCADE");


Circle(1) = {0, 0, 0, 29/2, 0, 2*Pi};
Curve Loop(1) = {1};
Plane Surface(1) = {1};


Circle(2) = {6, 0, 0, 9, 0, 2*Pi};
Curve Loop(2) = {2};
Plane Surface(2) = {2};

Rectangle(3) = {-4, -13/2, 0, 20, 13, 0};

BooleanFragments{ Surface{2}; Surface{3}; Delete; }{ }
//+
Recursive Delete {
  Surface{3}; Surface{5}; Surface{4}; 
}
//+
//+
BooleanFragments{ Surface{2}; Surface{1}; Surface{6}; Delete; }{ }
//+
Recursive Delete {
  Surface{1}; Surface{3}; Surface{2}; Surface{6}; Surface{4}; 
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, 17*Pi/180} {
  Surface{5}; 
}
//+
Rectangle(6) = {-110/2, -11/2, 0, 110, 11, 0};
//+
Translate {55/2+26.75, 0, 0} {
  Surface{5}; 
}

Circle(100) = {0, 0, 0, 32/2, 0, 2*Pi};
Curve Loop(100) = {100};
Plane Surface(100) = {100};


Circle(101) = {6, 0, 0, 13, 0, 2*Pi};
Curve Loop(101) = {101};
Plane Surface(101) = {101};

Rectangle(103) = {-6, -17/2, 0, 25, 17, 0};

//+
BooleanFragments{ Surface{100}; Surface{101}; Surface{103}; Delete; }{ }
//+
Recursive Delete {
  Surface{7}; Surface{14}; Surface{15}; Surface{18}; Surface{17}; Surface{16}; 
}
//+
Recursive Delete {
  Surface{13}; 
}
//+
BooleanUnion{ Surface{8}; Delete; }{ Surface{10}; Surface{11}; Surface{12}; Surface{9}; Delete; }

Rotate {{0, 0, 1}, {0, 0, 0}, 197*Pi/180} {
  Surface{8}; 
}

Translate {-55/2-28.25, 0, 0} {
  Surface{8}; 
}//+
BooleanFragments{ Surface{8}; Surface{6}; Surface{5}; Delete; }{ }
//+
Recursive Delete {
  Surface{7}; Surface{4}; 
}
//+
BooleanUnion{ Surface{2}; Delete; }{ Surface{3}; Surface{1}; Delete; }
//+
BooleanUnion{ Surface{8}; Delete; }{ Surface{6}; Surface{9}; Delete; }


x110 = (55/2+26.75) - Sqrt((23+14.5)^2 - (23+11/2)^2);
y110 = 23+11/2;


Circle(110) = {x110, y110, 0, 23, 0, 2*Pi};
Curve Loop(110) = {110};
Plane Surface(110) = {110};

Circle(111) = {x110, -y110, 0, 23, 0, 2*Pi};
Curve Loop(111) = {111};
Plane Surface(111) = {111};

Rectangle(112) = {-110/2, -11/2-11, 0, 110, 11, 0};
Rectangle(113) = {-110/2, 11/2, 0, 110, 11, 0};

x120 = -((55/2+28.25) - Sqrt((16+47)^2 - (47+11/2)^2));
y120 = 11/2 + 47;


Circle(120) = {x120, y120, 0, 47, 0, 2*Pi};
Curve Loop(120) = {120};
Plane Surface(120) = {120};

Circle(121) = {x120, -y120, 0, 47, 0, 2*Pi};
Curve Loop(121) = {121};
Plane Surface(121) = {121};//+
BooleanFragments{ Surface{113}; Surface{8}; Surface{110}; Surface{120}; Surface{2}; Surface{112}; Surface{121}; Surface{111}; Delete; }{ }
//+
Recursive Delete {
  Surface{23}; Surface{22}; Surface{21}; Surface{39}; Surface{38}; Surface{37}; Surface{29}; Surface{30}; Surface{33}; Surface{34}; Surface{28}; Surface{9}; Surface{11}; Surface{10}; Surface{13}; Surface{15}; 
}
//+
Recursive Delete {
  Surface{6}; Surface{32}; Surface{16}; Surface{36}; 
}
//+
BooleanUnion{ Surface{5}; Delete; }{ Surface{8}; Surface{31}; Surface{12}; Surface{35}; Delete; }
//+
BooleanUnion{ Surface{26}; Delete; }{ Surface{7}; Surface{27}; Surface{25}; Surface{24}; Delete; }
//+
BooleanUnion{ Surface{18}; Delete; }{ Surface{14}; Surface{17}; Surface{20}; Surface{19}; Delete; }

Circle(129) = {55/2, 0, 0, 3.25, 0, 2*Pi};
Curve Loop(129) = {129};
Plane Surface(129) = {129};
//+
Circle(130) = {-55/2, 0, 0, 3.25, 0, 2*Pi};
Curve Loop(130) = {130};
Plane Surface(130) = {130};
//+
Rectangle(24) = {-55/2, -6.5/2, 0, 55, 6.5, 0};
//+//+
BooleanFragments{ Surface{5}; Surface{24}; Surface{130}; Surface{129}; Delete; }{ }
//+
BooleanUnion{ Surface{28}; Delete; }{ Surface{30}; Surface{29}; Surface{31}; Surface{32}; Delete; }
//+
Physical Surface("Heads", 123) = {26, 18};
//+
Physical Surface("Handle_big", 124) = {27};
//+
Physical Surface("Handle_thin", 125) = {28};
//+
Coherence;

//+
Physical Curve("Force", 129) = {119};
//+
Physical Curve("Nut", 130) = {126, 123};


//+
Coherence;
