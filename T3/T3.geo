SetFactory("OpenCASCADE");

Point(1)={0,200,0,1};
Point(2)={200,0,0,1};
Point(3)={600,0,0,1};
Point(4)={1000,0,0,1};
Point(5)={1000 ,1000,0,1};
Point(6)={600 ,1000,0,1};
Point(7)={200 ,1000,0,1};
Point(8)={0,1000,0,1};
Point(9)={-200 ,1000,0,1};
Point(10)={-600 ,1000,0,1};
Point(11)={-1000 ,1000,0,1};
Point(14)={-200,0,0,1};
Point(13)={-600,0,0,1};
Point(12)={-1000,0,0,1};

Point(17)={-1000,200,0,1};
Point(18)={1000,200,0,1};

Point(15)={0,800,0,1};
Point(16)={0,700,0,1};

Point(19)={-600,200,0,1};
Point(20)={-200,200,0,1};
Point(21)={600,200,0,1};
Point(22)={200,200,0,1};

//+

Line(1)={1,2};
//+
Line(2)={2,3};
//+
Line(3)={3,4};
//+
Line(4)={4,18};
//+
Line(5)={18,5};
//+
Line(6)={5,6};
//+
Line(7)={6,7};
//+
Line(8)={7,8};
//+
Line(9)={8,9};
//+
Line(10)={9,10};
//+
Line(11)={10,11};
//+
Line(12)={11,17};
//+
Line(13)={17,12};
//+
Line(14)={12,13};
//+
Line(15)={13,14};
//+
Line(16)={14,1};
//+
Line(17)= {10,19};
Line(18)= {19,13};
Line(19)= {9,20};
Line(20)= {20,14};
Line(21)={7,22};
Line(22)={22,2};
Line(23)={6,21};
Line(24)={21,3};
Line(25)={17,19};
Line(26)={19,20};
Line(27)={20,1};
Line(28)={1,22};
Line(29)={22,21};
Line(30)={21,18};
Line(31)={8,1};






Curve Loop(10) = { 1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16};
//+
Plane Surface(20) = {10};

Physical Surface("Steel",201)={20};
//+
Physical Curve("BC_1", 123) = {12,13,5,4};
//+
Physical Point(202) = {15};

Physical Point(203) = {16};
//+
Transfinite Curve {1,22,28} = 60;
//+
Transfinite Curve {28,21,8,31} = 60;
//+
Transfinite Curve {27,16,20} = 60 ;
//+
Transfinite Curve {27,31,9,19} = 50 ;


//+
Transfinite Curve {29,23,7,21} = 30 ;
//+
Transfinite Curve {29,24,2,22} = 30 ;
//+
Transfinite Curve {19,10,17,26} = 30;
//+
Transfinite Curve {26,20,15,18} = 30;

//+
Transfinite Curve {17,11,12,25} = 10 ;
//+
Transfinite Curve {25,18,14,13} = 10 ;
//+
Transfinite Curve {6,23,30,5} = 10 ;
//+
Transfinite Curve {30,4,3,24} = 10 ;
