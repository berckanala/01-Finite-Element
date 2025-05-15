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


//+

Line(1)={1,2};
//+
Line(2)={2,3};
//+
Line(3)={3,4};
//+
Line(4)={4,5};
//+
Line(5)={5,6};
//+
Line(6)={6,7};
//+
Line(7)={7,8};
//+
Line(8)={8,9};
//+
Line(9)={9,10};
//+
Line(10)={10,11};
//+
Line(11)={11,12};
//+
Line(12)={12,13};
//+
Line(13)={13,14};
//+
Line(14)={14,1};
//+
Line(15)= {10,13};
Line(16)= {14,9};
Line(17)= {2,7};
Line(18)= {3,6};
Line(19)={17,18};





Curve Loop(10) = { 1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14};
//+
Plane Surface(20) = {10};

Physical Surface("Steel",201)={20};
//+
Physical Curve("BC_1", 123) = {11,4};
//+
Physical Point(202) = {15};

Physical Point(203) = {16};

//+
Transfinite Curve {1,17,7,8,16,14} = 6 Using Progression 1;
//+
Transfinite Curve {10, 11, 12, 15} = 10 Using Progression 1;
