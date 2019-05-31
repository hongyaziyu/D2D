a=[10,9,5;7,8,9;9,7,9];
 [assignmentSum, cost] = munkres(-a);
 [assignmentMin, dummyMin ] = maxMin( a );
 