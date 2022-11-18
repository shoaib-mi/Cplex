import cplex

let A = @[@[12, -51, 4],@[-51, 167, -68],@[4, -68, -41]].matrix
let C = @[1,2,3].row_vec
let B = @[@[12, 0, 0],@[0, 167, 0],@[0, 0, -41]].matrix
#let q = Evecs(A,decimal_tolerance=10,steps=1000)
#let lan = Evals(A)
#echo lan
#echo round(q.transpose*A*q,3)
echo round(diagonalize(A),3)
#var vecs = q.transpose
#echo vecs[0]
#echo A*vecs[0].transpose
#echo lan[0]*vecs[0].transpose

#echo vecs[0]*vecs[0].transpose