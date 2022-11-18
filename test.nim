import math,cplex

var 
    x:Complex = 2.1.complex
    y:Complex = 1.0.complex
    z:Complex = [1.0,2.0]

echo x+y+z*conj(z)
echo x-y
echo y-x
echo z.type
echo z.complex
echo z/(-2.5)
echo abs(z)^2
echo conj(z)
echo z*conj(z)
echo z + conj(z)
echo complex(2,8)
echo z
echo z.complex

echo round(exp(2+I*Pi))
echo exp(2)*exp(I*Pi)
echo round(0.23333,1)
echo sin(Pi/4*I)
echo sinh(Pi/4*I)

#[ 
var 
    x:Complex
    y: seq[seq[Complex]]

proc determinant(matrix:seq[seq[Complex]]):seq[seq[Complex]] =
    result = matrix

y = @[@[1.0 + I, 2.0 + I],@[3.0 + I, 4.0 + I]]
var z = determinant(y)
echo z
]#