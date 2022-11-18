import math

type
    Complex* = array[2,float]
    RowVec* = seq[Complex]
    ColVec* = seq[array[1,Complex]]
    Matrix* = seq[RowVec]

const
    I*: Complex = [0.0,1.0]

# ----------------------inside complex operation----------------------
proc complex*(a:int or float): Complex = [a.float,0.0]
proc complex*(a:Complex): Complex = a
proc complex*(a,b:int or float): Complex = [a.float,b.float]
proc complex*(a:array[2,int or float]): Complex = [a[0].float,a[1].float]

#--------------------inside RowVec operation-------------------
proc row_vec*(a: seq[int or float or Complex]): RowVec = 
    for i in 0..len(a)-1:
        result.add(a[i].complex)

proc col_vec*(a: seq[int or float or Complex]): ColVec = 
    for i in 0..len(a)-1:
        result.add([a[i].complex])

proc matrix*(a: seq[seq[int or float or Complex]]): Matrix = 
    for i in 0..len(a)-1:
        result.add(@[])
        for j in 0..len(a[i])-1:
            result[i].add(a[i][j].complex)

#-----------------normal arithmetic operations-----------------
proc abs*(a:Complex): float = sqrt(a[0]^2+a[1]^2)
proc conj*(a:Complex): Complex = [a[0],-a[1]]
proc round*(a:Complex,b:int=0):Complex = [round(a[0],b),round(a[1],b)]
proc `+`*(a:Complex): Complex = a
proc `+`*(a,b:Complex): Complex = [a[0]+b[0],a[1]+b[1]]
proc `+`*(a,b:int or float): float = a.float + b.float 
proc `+`*(a:int or float, b:Complex): Complex = a.complex + b
proc `+`*(a:Complex, b:int or float): Complex = a + b.complex
proc `-`*(a:Complex): Complex = [-a[0],-a[1]]
proc `-`*(a,b:Complex): Complex = [a[0]-b[0],a[1]-b[1]]
proc `-`*(a:int or float, b:Complex): Complex = a.complex - b 
proc `-`*(a:Complex, b:int or float): Complex = a - b.complex 
proc `*`*(a,b:Complex): Complex = [a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0]]
proc `*`*(a:Complex, b:int or float): Complex = [a[0]*b.float,a[1]*b.float]
proc `*`*(a:int or float, b:Complex): Complex = [a.float*b[0],a.float*b[1]]
proc `*`*(a,b:int or float): float = a.float*b.float 
proc `/`*(a:Complex,b:int or float): Complex = [a[0]/b.float,a[1]/b.float]
proc `/`*(a:int or float,b:Complex): Complex = a.float * conj(b)/abs(b)^2
proc `/`*(a,b:Complex): Complex = a*conj(b)/abs(b)^2
proc `^`*(a,b:int or float):float = exp(b.float*ln(a.float))
proc `^`*(a:Complex,b:int or float): Complex = 
    let 
        theta = arctan(a[1]/a[0])
        r = abs(a)
    result = r^b.float*exp(I*theta*b.float)

# --------------------exponential function--------------------------
proc exp*(a:int): float = exp(a.float)
proc exp*(a:Complex): Complex = exp(a[0])*(cos(a[1])+I*sin(a[1]))
proc sqrt*(a:Complex): Complex =
    let 
        theta = arctan(a[1]/a[0])
        r = abs(a)
    result = sqrt(r)*exp(I*theta/2)

proc ln*(a:Complex): Complex = 
    var 
        theta = arctan(a[1]/a[0])
        r = abs(a)
    if a[0]<0 and a[1]>=0:
        theta = theta + Pi
    elif a[0]<0 and a[1]<0:
        theta = theta - Pi
    result = ln(r) + I*theta

# --------------------tigonometric functions--------------------------
proc sin*(a:Complex): Complex = sin(a[0])*cosh(a[1])+I*cos(a[0])*sinh(a[1])
proc cos*(a:Complex): Complex = cos(a[0])*cosh(a[1])-I*sin(a[0])*sinh(a[1])
proc tan*(a:Complex): Complex = sin(a)/cos(a)
proc cot*(a:Complex): Complex = cos(a)/sin(a)
proc sec*(a:Complex): Complex = 1/cos(a)
proc csc*(a:Complex): Complex = 1/sin(a)

proc sinh*(a:Complex): Complex = -I*sin(I*a)
proc cosh*(a:Complex): Complex = cos(I*a)
proc tanh*(a:Complex): Complex = -I*tan(I*a)
proc coth*(a:Complex): Complex = I*cot(I*a)

proc atan*(a:Complex): Complex = ln((I-a)/(I+a))/(2*I)
proc acot*(a:Complex): Complex = ln((a+I)/(a-I))/(2*I)
proc asin*(a:Complex): Complex = 
    let 
        z = 1 - a^2
        theta = arctan(z[1]/z[0])
    result = ( ln( I*a + sqrt(abs(z)) * exp(I*theta/2) ) )/I

proc acos*(a:Complex): Complex =
    let 
        z = 1 - a^2
        theta = arctan(z[1]/z[0])
    result = ( ln( a + I*sqrt(abs(z)) * exp(I*theta/2) ) )/I

#------------------------matrix manipulation-------------------
proc zeros*(n:int): RowVec =
    for i in 0..n-1:
        result.add(0.complex)

proc zeros*(n:int,m:int): Matrix =
    for i in 0..n-1:
        result.add(@[])
        for j in 0..m-1:
            result[i].add(0.complex)

proc transpose*(vec: RowVec): ColVec =
    for i in 0..len(vec)-1:
        result.add([vec[i]])

proc transpose*(vec: ColVec): RowVec =
    for i in 0..len(vec)-1:
        result.add(vec[i][0])

proc conj_transpose*(vec: RowVec): ColVec =
    for i in 0..len(vec)-1:
        result.add([conj(vec[i])])

proc conj_transpose*(vec: ColVec): RowVec =
    for i in 0..len(vec)-1:
        result.add(conj(vec[i][0]))


proc transpose*(matrix: Matrix): Matrix =
    for i in 0..len(matrix)-1:
        result.add(@[])
        for j in 0..len(matrix[i])-1:
            result[i].add(matrix[j][i])

proc conj_transpose*(matrix: Matrix): Matrix =
    for i in 0..len(matrix)-1:
        result.add(@[])
        for j in 0..len(matrix[i])-1:
            result[i].add(conj(matrix[j][i]))


proc `+`*(a,b: RowVec): RowVec =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add(a[i] + b[i])

proc `+`*(a,b: ColVec): ColVec =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add([a[i][0] + b[i][0]])

proc `-`*(a,b: ColVec): ColVec =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add([a[i][0] - b[i][0]])

proc `-`*(a,b: RowVec): RowVec =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add(a[i] - b[i])

proc `*`*(a:RowVec, b:int or float or Complex):RowVec =
    for i in 0..len(a)-1:
        result.add(b.complex*a[i])

proc `*`*(a: int or float or Complex, b: RowVec): RowVec =
    for i in 0..len(b)-1:
        result.add(a.complex*b[i])

proc `*`*(a:ColVec, b:int or float or Complex): ColVec =
    for i in 0..len(a)-1:
        result.add([b.complex*a[i][0]])

proc `*`*(a: int or float or Complex, b: ColVec): ColVec =
    for i in 0..len(b)-1:
        result.add([a.complex*b[i][0]])

proc `/`*(a: RowVec,b: int or float or Complex): RowVec =
    if b.complex == complex(0,0):
        echo "divide by zero"
    else:
        for i in 0..len(a)-1:
            result.add(a[i]/b.complex)

proc `/`*(a: ColVec,b: int or float or Complex): ColVec =
    if b.complex == complex(0,0):
        echo "divide by zero"
    else:
        for i in 0..len(a)-1:
            result.add([a[i][0]/b.complex])

proc `/`*(matrix:Matrix, a: int or float or Complex): Matrix =
    result = matrix
    for i in 0..len(matrix)-1:
        for j in 0..len(matrix)-1:
            result[i][j] = matrix[i][j]/a.complex

proc `/`*(a: int or float or Complex, matrix: Matrix): Matrix =
    result = matrix
    for i in 0..len(matrix)-1:
        for j in 0..len(matrix)-1:
            result[i][j] = matrix[i][j]/a.complex

proc `-`*(m1, m2: Matrix): Matrix =
    result = m1
    for i in 0..len(m1)-1:
        for j in 0..len(m1)-1:
            result[i][j] = m1[i][j]-m2[i][j]

proc sign*(a:int or float): int =
    if a<0:
        result = -1
    else:
        result = 1

proc ones*(n:int): Matrix =
    result = zeros(n,n)
    for i in 0..n-1:
        result[i][i] = 1.complex


proc remove_row*(matrix: Matrix,row_num:int): Matrix =
    result = matrix
    result.delete(row_num)

proc remove_column*(matrix: Matrix,column_num:int): Matrix =
    result = matrix
    var i = 0
    while i<len(result):
        result[i].delete(column_num)
        i += 1

proc det*(matrix: Matrix): Complex =
    if len(matrix) == 1:
        result = matrix[0][0]
    else:
        let first_row = matrix[0]
        var 
            tmp = matrix.remove_row(0)
            sum = complex(0,0)
        for i in 0..len(first_row)-1:
            let reduced_matrix = tmp.remove_column(i)
            sum = sum + (-1)^i*first_row[i]*det(reduced_matrix)
        result = sum

proc abs*(vec: RowVec): float =
    result = 0
    for i in 0..len(vec)-1:
        result = result + abs(vec[i])^2
    return sqrt(result)

proc abs*(vec: ColVec): float =
    result = 0
    for i in 0..len(vec)-1:
        result = result + abs(vec[i][0])^2
    return sqrt(result)

proc `*`*(vec1: RowVec, vec2: ColVec): Complex =
    if len(vec1) != len(vec2):
        echo "wrong dimensions!"
    else:
        for i in 0..len(vec1)-1:
            result = result + vec1[i]*vec2[i][0]

proc `*`*(vec1: ColVec, vec2: RowVec): Matrix =
    if len(vec1) != len(vec2):
        echo "wrong dimensions!"
    else:
        for i in 0..len(vec1)-1:
            result.add(@[])
            for j in 0..len(vec2)-1:
                result[i].add(vec1[i][0]*vec2[j]) 

proc `*`*(m1,m2: Matrix): Matrix =
    result = zeros(len(m1),len(m2))
    let m2t = m2.transpose
    for i in 0..len(m1)-1:
        for j in 0..len(m2)-1:
            result[i][j] = result[i][j] + m1[i]*(m2t[j]).transpose

proc insert_row*(matrix: Matrix,row: RowVec, index: int): Matrix =
    for i in 0..len(matrix)-1:
        if i == index:
            result.add(row)
        result.add(matrix[i])

proc QR_Decom*(A: Matrix): array[2,Matrix] =
    var
        e1,u,v: ColVec
        norm: float
        newA,a,H,tmp: Matrix
        Hk_s: seq[Matrix]

    for i in 0..len(A)-2:
        if i == 0:
            newA = A
        else:
            newA = newA.remove_row(0).remove_column(0)
        a = transpose(newA)
        e1 = zeros(len(a)).transpose
        norm = abs(a[0])
        e1[0][0] = 1.complex
        u = a[0].transpose + (sign((a[0][0])[0])*norm)*e1
        v = u/u[0][0]
        H = ones(len(a)) - 2*v*v.transpose/(v.transpose*v)
        newA = H*newA
        if i == 0:
            Hk_s.add(H)
        else:
            tmp = ones(len(A))
            for j in i..len(A)-1:
                for k in i..len(A)-1:
                    tmp[j][k] = H[j-i][k-i]
            Hk_s.add(tmp)

    var R,Q = ones(len(A))
    for i in 0..len(Hk_s)-1:
        R = Hk_s[i]*R
        Q = Q*Hk_s[i]
    R = R*A
    # note that A = Q*R
    result[0] = Q
    result[1] = R

proc round*(matrix: Matrix, decimal: int = 0): Matrix =
    result = matrix
    for i in 0..len(matrix)-1:
        for j in 0..len(matrix)-1:
            result[i][j] = round(result[i][j],decimal)

proc upper_triangular_Q*(matrix:Matrix): bool =
    result = true
    for i in 0..len(matrix)-2:
        for j in i+1..len(matrix)-1:
            if matrix[j][i] != 0.complex:
                result = false
                break

proc diagonal_Q*(matrix:Matrix): bool =
    result = true
    for i in 0..len(matrix)-1:
        for j in 0..len(matrix)-1:
            if i != j and matrix[i][j] != 0.complex:
                result = false
                break

proc Evals*(A:Matrix, decimal_tolerance=8, steps:int=100): RowVec =
    var 
        a,previous,Q,R: Matrix
        qr: array[2,Matrix]
    a = A
    var count = 0
    while true:
        previous = a
        qr = QR_Decom(a)
        Q = qr[0]
        R = qr[1]
        a = R*Q
        if upper_triangular_Q(round(a-previous,decimal_tolerance)) or count == steps:
            if count < steps:
                echo "the algorithm converged after ", count, " steps"
            else:
                echo "the algorithm did not converge after ", count, " steps"
            break
        count += 1
    for i in 0..len(A)-1:
        result.add(a[i][i])

proc Evecs*(A:Matrix,decimal_tolerance=8, steps=100): Matrix =
    var 
        a,previous,Q,R: Matrix
        qr: array[2,Matrix]
        all_Q_s: seq[Matrix]
    a = A
    var count = 0
    while true:
        previous = a
        qr = QR_Decom(a)
        Q = qr[0]
        all_Q_s.add(Q)
        R = qr[1]
        a = R*Q
        if diagonal_Q(round(a-previous,decimal_tolerance)) or count == steps:
            if count < steps:
                echo "the algorithm converged after ", count, " steps"
            else:
                echo "the algorithm did not converge after ", count, " steps"
            break
        count += 1
    result = ones(len(A))
    for i in 0..len(all_Q_s)-1:
        result = result*all_Q_s[i]

proc `*`*(vec:RowVec, matrix:Matrix): RowVec =
    for i in 0..len(matrix)-1:
        var tmp = 0.complex
        for j in 0..len(matrix)-1:
            tmp = tmp + vec[j]*matrix[j][i]
        result.add( tmp )

proc `*`*(matrix:Matrix,vec:ColVec): ColVec =
    for i in 0..len(matrix)-1:
        var tmp = 0.complex
        for j in 0..len(matrix)-1:
            tmp = tmp + matrix[i][j]*vec[j][0]
        result.add( [tmp] )

proc diagonalize*(matrix: Matrix): Matrix =
    let q = Evecs(matrix)
    result = q.transpose*matrix*q