import strutils, mathexpr, math

type
    Complex* = object
        re*, im*: float
    RVecC* = seq[Complex]
    CVecC* = seq[array[1,Complex]]
    MatrixC* = seq[RVecC]
    
    RVec* = seq[float]
    CVec* = seq[array[1,float]]
    Matrix* = seq[RVec]

const
    I* = Complex(im:1.0)

proc abs*(a: int or float): float = 
    if a < 0:
        result = -a
    else:
        result = a

proc sign*(a:int or float): int =
    if a<0: result = -1
    else: result = 1

#----------------------------------------PreDeclaration of some functions----------------------------------
proc ln*(a:Complex): Complex
proc exp*(a:Complex): Complex
proc evalc*(input: string): Complex
proc remove_row*(mat: Matrix,row_num:int): Matrix
proc remove_row*(mat: MatrixC,row_num:int): MatrixC
proc remove_column*(mat: Matrix,column_num:int): Matrix
proc remove_column*(mat: MatrixC,column_num:int): MatrixC


# ----------------------inside complex operation----------------------
proc complex*(a:int or float): Complex = 
    result.re = a.float
proc complex*(a:Complex): Complex = a
proc complex*(a,b:int or float): Complex = 
    result.re = a.float
    result.im = b.float
proc complex*(a:array[2,int or float]): Complex = 
    result.re = a.re.float
    result.im = a.im.float

#-----------------normal arithmetic operations-----------------
proc abs*(a:Complex): float = sqrt(a.re^2+a.im^2)
proc norm*(z: Complex): float = abs(z)^2 
proc conj*(a:Complex): Complex = 
    result.re = a.re
    result.im = -a.im

proc round*(a:Complex,b:int=0):Complex = 
    result.re = round(a.re,b)
    result.im = round(a.im,b)

proc `+`*(a:Complex): Complex = a
proc `+`*(a,b:Complex): Complex = 
    result.re = a.re + b.re
    result.im = a.im + b.im
proc `+`*(a,b:int or float): float = a.float + b.float 
proc `+`*(a:int or float, b:Complex): Complex = a.complex + b
proc `+`*(a:Complex, b:int or float): Complex = a + b.complex
proc `-`*(a:Complex): Complex = 
    result.re = -a.re
    result.im = -a.im
proc `-`*(a,b:Complex): Complex = 
    result.re = a.re - b.re
    result.im = a.im - b.im
proc `-`*(a:int or float, b:Complex): Complex = a.complex - b 
proc `-`*(a:Complex, b:int or float): Complex = a - b.complex 
proc `*`*(a,b:Complex): Complex = 
    result.re = a.re * b.re - a.im * b.im
    result.im = a.re * b.im + a.im * b.re
proc `*`*(a:Complex, b:int or float): Complex = 
    result.re = a.re * b.float
    result.im = a.im * b.float
proc `*`*(a:int or float, b:Complex): Complex = 
    result.re = a.float * b.re
    result.im = a.float * b.im
proc `*`*(a,b:int or float): float = a.float * b.float 
proc `/`*(a:Complex,b:int or float): Complex = 
    result.re = a.re / b.float 
    result.im = a.im / b.float
proc `/`*(a:int or float,b:Complex): Complex = a.float * conj(b)/abs(b)^2
proc `/`*(a,b:Complex): Complex = a*conj(b)/abs(b)^2

proc `^`*(a:int or float, b: int): float = 
    var s = 1
    if a < 0 and b mod 2 != 0: s = -1
    if a != 0:
        result = s * exp(b.float*ln(abs(a.float)))
    else:
        if b <= 0: result = NaN

proc `^`*(a:int or float, b:float): Complex = 
    var s = 0
    if a < 0: s = 1
    if a != 0: 
        result = exp(b.float*(2*s*ln(I) + ln(abs(a.float))))
    else:
        if b <= 0: result = NaN.complex

proc `^`*(a,b:Complex): Complex = exp(b*ln(a))
proc `^`*(a: int or float, b:Complex): Complex = exp(b*ln(a.float))
proc `^`*(a:Complex,b:int or float): Complex = exp(b.float*ln(a))
    #let 
    #    theta = arctan(a.im/a.re)
    #    r = abs(a)
    #result = r^b.float*exp(I*theta*b.float)

proc `+=`*(a: var Complex, b: int or float or Complex) = a = a + b

# --------------------exponential function--------------------------
proc exp*(a:int): float = exp(a.float)
proc exp*(a:Complex): Complex = exp(a.re)*(cos(a.im)+I*sin(a.im))

proc sqrt*(a:int): float = sqrt(a.float)

proc sqrt*(a:Complex): Complex =
    let 
        theta = arctan(a.im/a.re)
        r = abs(a)
    result = sqrt(r)*exp(I*theta/2)

proc ln*(a:Complex): Complex = 
    var 
        theta = arctan(a.im/a.re)
        r = abs(a)
    if a.re<0 and a.im>=0:
        theta = theta + Pi
    elif a.re<0 and a.im<0:
        theta = theta - Pi
    result = ln(r) + I*theta

# --------------------tigonometric functions--------------------------
proc sin*(a:Complex): Complex = sin(a.re)*cosh(a.im)+I*cos(a.re)*sinh(a.im)
proc cos*(a:Complex): Complex = cos(a.re)*cosh(a.im)-I*sin(a.re)*sinh(a.im)
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
        theta = arctan(z.im/z.re)
    result = ( ln( I*a + sqrt(abs(z)) * exp(I*theta/2) ) )/I

proc acos*(a:Complex): Complex =
    let 
        z = 1 - a^2
        theta = arctan(z.im/z.re)
    result = ( ln( a + I*sqrt(abs(z)) * exp(I*theta/2) ) )/I

#------------------------Real matrix manipulation-------------------
proc rvec*(a: seq[int or float]): RVec = 
    for i in 0..len(a)-1:
        result.add(a[i].float)

proc cvec*(a: seq[int or float]): CVec = 
    for i in 0..len(a)-1:
        result.add([a[i].float])

proc matrix*(a: seq[seq[int or float]]): Matrix = 
    for i in 0..len(a)-1:
        result.add(@[])
        for j in 0..len(a[i])-1:
            result[i].add(a[i][j].float)

proc shape*(vec: RVec or CVec or RVecC or CVecC): int = len(vec)

proc shape*(mat: Matrix or MatrixC): tuple = (len(mat), len(mat[0]))

proc abs*(vec: RVec): float =
    result = 0
    for i in 0..len(vec)-1:
        result += abs(vec[i])^2
    result = sqrt(result)

proc abs*(vec: CVec): float =
    result = 0
    for i in 0..len(vec)-1:
        result += abs(vec[i][0])^2
    result = sqrt(result)

proc zeros*(n:int): RVec =
    for i in 0..n-1:
        result.add(0.0)

proc zeros*(n:int,m:int): Matrix =
    for i in 0..n-1:
        result.add(@[])
        for j in 0..m-1:
            result[i].add(0.0)

proc ones*(n:int): RVec =
    result = zeros(n)
    for i in 0..n-1:
        result[i] = 1.0

proc ones*(n,m:int): Matrix =
    result = zeros(n,m)
    for i in 0..n-1:
        for j in 0..m-1:
            result[i][j] = 1.0

proc eye*(n:int): Matrix =
    result = zeros(n,n)
    for i in 0..n-1:
        result[i][i] = 1.0

proc norm*(vec: RVec or CVec): float = 
    for i in 0..vec.len()-1:
        result += vec[i]^2

proc transpose*(vec: RVec): CVec =
    for i in 0..len(vec)-1:
        result.add([vec[i]])

proc transpose*(vec: CVec): RVec =
    for i in 0..len(vec)-1:
        result.add(vec[i][0])

proc transpose*(mat: Matrix): Matrix =
    for i in 0..len(mat)-1:
        result.add(@[])
        for j in 0..len(mat[i])-1:
            result[i].add(mat[j][i])

proc round*(vec: RVec, decimal: int = 0): RVec =
    result = vec
    for i in 0..len(vec)-1:
        result[i] = round(vec[i], decimal)

proc round*(vec: CVec, decimal: int = 0): CVec =
    result = vec
    for i in 0..len(vec)-1:
        result[i][0] = round(vec[i][0], decimal)

proc round*(mat: Matrix, decimal: int = 0): Matrix =
    result = mat
    for i in 0..len(mat)-1:
        for j in 0..len(mat[i])-1:
            result[i][j] = round(mat[i][j], decimal)

proc mat2rvec*(mat: Matrix): RVec =
    for i in 0..mat.len() - 1:
        for j in 0..mat[i].len() - 1:
            result.add(mat[i][j])

proc reshape*(vec: RVec, rows, cols: int): Matrix =
    result = zeros(rows,cols)
    var tmp = vec
    if rows*cols != tmp.len():
        echo "Incompattible shape"
    else:
        var k = 0
        for i in 0..rows-1:
            for j in 0..cols-1:
                result[i][j] = tmp[k]
                inc k

proc reshape*(mat: Matrix, rows, cols: int): Matrix =
    result = zeros(rows,cols)
    var tmp = mat2rvec(mat)
    if rows*cols != tmp.len():
        echo "Incompattible shape"
    else:
        var k = 0
        for i in 0..rows-1:
            for j in 0..cols-1:
                result[i][j] = tmp[k]
                inc k

proc det*(mat: Matrix): float =
    if mat.len == 1:
        result = mat[0][0]
    else:
        let first_row = mat[0]
        var 
            tmp = mat.remove_row(0)
            sum = 0.0
        for i in 0..len(first_row)-1:
            let reduced_MatrixC = tmp.remove_column(i)
            sum = sum + (-1)^i*first_row[i]*det(reduced_MatrixC)
        result = sum

proc remove_row*(mat: Matrix,row_num:int): Matrix =
    result = mat
    result.delete(row_num)

proc remove_column*(mat: Matrix,column_num:int): Matrix =
    result = mat
    var i = 0
    while i<len(result):
        result[i].delete(column_num)
        inc i

proc insert_row*(mat: Matrix,row: RVec, index: int): Matrix =
    for i in 0..mat.len-1:
        if i == index:
            result.add(row)
        result.add(mat[i])

proc str2rvecf*(input: string): RVec = 
    var tmp = input[1..^2].split(" ")
    for i in 0..tmp.len()-1:
        result.add( tmp[i].parseFloat )

proc evalf*(input: string): float =
    let e = newEvaluator()
    result = e.eval(input)

#---------------------------------- Real matrix operations-----------------------------------
proc `+`*(a,b: RVec): RVec =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add(a[i] + b[i])

proc `+`*(a,b: CVec): CVec =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add([a[i][0] + b[i][0]])

proc `+`*(m1, m2: Matrix): Matrix =
    result = m1
    for i in 0..len(m1)-1:
        for j in 0..len(m1)-1:
            result[i][j] = m1[i][j] + m2[i][j]

proc `-`*(a,b: CVec): CVec =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add([a[i][0] - b[i][0]])

proc `-`*(a,b: RVec): RVec =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add(a[i] - b[i])

proc `-`*(m1, m2: Matrix): Matrix =
    result = m1
    for i in 0..len(m1)-1:
        for j in 0..len(m1)-1:
            result[i][j] = m1[i][j]-m2[i][j]

proc `*`*(a:RVec, b:int or float):RVec =
    for i in 0..len(a)-1:
        result.add(b.float * a[i])

proc `*`*(a: int or float, b: RVec): RVec =
    for i in 0..len(b)-1:
        result.add(a.float * b[i])

proc `*`*(a:CVec, b:int or float): CVec =
    for i in 0..len(a)-1:
        result.add([b.float * a[i][0]])

proc `*`*(a: int or float, b: CVec): CVec =
    for i in 0..len(b)-1:
        result.add([a.float * b[i][0]])

proc `*`*(vec1: RVec, vec2: CVec): float =
    if len(vec1) != len(vec2):
        echo "Error!!, Wrong Dimensions!"
    else:
        for i in 0..len(vec1)-1:
            result += vec1[i]*vec2[i][0]

proc `*`*(vec1: CVec, vec2: RVec): Matrix =
    result = zeros(len(vec1),len(vec2))
    for i in 0..len(vec1)-1:
        for j in 0..len(vec2)-1:
            result[i][j] = vec1[i][0] * vec2[j]

proc `*`*(m1,m2: Matrix): Matrix =
    result = zeros(len(m1),len(m2))
    let m2t = m2.transpose
    for i in 0..len(m1)-1:
        for j in 0..len(m2)-1:
            result[i][j] += m1[i]*(m2t[j]).transpose

proc `*`*(mat:Matrix, vec:CVec): CVec =
    if mat[0].len != vec.len:
        echo "Incompattible shapes"
    else:
        for i in 0..len(mat)-1:
            var tmp = 0.0
            for j in 0..len(mat)-1:
                tmp += mat[i][j]*vec[j][0]
            result.add( [tmp] )

proc `*`*(vec:RVec, mat:Matrix): RVec =
    for i in 0..len(mat)-1:
        var tmp = 0.0
        for j in 0..len(mat)-1:
            tmp += vec[j]*mat[j][i]
        result.add( tmp )

proc `/`*(a: RVec,b: int or float): RVec =
    if b.float == 0.0:
        echo "Error!!, divide by zero"
    else:
        for i in 0..len(a)-1:
            result.add(a[i]/b.float)

proc `/`*(a: CVec,b: int or float): CVec =
    if b.float == 0.0:
        echo "Error!!, divide by zero"
    else:
        for i in 0..len(a)-1:
            result.add([a[i][0]/b.float])

proc `/`*(mat:Matrix, a: int or float): Matrix =
    result = mat
    for i in 0..len(mat)-1:
        for j in 0..len(mat)-1:
            result[i][j] = mat[i][j]/a.float

proc `/`*(a: int or float, mat: Matrix): Matrix =
    result = mat
    for i in 0..len(mat)-1:
        for j in 0..len(mat)-1:
            result[i][j] = mat[i][j]/a.float

#------------------------Complex matrix manipulation-------------------
proc rvecc*(a: seq[int or float or Complex]): RVecC = 
    for i in 0..len(a)-1:
        result.add(a[i].complex)

proc cvecc*(a: seq[int or float or Complex]): CVecC = 
    for i in 0..len(a)-1:
        result.add([a[i].complex])

proc matrixc*(a: seq[seq[int or float or Complex]]): MatrixC = 
    for i in 0..len(a)-1:
        result.add(@[])
        for j in 0..len(a[i])-1:
            result[i].add(a[i][j].complex)

proc abs*(vec: RVecC): float =
    result = 0
    for i in 0..len(vec)-1:
        result += abs(vec[i])^2
    result = sqrt(result)

proc abs*(vec: CVecC): float =
    result = 0
    for i in 0..len(vec)-1:
        result += abs(vec[i][0])^2
    result = sqrt(result)

proc zerosc*(n:int): RVecC =
    for i in 0..n-1:
        result.add(0.complex)

proc zerosc*(n:int,m:int): MatrixC =
    for i in 0..n-1:
        result.add(@[])
        for j in 0..m-1:
            result[i].add(0.complex)

proc onesc*(n:int): MatrixC =
    result = zerosc(n,n)
    for i in 0..n-1:
        result[i][i] = 1.complex

proc onesc*(n,m:int): MatrixC =
    result = zerosc(n,m)
    for i in 0..n-1:
        for j in 0..m-1:
            result[i][j] = 1.complex

proc transpose*(vec: RVecC): CVecC =
    for i in 0..len(vec)-1:
        result.add([vec[i]])

proc transpose*(vec: CVecC): RVecC =
    for i in 0..len(vec)-1:
        result.add(vec[i][0])

proc conj_transpose*(vec: RVecC): CVecC =
    for i in 0..len(vec)-1:
        result.add([conj(vec[i])])

proc conj_transpose*(vec: CVecC): RVecC =
    for i in 0..len(vec)-1:
        result.add(conj(vec[i][0]))

proc transpose*(MatrixC: MatrixC): MatrixC =
    for i in 0..len(MatrixC)-1:
        result.add(@[])
        for j in 0..len(MatrixC[i])-1:
            result[i].add(MatrixC[j][i])

proc conj_transpose*(MatrixC: MatrixC): MatrixC =
    for i in 0..len(MatrixC)-1:
        result.add(@[])
        for j in 0..len(MatrixC[i])-1:
            result[i].add(conj(MatrixC[j][i]))

proc round*(mat: MatrixC, decimal: int = 0): MatrixC =
    result = mat
    for i in 0..len(mat)-1:
        for j in 0..len(mat[i])-1:
            result[i][j] = round(result[i][j],decimal)

proc round*(vec: CVecC, decimal: int = 0): CVecC =
    result = vec
    for i in 0..len(vec)-1:
        result[i][0] = round(result[i][0],decimal)

proc round*(vec: RVecC, decimal: int = 0): RVecC =
    result = vec
    for i in 0..len(vec)-1:
        result[i] = round(result[i],decimal)

proc reshape*(vec: RVecC, rows, cols: int): MatrixC =
    result = zerosc(rows,cols)
    var tmp = vec
    if rows*cols != tmp.len():
        echo "Incompattible shape"
    else:
        var k = 0
        for i in 0..rows-1:
            for j in 0..cols-1:
                result[i][j] = tmp[k]
                inc k

proc det*(mat: MatrixC): Complex =
    if mat.len == 1:
        result = mat[0][0]
    else:
        let first_row = mat[0]
        var 
            tmp = mat.remove_row(0)
            sum = complex(0,0)
        for i in 0..len(first_row)-1:
            let reduced_MatrixC = tmp.remove_column(i)
            sum = sum + (-1)^i*first_row[i]*det(reduced_MatrixC)
        result = sum

proc mat2rvecc*(mat: MatrixC): RVecC =
    for i in 0..mat.len() - 1:
        for j in 0..mat[i].len() - 1:
            result.add(mat[i][j])

proc reshape*(mat: MatrixC, rows, cols: int): MatrixC =
    result = zerosc(rows,cols)
    var tmp = mat2rvecc(mat)
    if rows*cols != tmp.len():
        echo "Incompattible shape"
    else:
        var k = 0
        for i in 0..rows-1:
            for j in 0..cols-1:
                result[i][j] = tmp[k]
                inc k

proc remove_row*(mat: MatrixC,row_num:int): MatrixC =
    result = mat
    result.delete(row_num)

proc remove_column*(mat: MatrixC,column_num:int): MatrixC =
    result = mat
    var i = 0
    while i<len(result):
        result[i].delete(column_num)
        inc i

proc insert_row*(mat: MatrixC,row: RVecC, index: int): MatrixC =
    for i in 0..mat.len-1:
        if i == index:
            result.add(row)
        result.add(mat[i])

proc str2rvecc*(input: string): RVecC = 
    var tmp = input[1..^2].split(" ")
    for i in 0..tmp.len()-1:
        result.add( evalc(tmp[i]) )

proc str2cvecc*(input: string): CVecC = 
    var tmp = input[1..^2].split(" ")
    for i in 0..tmp.len()-1:
        result.add( [evalc(tmp[i])] )

proc evalc*(input: string): Complex =
    let e = newEvaluator()
    e.addVar("I", 0.0)
    e.addVar("i", 0.0)
    result.re = e.eval(input)
    e.addVar("I", 1.0)
    e.addVar("i", 1.0)
    result.im = e.eval(input) - result.re

proc evalc*(input: seq[string]): seq[Complex] =
    result = @[]
    var element: Complex
    for i in 0..input.len()-1:
        let e = newEvaluator()
        e.addVar("I", 0.0)
        e.addVar("i", 0.0)
        element.re = e.eval(input[i])
        e.addVar("I", 1.0)
        e.addVar("i", 1.0)
        element.im = e.eval(input[i]) - element.re
        result.add( element )

proc trace*(mat: Matrix): float =
    result = 0.0
    for i in 0..mat.len()-1:
        result += mat[i][i]

proc trace*(mat: MatrixC): Complex =
    result = 0.complex
    for i in 0..mat.len()-1:
        result += mat[i][i]

#-------------------------------------Complex matrix operations-----------------------
proc `+`*(a,b: RVecC): RVecC =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add(a[i] + b[i])

proc `+`*(a,b: CVecC): CVecC =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add([a[i][0] + b[i][0]])

proc `-`*(a,b: CVecC): CVecC =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add([a[i][0] - b[i][0]])

proc `-`*(a,b: RVecC): RVecC =
    if len(a) != len(b):
        echo "wrong dimensions!"
    else:
        for i in 0..len(a)-1:
            result.add(a[i] - b[i])

proc `-`*(m1, m2: MatrixC): MatrixC =
    result = m1
    for i in 0..len(m1)-1:
        for j in 0..len(m1)-1:
            result[i][j] = m1[i][j]-m2[i][j]

proc `*`*(a:RVecC, b:int or float or Complex):RVecC =
    for i in 0..len(a)-1:
        result.add(b.complex*a[i])

proc `*`*(a: int or float or Complex, b: RVecC): RVecC =
    for i in 0..len(b)-1:
        result.add(a.complex*b[i])

proc `*`*(a:CVecC, b:int or float or Complex): CVecC =
    for i in 0..len(a)-1:
        result.add([b.complex*a[i][0]])

proc `*`*(a: int or float or Complex, b: CVecC): CVecC =
    for i in 0..len(b)-1:
        result.add([a.complex*b[i][0]])

proc `*`*(mat:MatrixC, vec:CVecC): CVecC =
    if mat[0].len != vec.len:
        echo "Incompattible shapes"
    else:
        for i in 0..len(mat)-1:
            var tmp = 0.complex
            for j in 0..len(mat)-1:
                tmp += mat[i][j]*vec[j][0]
            result.add( [tmp] )

proc `*`*(vec:RVecC, mat:MatrixC): RVecC =
    if mat.len != vec.len:
        echo "Incompattible shapes"
    else:
        for i in 0..len(mat)-1:
            var tmp = 0.complex
            for j in 0..len(mat)-1:
                tmp += vec[j]*mat[j][i]
            result.add( tmp )

proc `*`*(vec1: RVecC, vec2: CVecC): Complex =
    if len(vec1) != len(vec2):
        echo "Incompattible shapes"
    else:
        for i in 0..len(vec1)-1:
            result = result + vec1[i]*vec2[i][0]

proc `*`*(vec1: CVecC, vec2: RVecC): MatrixC =
    if len(vec1) != len(vec2):
        echo "wrong dimensions!"
    else:
        for i in 0..len(vec1)-1:
            result.add(@[])
            for j in 0..len(vec2)-1:
                result[i].add(vec1[i][0]*vec2[j]) 

proc `*`*(m1,m2: MatrixC): MatrixC =
    result = zerosc(len(m1),len(m2))
    let m2t = m2.transpose
    for i in 0..len(m1)-1:
        for j in 0..len(m2)-1:
            result[i][j] = result[i][j] + m1[i]*(m2t[j]).transpose

proc `/`*(a: RVecC,b: int or float or Complex): RVecC =
    if b.complex == complex(0,0):
        echo "divide by zero"
    else:
        for i in 0..len(a)-1:
            result.add(a[i]/b.complex)

proc `/`*(a: CVecC,b: int or float or Complex): CVecC =
    if b.complex == complex(0,0):
        echo "divide by zero"
    else:
        for i in 0..len(a)-1:
            result.add([a[i][0]/b.complex])

proc `/`*(MatrixC:MatrixC, a: int or float or Complex): MatrixC =
    result = MatrixC
    for i in 0..len(MatrixC)-1:
        for j in 0..len(MatrixC)-1:
            result[i][j] = MatrixC[i][j]/a.complex

proc `/`*(a: int or float or Complex, MatrixC: MatrixC): MatrixC =
    result = MatrixC
    for i in 0..len(MatrixC)-1:
        for j in 0..len(MatrixC)-1:
            result[i][j] = MatrixC[i][j]/a.complex


proc upper_triangular_Q*(MatrixC:MatrixC): bool =
    result = true
    for i in 0..len(MatrixC)-2:
        for j in i+1..len(MatrixC)-1:
            if MatrixC[j][i] != 0.complex:
                result = false
                break

proc diagonal_Q*(MatrixC:MatrixC): bool =
    result = true
    for i in 0..len(MatrixC)-1:
        for j in 0..len(MatrixC)-1:
            if i != j and MatrixC[i][j] != 0.complex:
                result = false
                break

#[ 
proc QR_Decom*(A: MatrixC): array[2,MatrixC] =
    var
        e1,u,v: CVecC
        norm: float
        newA,a,H,tmp: MatrixC
        Hk_s: seq[MatrixC]

    for i in 0..len(A)-2:
        if i == 0:
            newA = A
        else:
            newA = newA.remove_row(0).remove_column(0)
        a = transpose(newA)
        e1 = zeros(len(a)).transpose
        norm = abs(a[0])
        e1[0][0] = 1.complex
        u = a[0].transpose + (sign((a[0][0]).re)*norm)*e1
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

 
proc Evals*(A:MatrixC, decimal_tolerance=8, steps:int=100): RVecC =
    var 
        a,previous,Q,R: MatrixC
        qr: array[2,MatrixC]
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

proc Evecs*(A:MatrixC,decimal_tolerance=8, steps=100): MatrixC =
    var 
        a,previous,Q,R: MatrixC
        qr: array[2,MatrixC]
        all_Q_s: seq[MatrixC]
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
 
proc diagonalize*(MatrixC: MatrixC): MatrixC =
    let q = Evecs(MatrixC)
    result = q.transpose*MatrixC*q
]#

proc column*(mat: Matrix, col: int): CVec =
    result = zeros(mat.len).cvec
    for i in 0..mat.len-1:
        result[i] = [mat[i][col]]

proc column*(mat: MatrixC, col: int): CVecC =
    result = zeros(mat.len).cvecc
    for i in 0..mat.len-1:
        result[i] = [mat[i][col]]

# thesis codes
proc rep*(a: Complex): string = 
    if a.im < 0:
        result = $a.re & $a.im & "i"
    else:
        result = $a.re & "+" & $a.im & "i"

proc `^`*(mat: Matrix, x: int): Matrix = 
    result = mat
    for i in 1..x-1:
        result = result * mat 

proc `^`*(mat: MatrixC, x: int): MatrixC = 
    result = mat
    for i in 1..x-1:
        result = result * mat 

#[ 
proc fisher*(eval: Matrix, evec: MatrixC, d_a_eval, d_b_eval: Matrix, d_a_evec, d_b_evec: MatrixC): float =
    var 
        decimal = 8
        fisher1 = 0.0
        fisher2 = 0.0
        fisher3 = 0.0
    for i in 0..2:
        if round(eval[i][i], decimal) != 0:
            fisher1 += d_a_eval[i][i] * d_b_eval[i][i] / round(eval[i][i], decimal)
            fisher2 += eval[i][i] * ( conj_transpose(d_a_evec.column(i)) * d_b_evec.column(i) ).re
            for j in 0..2:
                if round(eval[j][j], decimal) != 0:
                    fisher3 += eval[i][i] * eval[j][j] * ( conj_transpose(d_a_evec.column(i)) * evec.column(j) * conj_transpose(evec.column(j)) * d_b_evec.column(i) ).re / round( eval[i][i] + eval[j][j], decimal )
    result = fisher1 + 4 * fisher2 - 8 * fisher3

]#