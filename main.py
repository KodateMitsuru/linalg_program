from sympy import Rational, Matrix, eye, zeros, sqrt

output = []

IDNum = input("请输入你的学号: ")
while len(IDNum) != 12:
    if len(IDNum) != 12:
        print("学号长度不对，请重新输入。")
    IDNum = input("请输入你的学号: ")

output.append(f"学号是: {IDNum}\n")
output.append("\n")

output.append("计算结果如下:\n")
output.append("\n")

IDNum = list(map(int, list(IDNum)))
A = zeros(13, 13)
power = 1
it = 0
length = len(IDNum)
for i in range(13):
    for j in range(13):
        A[i, j] = Rational(IDNum[it])
        it += 1
        if it == length:
            it = 0
            power += 1
            IDNum = list(map(int, list(''.join(map(str, map(pow, IDNum, [power] * length))))))
            length = len(IDNum)

# 计算行列式
det = A.det()
output.append(f"矩阵的行列式是: {det}\n")
output.append("\n")

# 将矩阵转换为规范的阶梯形矩阵

H, pivot_columns= A.rref()
rank = len(pivot_columns)
output.append("矩阵的规范阶梯形是:\n")
output.append(f"H =\n")
for row in H.tolist():
    formatted_row = [str(frac) for frac in row]
    output.append(f"{formatted_row}\n")
output.append("\n")

# 求方程AX = 0的通解
rows, cols = A.shape
if rank == cols:
    output.append("无通解。\n")
else:
    mat_free_variables = cols - rank
    E = eye(mat_free_variables)
    solve_space = (-H[:rank, rank:]).col_join(E)
    output.append("方程AX = 0的通解是:\n")
    for i in range(mat_free_variables):
        solution = solve_space[:, i]
        formatted_solution = [str(frac) for frac in solution]
        output.append(f"n{i+1} = {formatted_solution}^T\n")
    output.append(f"∑kini，i∈[1, {mat_free_variables}]，ki∈R\n")
output.append("\n")

# 找到极大无关组
independent_matrix = A[:, pivot_columns]
output.append("极大无关组是:\n")
for i in range(len(pivot_columns)):
    column = independent_matrix[:, i]
    formatted_column = [str(frac) for frac in column]
    output.append(f"n{i+1} = {formatted_column}^T\n")
output.append("\n")

# 施行 Schmidt 正交化
def schmidt_orthogonalization(matrix):
    num_columns = matrix.shape[1]
    orthogonal_matrix = zeros(matrix.shape[0], num_columns)
    
    for i in range(num_columns):
        v = matrix[:, i]
        for j in range(i):
            proj = (v.dot(orthogonal_matrix[:, j]) / orthogonal_matrix[:, j].dot(orthogonal_matrix[:, j])) * orthogonal_matrix[:, j]
            v -= proj
        orthogonal_matrix[:, i] = v / v.norm()
    
    return orthogonal_matrix

orthogonal_matrix = schmidt_orthogonalization(independent_matrix)
output.append("标准正交组是:\n")
i = 0
for row in orthogonal_matrix.tolist():
    i += 1
    if i > len(pivot_columns):
        break
    formatted_row = [str(frac).replace('sqrt', '√').replace('*', '') for frac in row]
    output.append(f"n{i} = {formatted_row}^T\n")
output.append("\n")

# 计算A + AT的正负惯性指数
A_AT = A + A.T
eigenvalues = A_AT.eigenvals()
pos_inertia = sum(1 for val in eigenvalues if val > 0)
neg_inertia = sum(1 for val in eigenvalues if val < 0)

output.append(f"A + AT的正惯性指数是: {pos_inertia}\n")
output.append("\n")
output.append(f"A + AT的负惯性指数是: {neg_inertia}\n")

# 将输出写入文件
with open("output.txt", "w", encoding="utf-8") as f:
    f.writelines(output)