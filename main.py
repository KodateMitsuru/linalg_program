from sympy import *
import subprocess
init_printing()
output = []

IDNum = input("请输入你的学号: ")
while len(IDNum) != 12:
    if len(IDNum) != 12:
        print("学号长度不对，请重新输入。")
    IDNum = input("请输入你的学号: ")
output.append(r'''
\documentclass[UTF8]{ctexart}
\usepackage{amsmath}
\begin{document}
\begin{sloppypar}
\title{\textbf{\huge %(school)s} \\ \textbf{\normalsize %(title)s}}
\author{name \\ 学号: %(IDNum)s}
\date{\today}
\maketitle
\pagenumbering{roman}
\tableofcontents
\newpage
\pagenumbering{arabic}
\section{计算行列式|A|}
''' % {'school': '上海交通大学', 'title': '线性代数编程作业', 'IDNum': IDNum})

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
output.append(f"$\\left | \\textbf{{A}} \\right | = {det}$\\\\\n")
output.append("\n")

# 将矩阵转换为规范的阶梯形矩阵

H, pivot_columns= A.rref()
rank = len(pivot_columns)
output.append("\\section{做行初等变换求\\textbf{A}的规范的阶梯形矩阵\\textbf{H}和方程\\textbf{AX = 0}的通解}\n")
output.append("\\subsection{阶梯形矩阵\\textbf{H}}\n")
output.append("矩阵的规范阶梯形是:\\\\\n")
output.append(f"$\\textbf{{H}} = {latex(H)}$\\\\\n")

# 求方程AX = 0的通解
output.append("\\subsection{方程\\textbf{AX = 0}的通解}\n")
rows, cols = A.shape
if rank == cols:
    output.append("无通解。\\\\\n")
else:
    mat_free_variables = cols - rank
    E = eye(mat_free_variables)
    solve_space = (-H[:rank, rank:]).col_join(E)
    output.append("方程\\textbf{AX = 0}的通解是:\\\\\n")
    for i in range(mat_free_variables):
        solution = solve_space[:, i]
        formatted_solution = latex(solution, mat_str="matrix")
        output.append(f"\\mbox{{$\\textbf{{n}}_{{{i+1}}} = {formatted_solution}$}}\n")
    output.append("\n")
    output.append(f"$\\sum\\limits_{{i=1}}^{{{mat_free_variables}}} \\left ( k_i\\textbf{{n}}_i \\right ) ,i\\in \\left [ 1,{mat_free_variables} \\right ] ,k_i\\in R$\\\\\n")
output.append("\n")

output.append("\\section{求\\textbf{A}的列向量组的一个极大无关组及Schmidt正交化求标准正交组}\n")
# 找到极大无关组
output.append("\\subsection{极大无关组}\n")
independent_matrix = A[:, pivot_columns]
output.append("极大无关组是:\\\\\n")

for i in range(len(pivot_columns)):
    column = independent_matrix[:, i]
    formatted_column = latex(column, mat_str="matrix")
    output.append(f"\\mbox{{$\\textbf{{n}}_{{{i+1}}} = {formatted_column}$}}\n")
output.append("\\newline\n")

# 施行 Schmidt 正交化
output.append("\\subsection{标准正交组}\n")
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
output.append("标准正交组是:\\\\\n")
for i in range(len(pivot_columns)):
    formatted_col = latex(orthogonal_matrix[:, i], mat_str="matrix")
    output.append(f"\\mbox{{$\\textbf{{n}}_{{{i+1}}} = {formatted_col}$}}\n")
output.append("\\\\\n")

# 计算A + AT的正负惯性指数
output.append("\\section{计算A + AT的正负惯性指数}\n")
A_AT = A + A.T
eigenvalues = A_AT.eigenvals()
pos_inertia = sum(1 for val in eigenvalues if val > 0)
neg_inertia = sum(1 for val in eigenvalues if val < 0)

output.append(f"\\textbf{{A + AT}}的正惯性指数是: ${pos_inertia}$\\\\\n")
output.append("\n")
output.append(f"\\textbf{{A + AT}}的负惯性指数是: ${neg_inertia}$\\\\\n")

output.append(r'''
\end{sloppypar}
\end{document}
''')

# 将输出写入文件
with open("output.tex", "w", encoding="utf-8") as f:
    f.writelines(output)
    
# 使用 latexindent 格式化 LaTeX 文件
subprocess.run(['latexindent', '-w', 'output.tex', '-l'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# 使用 xelatex 编译 LaTeX 文档两次
for _ in range(2):
    proc = subprocess.Popen(['xelatex', 'output.tex'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        print(stderr.decode('utf-8'))