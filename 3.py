from math import sin, pi
E = 0.001
N = 1000
TRAPEZOIDAL_OMEGA = 1/3
x_g = [-0.5773503,
0.5773503]
c_g = [1, 1]
def F(x, y):
    return 2*x*sin(x*y)

def g_gauss(x):
    return calculate_gauss(lambda y: F(x,y), 0, pi/2, N)

def g_trap(x):
    return calculate_trapezoidal(lambda y: F(x,y), 0, pi/2, E, TRAPEZOIDAL_OMEGA)

def runge_check(I2n, In, omega, E):
    return omega*abs(I2n - In) < E

def calculate_cells(f, a, b, c, d, n, m):
    step1 = (b-a)/n
    step2 = (d-c)/m
    x_i1 = a
    x_i2 = x_i1+step1
    sum_I = 0
    while x_i1 < b:
        y_i1 = c
        y_i2 = y_i1+step2
        while y_i1 < d:
            sum_I += f((x_i2+x_i1)/2, (y_i2+y_i1)/2)
            y_i1 = y_i2
            y_i2 += step2
        x_i1 = x_i2
        x_i2 += step1
        
    return sum_I*step1*step2

def sum_integral(f, a, b, n):
    step = (b-a)/n
    x_i1 = a
    x_i2 = x_i1+step
    sum_I = 0
    while x_i1 < b:
        sum_I += (f(x_i1)+f(x_i2) )/2 * step
        x_i1 = x_i2
        x_i2 += step
    return sum_I

def calculate_trapezoidal(f, a, b, E, omega):
    n = 2
    In = sum_integral(f, a, b, n)
    n*=2
    I2n = sum_integral(f, a, b, n)
    while not runge_check(I2n, In, omega, E):
        In = I2n
        n *= 2
        I2n = sum_integral(f, a, b, n)  
    return I2n

def calculate_gauss(f, a, b, n):
    step = (b-a)/n
    x = a+step
    sum_I = 0   
    while x < b:
        sum_I += gauss_method(f, x, step)
        x += step*2
    return step*sum_I


def sum_integralG(f, a, b, n):
    step = (b-a)/n
    x = a+step
    sum_I = 0   
    while x < b:
        sum_I += gauss_method(f, x, step)
        x += step
    return step*sum_I/2

def gauss_method(f, xj_1, step):
    sum = 0
    xj_0 = xj_1-step
    for i in range(len(c_g)):
        xi = xj_0 + (x_g[i]+1)*(xj_0-xj_1)/2
        sum += c_g[i]*f(xi)
    return sum
        
result_gauss_y = []

result = calculate_cells(F, 0, 1, 0, pi/2, 20, 20)
print("calculate_cells: ", result)
result = calculate_trapezoidal(g_trap, 0, 1, E, TRAPEZOIDAL_OMEGA)
print("calculate_trapezoidal: ", result)
result = calculate_gauss(g_gauss, 0, 1, N)
print("calculate_gauss: ", result)