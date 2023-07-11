from math import sin, pi
E = 0.0001
GAUSS_OMEGA = 1/E
RECTANGLES_OMEGA = 1/3 
TRAPEZOIDAL_OMEGA = 1/3
SIMPSON_OMEGA = 1/15
x_g = [-0.9324700,
    -0.6612094,
    -0.2386142,
    0.2386142,
    0.6612094,
    0.9324700]
c_g = [0.1713245,
    0.3607616,
    0.4679140,
    0.4679140,
    0.3607616,
    0.1713245]

def F(x):
    return sin(x)/(2+sin(x))

def runge_check(I2n, In, omega, E):
    return omega*abs(I2n - In) < E

def sum_integralRT(f, a, b, n, method):
    step = (b-a)/n
    x_i1 = a
    x_i2 = x_i1+step
    sum_I = 0
    while x_i1 < b:
        sum_I += method(f, x_i1, x_i2, step)
        x_i1 = x_i2
        x_i2 += step
    return sum_I

def sum_integralS(f, a, b, n, method):
    step = (b-a)/n
    x = a+step
    sum_I = 0
    while x < b:
        sum_I += method(f, x-step, x, x+step)
        x += step*2
   
    return sum_I*step/3

def sum_integralG(f, a, b, n):
    step = (b-a)/n
    x = a+step
    sum_I = 0   
    while x < b:
        sum_I += gauss_method(f, x, step)
        x += step
    return step*sum_I/2

def rectangles_method(f, x_i1, x_i2, step):
    return f( (x_i1+x_i2)/2 ) * step

def trapezoidal_method(f, x_i1, x_i2, step):
    return (f(x_i1)+f(x_i2) )/2 * step

def simpson_method(f, x_i1, x_i2, x_i3):
    return f(x_i1) + 4 * f(x_i2) + f(x_i3) 

def gauss_method(f, xj_1, step):
    sum = 0
    xj_0 = xj_1-step
    for i in range(len(c_g)):
        xi = xj_0 + (x_g[i]+1)*(xj_0-xj_1)/2
        sum += c_g[i]*f(xi)
    return sum
        
def calculate(f, a, b, E, omega, method):
    sum_integral = sum_integralRT
    if method == simpson_method:
        sum_integral = sum_integralS
    n = 2
    In = sum_integral(f, a, b, n, method)
    n*=2
    I2n = sum_integral(f, a, b, n, method)
    steps = 1
    while not runge_check(I2n, In, omega, E):
        In = I2n
        n *= 2
        I2n = sum_integral(f, a, b, n, method)
        steps += 1    
    return I2n, steps

def calculate_gauss(f, a, b, E):
    n = 2
    In = sum_integralG(f, a, b, n)
    n*=2
    I2n = sum_integralG(f, a, b, n)
    steps = 1
    while not runge_check(I2n, In, GAUSS_OMEGA, E):
        In = I2n
        n *= 2
        I2n = sum_integralG(f, a, b, n)
        steps += 1    
    return I2n, steps

result, steps = calculate_gauss(F, 0, pi/2, E)
print("gauss_method:", result, "| steps:", steps*2)
print()
result, steps = calculate(F, 0, pi/2, E, RECTANGLES_OMEGA, rectangles_method)
print("rectangles_method:", result, "| steps:", steps)
print()
result, steps = calculate(F, 0, pi/2, E, TRAPEZOIDAL_OMEGA, trapezoidal_method)
print("trapezoidal_method:", result, "| steps:", steps)
print()
result, steps = calculate(F, 0, pi/2, E, SIMPSON_OMEGA, simpson_method)
print("simpson_method:", result, "| steps:", steps)