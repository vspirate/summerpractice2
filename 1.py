def calculate(a):
    k = 1
    pow_of_two = 2
    fac = 1
    even_or_odd = -1
    o = even_or_odd/(pow_of_two*fac)
    sum = o
    while abs(o) > a:
        k += 1
        fac *= k
        pow_of_two *= 2
        even_or_odd *= -1
        o = even_or_odd/(pow_of_two*fac)
        sum += o
    
    return sum

print(calculate(0.001))