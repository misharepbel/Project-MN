from math import prod
from matplotlib import pyplot as plt
import numpy as np

# Functions to generate variables for the functions or the nodes for the interpolation
def generate_xes(a, b, n):
    #vars for the function graph
    x=np.linspace(a,b,1000)
    #classic interpolation nodes
    x_int=(np.linspace(a,b,n+1))
    #Chebyshev nodes
    x_ch=[((b-a)/2)*np.cos(np.pi*(2*i+1)/(2*n+2))+((a+b)/2) for i in range(n+1)]
    res = { "function": x, "interpolation": x_int, "chebyshev": x_ch }
    return res

# The function we are approximating via interpolation
def test_function(x, f):
    return([f(i) for i in x])

# Newton interpolation function
def interpolate(x, y):
    newton_table = [[y[len(y)-i]]+[0]*(i-1) for i in range(len(y), 0, -1)]
    for i in range(len(y)-2, -1, -1):
        for j in range(1, len(newton_table[i])):
            newton_table[i][j]=(newton_table[i+1][j-1]-newton_table[i][j-1])/(x[i+j]-x[i])
    return [lambda z: newton_table[0][0]+sum([newton_table[0][i]*prod([(z-x[j]) for j in range(i)]) for i in range(1, len(newton_table[0]))]), newton_table[0]]

# Function for the interpolation error
def compute_error(x_test, f_true, f_int):
    f_true = test_function(x_test, f_true)
    f_int = test_function(x_test, f_int)
    return (np.abs(np.array(f_true) - np.array(f_int)))

# main 

fout = open("nodes.txt", "w")
fout2 = open("errors.txt", "w")
fout3 = open("x_and_ys.txt", "w")
fields = generate_xes(-5.3, 4.7, 10)
fout.write("Different interpolation points:\n")
fout.write("interpolation:\n"+str(fields["interpolation"]))
fout.write("\nChebyshev:\n"+str(fields["chebyshev"]))
fu = lambda x: 25/(1+3*x**2)
fsi = lambda x: np.sin(1/(x**2+np.exp(x)))

interpol = interpolate(fields["interpolation"], test_function(fields["interpolation"], fsi))
cheburashka = interpolate(fields["chebyshev"], test_function(fields["chebyshev"], fsi))

fout3.write(f'Int:\nx:\n{fields["interpolation"]}\ny:\n{interpol[1]}\n')
fout3.write(f'Cheb:\nx:\n{fields["chebyshev"]}\ny:\n{cheburashka[1]}\n')

err_int = compute_error(fields["interpolation"], fsi, interpol[0])
err_cheb = compute_error(fields["chebyshev"], fsi, cheburashka[0])

fout2.write("Interpolation:\n")
fout2.write(str(err_int))
fout2.write(f"\nsum: {sum(err_int)}")
fout2.write("\nChebyshev:\n")
fout2.write(str(err_cheb))
fout2.write(f"\nsum: {sum(err_cheb)}")

plt.plot(fields["function"], test_function(fields["function"], fsi))
plt.plot(fields["function"], test_function(fields["function"], interpol[0]), color="red")
plt.plot(fields["function"], test_function(fields["function"], cheburashka[0]), color="black")
plt.gca().spines['left'].set_position(('zero'))
plt.gca().spines['right'].set_color('none')
plt.legend(["Original function", "Equidistant nodes", "Chebyshev nodes"])
plt.grid()
plt.show()