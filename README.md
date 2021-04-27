# Python library that calculates metric related things on Manifolds.

It calculates:
- Chrystoffel Symbols
- Parallel transport equations
- Metric coefficients

## Latex

It is also possible to get your result in latex, however
to avoid messy logic with "\\" in python, 
those have being replaced by "~Ñ".

You will have to replace them afterwards in your tex document.


## Example

```
myLobat = Lobatchevsky()
print(myLobat)
```

Will produce:
```
Manifold name: Lobatchevsky Plane


Coordinates
x_1 : x1
x_2 : x2

Metric
g_11 : y**(-2)
g_12 : 0
g_21 : 0
g_22 : y**(-2)

Chrystoffel Symbols
Gamma_11^1 : 0
Gamma_12^1 : -1.0/y
Gamma_21^1 : -1.0/y
Gamma_22^1 : 0
Gamma_11^2 : 1.0/y
Gamma_12^2 : 0
Gamma_21^2 : 0
Gamma_22^2 : -1.0/y

Geodesics equations
ddx1 - 2.0*dx1*dx2/x2 
(ddx2*x2 + 1.0*dx1**2 - 1.0*dx2**2)/x2 
```

We can also calculate the symbols of a sphere:
```
Sphere = S2(1)
print(Sphere)
```

```
Manifold name: 2-Sphere


Coordinates
x_1 : x1
x_2 : x2

Metric
g_11 : 1
g_12 : 0
g_21 : 0
g_22 : cos(theta)**2

Chrystoffel Symbols
Gamma_11^1 : 0
Gamma_12^1 : 0
Gamma_21^1 : 0
Gamma_22^1 : 0.5*sin(2*theta)
Gamma_11^2 : 0
Gamma_12^2 : -1.0*tan(theta)
Gamma_21^2 : -1.0*tan(theta)
Gamma_22^2 : 0

Geodesics equations
ddx1 + 0.5*dx2**2*sin(2*x1) 
ddx2 - 2.0*dx1*dx2*tan(x1) 
```

## TO-DO
- Add examples
- Add parallel equations solver (I actually have this somewhere)
- Make it more modular: Separate classes into their own files
