# Python library that calculates metric related things on Manifolds.

It calculates:
- Chrystoffel Symbols
- Parallel transport equations
- Metric coefficients


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

## TO-DO
-Add examples
-Add parallel equations solver (I actually have this somewhere)
