# About PV_TSSOS
PV_TSSOS is a Julia package for both polynomial optimization and polynomial systems by combining [TSSOS](https://github.com/wangjie212/TSSOS) and the SDP hierarchy relying on [Putinar-Vasilescu's Positivstellensatz](https://arxiv.org/abs/1911.11428). 
## New advances
- Obtaining an approximate global optimizer of a general polynomial optimization 

     ![POP](https://github.com/maihoanganh/PV_TSSOS/blob/master/images/POP.gif)

with guarantee in theory.
- Obtaining generically a feasible solution of a basic semialgebraic set 

    ![POP](https://github.com/maihoanganh/PV_TSSOS/blob/master/images/S.gif)

with possibly uncountably many solutions.
## Installation
To use PV_TSSOS in Julia, run:
```ruby
pkg> add https://github.com/maihoanganh/PV_TSSOS
```
# Dependency
MOSEK (SDP solver)

# Usage
The following examples are performed on WINDOW 10, Julia 1.1.1, and MOSEK 8.0.
## 1. Polynomial optimization
### 1.1. Compact case
```ruby
using DynamicPolynomials
using PV_TSSOS

k=2 # relaxed order
r=1 # sparse order

@polyvar x1 x2 # define polnomial variables
x=[x1;x2] 

f = x1^4+x2^4-x1*x2 # objective function
g = [x1;1-x1^2-x2^2] # inequality constraints
h = [x1+1/2-(x2+1/2)^3] # equality constraints

PV_TSSOS.block_compact_POP(x,f,g,h,k,r)
```
Feedback: ```opt_val=-0.125``` and ```sol=[0.499998, 0.500002]```.

### 1.2. Non-compact case
#### 1.2.1. Unconstrained case
```ruby
using DynamicPolynomials
using PV_TSSOS

eps = 1e-5 # small parameter
k=2 # relaxed order
r=1 # sparse order

@polyvar x1 x2 # define polnomial variables
x=[x1;x2] 

f = x1^2*x2^2*(x1^2+x2^2-1) # objective function
g = [] # inequality constraints
h = [] # equality constraints

opt_val=PV_TSSOS.block_noncompact_POP(x,f,g,h,eps,k,r) # take the optimal value

g=[opt_val-f]

k=7 # relaxed order
r=2 # sparse order

PV_TSSOS.adding_spherical_constraints(x,g,h,k,r) # compute minimizers
```
Feedback: ```opt_val=-0.0369``` and ```sol=[0.571385, 0.571382]```.

#### 1.2.2.  Constrained case
```ruby
using DynamicPolynomials
using PV_TSSOS

eps = 1e-5 # small parameter
k=5 # relaxed order
r=2 # sparse order

@polyvar x1 x2 x3 # define polnomial variables
x=[x1;x2;x3] 

f = x1+x2+x3 # objective function
g = [x1;x2;x3] # inequality constraints
h = [x1*x2*x3-1] # equality constraints

opt_val=PV_TSSOS.block_noncompact_POP(x,f,g,h,eps,k,r) # take the optimal value

g=[g;opt_val-f]

k=5 # relaxed order
r=2 # sparse order

PV_TSSOS.adding_spherical_constraints(x,g,h,k,r) # compute minimizers
```
Feedback: ```opt_val=3.000``` and ```sol=[1.00043, 0.999788, 0.999782]```.

## 2. Polynomial systems
```ruby
using DynamicPolynomials
using PV_TSSOS

k=5 # relaxed order
r=2 # sparse order

@polyvar x1 x2 x3 x4 # define polnomial variables
x=[x1;x2;x3;x4] 

g = [] #inequality constraints
h = [x1*x2^2 + x1*x3^2 - 1.1*x1 + 1;
    x2*x1^2 + x2*x3^2 - 1.1*x2 + 1;
    x3*x1^2 + x3*x2^2 - 1.1*x3 + 1] #equality constraints

PV_TSSOS.adding_spherical_constraints(x,g,h,k,r)
```
Feedback: ```sol=[-1.01992, -1.01992, -1.01992, 5.59501e-8]```.

# References
For more details, please refer to:
1. N. H. A. Mai, J.-B. Lasserre, and V. Magron. Positivity certificates and polynomial optimization on non-compact semialgebraic sets, 2019. Submitted.
https://arxiv.org/abs/1911.11428
2. J. Wang, V. Magron, and J.-B. Lasserre. TSSOS: a moment-SOS hierarchy that exploits term sparsity, 2019. Submitted. 
https://arxiv.org/abs/1912.08899
3. N. H. A. Mai, J. Wang, V. Magron, and J.-B. Lasserre. Exploiting term sparsity for polynomial optimization on non-compact semialgebraic sets, 2020. Forthcoming.
