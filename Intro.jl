### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ 6c937839-6ac4-4dec-878b-8a18d0f5d6bd
md"""
# Biodomain Problem Modeled as a syetem of PDEs


$u= u_i-u_e, \ u_e \; and \; \ v \;on \: \ Q\times\Omega= [0,\ T=30] \times [0, \ L=70]$

$\frac{\partial u}{\partial t} -
	\nabla\cdot(\sigma_i(\nabla u+ \nabla u_e))
     - \frac{1}{\epsilon} f(u,v) = 0$

$\nabla\cdot(\sigma_i\nabla u+(\sigma_i +\sigma_e)\nabla u_e)=0$

$\frac{\partial v}{\partial t} - \epsilon g(u,v)=0$.

with intial and boundary condition:

 $u$ and $v$ constatnt at the equiblrium value of the system $f(u,v)=g(u,v)=0$ and $u_e$ constant at $0$ , expect on the interval $[0,L/20]$ where we fix $u$ at a supercritical value, $u(x)=2$  

So we write :  

$u_e(0,x)=0 \:\; \forall x\in \Omega$  


$u(0,x)=2 \;\;\forall x \in [0,L/20]$  

$f(u(0,x_1),v(0,x))=g(u(0,x_1),v(0,x))=0 \;\;\forall x_1 \in [L/20,L] \;\; and\;\; \forall x \in [0,L]$

Where $f(u,v) = u - \frac{u^3}{3} - v$ and $g(u,v) = u + \beta -\gamma v$
and Neumann boundary conditions:

$\partial_x u(0)=\partial_x u(L)=0 \;\; and \;\; \partial_x u_e(0)=\partial_x u_e(L)=0$
since pure neumann boundary condition is ill-posed  to avoide solving a singular system:
$u_e(0)=0.$
"""

# ╔═╡ c7aaa1d6-e85a-4075-9191-a2cca74f163c
md"""
### Bidomain as a sysytem of reaction-diffusion system
In order to define bidomain problem as a system of $n$ coupled PDEs so that we can handle it to VoronoiFVM we define "**reaction**" $r(u)$, ""**storage"**" $s(u)$, ""**source"**" $f$ and ""**flux"**" $\vec{j}$ in vector form as follows:  

 $\partial_t s(u) + \nabla \cdot \vec{j}(u) + r(u) = f$  

 $u = \begin{bmatrix} u\\u_e\\v\end{bmatrix}$  

 $s(u) = \begin{bmatrix} u\\0\\v\end{bmatrix}$   

 $\vec{j}(u) = \begin{bmatrix} -\sigma_i(\nabla u+\nabla u_e)\\
\sigma_i\nabla u + (\sigma_i+\sigma_e)\nabla u_e
\\0\end{bmatrix}$  

 $r(u) = \begin{bmatrix} -f(u,v)/\epsilon\\0\\ -\epsilon g(u,v)\end{bmatrix}$  

 $f = 0$.  

 stoarge, recation and flux terms in our system is not dependent on $\vec{x}$ and $t$, but in general they can be.
"""

# ╔═╡ 90479d25-4018-4074-8e31-fe2ef9d881eb
md"""
### Discretization over Finite Volumes $\omega_k$:
#### Constructing control volumes
We construct and formulate excatly as it was done in the Lecture 
"""

# ╔═╡ 82660825-5862-49c7-be5f-d034b1e684f5
md"""
  Assume $\Omega\subset \mathbb{R}^d$ is a polygonal domain such that
  $\partial\Omega=\bigcup_{m\in\mathcal{G}} \Gamma_m$, where $\Gamma_m$ are planar such that $\vec{n}|_{\Gamma_m}=\vec{n}_m$.

  Subdivide $\Omega$ into into a finite number of
  _control volumes_ 
  $\bar\Omega= \bigcup_{k\in \mathcal N} \bar \omega_k$
  such that 
  
  -  $\omega_k$ are open  convex domains such that  $\omega_k\cap \omega_l = \emptyset$ if
    $\omega_k\neq \omega_l$ 
  -  $\sigma_{kl}=\bar \omega_k \cap \bar \omega_l$ are either empty,
    points or straight lines.
    If $|\sigma_{kl}|>0$ we say that $\omega_k$, $\omega_l$
    are neighbours.
  -  $\vec{n}_{kl}\perp \sigma_{kl}$: normal of $\partial\omega_k$ at $\sigma_{kl}$
  -  $\mathcal{N}_k = \{l \in \mathcal{N}: |\sigma_{kl}|>0\}$: set of neighbours of $\omega_k$
  -  $\gamma_{km}=\partial\omega_k \cap \Gamma_m$: boundary part of $\partial\omega_k$
  -  $\mathcal{G}_k= \{ m \in \mathcal{G}:  |\gamma_{km}|>0\}$: set of non-empty boundary parts of $\partial\omega_k$.

  
  $\Rightarrow$ $\partial\omega_k= \left(\cup_{l\in \mathcal{N}_k} \sigma_{kl}\right)\bigcup\left(\cup_{m\in \mathcal{G}_k} \gamma_{km}\right)$

"""

# ╔═╡ eba8c349-7097-401d-bb86-5c5ec529f197
md"""
  To each control volume $\omega_k$ assign a _collocation
    point_: $\vec{x}_k \in \bar \omega_k$ such that\\
  
  -  _Admissibility condition_:if $l\in \mathcal N_k$ then the
    line $\vec{x}_k\vec{x}_l$ is orthogonal to $\sigma_{kl}$
    
    -  For a given function $u:\Omega \to \mathbb{R}$ this will allow to associate its value $u_k=u(\vec{x}_k)$
      as the value of an unknown at $\vec{x}_k$.
    -  For two neigboring control volumes $\omega_k, \omega_l$ , this will allow to approximate
      $\vec\nabla u \cdot \vec{n}_{kl} \approx  \frac{u_l - u_k}{h_{kl}}$
    
    
  -  _Placement of boundary unknowns at the boundary_: if
    $\omega_k$ is situated at the boundary, i.e. for 
    $|\partial \omega_k \cap \partial \Omega| >0$,
    then $\vec{x}_k \in \partial\Omega$
    
    -  This will allow to apply boundary conditions in a direct manner

"""

# ╔═╡ b9d4522b-a589-4a96-befa-ea4ee5ab0fe2
md"""
### First Equation 
For $k\in\mathcal{N}$ Integrate over  each control volume $\omega_k$:
```math
\newcommand{\eps}{\varepsilon}
\newcommand{\pth}[1]{\left(#1\right)}


\begin{equation}
   \int_{\omega_k}\frac{\partial u}{\partial t}=\int_{\omega_k}\frac{1}{\eps}f(u,v)d\omega +\\
\int_{\omega_k}\nabla \cdot (\sigma_i \nabla u)d\omega + \int_{\omega_k}\nabla \cdot (\sigma_i \nabla u_e)d\omega
\end{equation}
```
By use of Gauss's thereom so the integral of divergence of flux over volume become integral of flux multiply by normal vector over bundary:
```math
\begin{align*}
    \int_{\omega_k}\frac{\partial u}{\partial t} =\int_{\omega_k}\frac{1}{\eps}f(u,v)d\omega +
    \int_{\partial\omega_k}\sigma_i \nabla u \cdot \vec{n}ds +
    \int_{\partial\omega_k}\sigma_i \nabla u_e \cdot \vec{n}ds
\end{align*}
```
Since the $\partial\Omega_k$ is either edges in domain boundary $\mathcal{G}_k$ or neighbours $\mathcal{N}_k$ 
```math
\begin{align*}
    \int_{\omega_k}\frac{\partial u}{\partial t} &=\int_{\omega_k}\frac{1}{\eps}f(u,v)d\omega + \sum_{l \in N_k}\int_{\sigma_{kl}} \sigma_i \nabla u \cdot \vec{n}_{kl}ds +
    \sum_{m \in \mathcal{G}_k}\int_{\gamma_{kl}} \sigma_i \nabla u \cdot \vec{n}_{m}ds \\
    &+ \sum_{l \in N_k}\int_{\sigma_{kl}} \sigma_i \nabla u_e \cdot \vec{n}_{kl}ds +
    \sum_{m \in \mathcal{G}_k}\int_{\gamma_{kl}} \sigma_i \nabla u_e \cdot \vec{n}_{m}ds
\end{align*}
```
Finite difference approximation of normal derivative:
```math
\begin{align*}
    h_{kl} = |x_k - x_l|\\
\sigma_i \nabla u \cdot \vec{n} \approx \sigma_i \frac{u_k - u_l}{h_{kl}}\\
    \sigma_i \nabla u_e \cdot \vec{n} \approx \sigma_i \frac{u_{e_k} - u_{e_l}}{h_{kl}}
\end{align*}
```
We also approximate integrals by the length of the edge times approximated flux: 
```math
\begin{align*}
|\omega_k|\frac{\partial u_k}{\partial t} &= \int_{w_k}\frac{1}{\eps}(u - \frac{u^3}{3} - v)d\omega + \sum_{l \in N_k} |\sigma_{kl}| \sigma_i  \frac{u_k - u_l}{h_{kl}}  + 
\sum_{m \in \mathcal{G}_k} |\gamma_{km}| \sigma_i  \frac{u_k - u_l}{h_{kl}}   \\
&+ \sum_{l \in N_k} |\sigma_{kl}| \sigma_i  \frac{u_{e_k} - u_{e_l}}{h_{kl}}  + 
\sum_{m \in \mathcal{G}_k} |\gamma_{km}| \sigma_i  \frac{u_{e_k} - u_{e_l}}{h_{kl}}  
\end{align*}
```
 Approximate reaction integral by multiplying $|\omega_k|$ and combining $u$ and $u_e$ sum:
```math
\begin{align*}
|\omega_k|\frac{\partial u_k}{\partial t} &= \frac{|\omega_k|}{\eps}\pth{u_k - \frac{u_k^3}{3} - v_k}+ \sum_{l \in N_k} |\sigma_{kl}| \sigma_i  \frac{u_k - u_l + u_{e_k} - u_{e_l}}{h_{kl}} \\
&+ \sum_{m \in \mathcal{G}_k} |\gamma_{km}| \sigma_i  \frac{u_k - u_l + u_{e_k} - u_{e_l}}{h_{kl}}  
\end{align*}
```

Then discretizing in time in forward Euler method yields:
```math
\begin{align}
\frac{|\omega_k|}{\Delta t}(u_k^{n+1} - u_k^n) = \sum_{l \in N_k} |\sigma_{kl}| \sigma_i \frac{u_k^{n} - u_l^{n} + u_{e_k}^n - u_{e_l}^{n}}{h_{kl}} 
&+ \frac{|\omega_k|}{\eps}\pth{u_k^{n} - \frac{(u_k^{n})^3}{3} - v_k^{n}}
\\+
\sum_{m \in \mathcal{G}_k} |\gamma_{km}| \sigma_i  \frac{u_k^{n} - u_l^{n} + u_{e_k}^{n} - u_{e_l}^{n}}{h_{kl}} 
\end{align}
```
As $u_k^n = u(\vec{x_k},n\Delta t)$
"""

# ╔═╡ af8bc723-a46a-4621-af80-d3a319da173b
md"""
#### Second Equation 2 

After Folloeing same procedure as Dicretization of Equation 1:
 Since there is no stoarge term in this equation
  $u_k=u_k^n = u(\vec{x_k},n\Delta t)$
```math
\begin{align}
\sum_{l \in N_k} |\sigma_{kl}| \sigma_i \frac{u_k - u_l}{h_{kl}}
+ \sum_{l \in N_k} |\sigma_{kl}| \pth{\sigma_i+\sigma_e}  \frac{u_{e_k} - u_{e_l}}{h_{kl}}  \\
+\sum_{m \in \mathcal{G}_k} |\gamma_{km}| \sigma_i \frac{u_k - u_l}{h_{kl}}  
+\sum_{m \in \mathcal{G}_k} |\gamma_{km}| (\sigma_i+\sigma_e)  \frac{u_{e_k} - u_{e_l}}{h_{kl}}=0
\end{align}
```
"""

# ╔═╡ 198ec4a6-27dc-41ff-be30-d4e1fe6c35d3
md"""
#### Discretization of equation 3
With 
```math
\begin{align*}
\frac{\partial v}{\partial t} - \eps(u + \beta - \varphi v) = 0
\end{align*}
```
We take the integral as before with respect to the volume $\omega_k$:
```math
\begin{align*}
\int_{\omega_k}\frac{\partial v}{\partial t} d\omega - \int_{\omega_k} \eps(u + \beta - \varphi v)d\omega = 0
\end{align*}
```
The integral is simply a multiplication by the area of the volume:
```math
\begin{align*}
|\omega_k|\frac{\partial v_k}{\partial t} - \eps |\omega_k| (u_k + \beta - \varphi v_k) = 0
\end{align*}
```
Then discretizing in time yields:
```math
\begin{align*}
\frac{|\omega_k|}{\Delta t}(v_k^n - v_k^{n-1}) - \eps |\omega_k| (u_k^{\theta} + \beta - \varphi v_k^{\theta}) = 0
\end{align*}
```

This concludes the space discretization.
"""

# ╔═╡ 2374c5d0-bdad-11ec-1718-05131c0e0731


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─6c937839-6ac4-4dec-878b-8a18d0f5d6bd
# ╟─c7aaa1d6-e85a-4075-9191-a2cca74f163c
# ╟─90479d25-4018-4074-8e31-fe2ef9d881eb
# ╟─82660825-5862-49c7-be5f-d034b1e684f5
# ╟─eba8c349-7097-401d-bb86-5c5ec529f197
# ╟─b9d4522b-a589-4a96-befa-ea4ee5ab0fe2
# ╟─af8bc723-a46a-4621-af80-d3a319da173b
# ╟─198ec4a6-27dc-41ff-be30-d4e1fe6c35d3
# ╠═2374c5d0-bdad-11ec-1718-05131c0e0731
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
