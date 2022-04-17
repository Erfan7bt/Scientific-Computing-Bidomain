### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ bb68c5f8-d906-474b-8492-bec212fb5def
using ShortCodes, LaTeXStrings

# ╔═╡ 3626452a-23be-44de-be1e-a0ebefd7d170
md"""
# Bidomain Model

The bidomain model is a system of partial differential equations used to model the propagation of electrical potential waves in the myocardium[1]. It is formed of an intracellular potential ``u_i`` and an extracellular potential ``u_e``.
"""

# ╔═╡ 6c937839-6ac4-4dec-878b-8a18d0f5d6bd
md"""
## Bidomain Problem Modeled as a system of PDEs


$u= u_i-u_e, \ u_e \; and \; \ v \;on \: \ Q\times\Omega= [0,\ T=30] \times [0, \ L=70]$

$\frac{\partial u}{\partial t} -
	\nabla\cdot(\sigma_i(\nabla u+ \nabla u_e))
     - \frac{1}{\epsilon} f(u,v) = 0$

$\nabla\cdot(\sigma_i\nabla u+(\sigma_i +\sigma_e)\nabla u_e)=0$

$\frac{\partial v}{\partial t} - \epsilon g(u,v)=0$

Where ``\sigma_i`` and ``\sigma_e`` are second order tensors representing the intracellular and extracellular tissue's electrical conductivity in each spatial direction, ``v`` is a lumped ionic variable, and ``\eps`` is a parameter linked to the ratio between the repolarisation rate and the tissue excitation rate[1].

We take the initial and boundary conditions:

 $u$ and $v$ constant at the equilibrium value of the system $f(u,v)=g(u,v)=0$ and $u_e$ constant at $0$ , except on the interval $[0,L/20]$ where we fix $u$ at a supercritical value, $u(x)=2$  

So we write :  

$u_e(0,x)=0 \:\; \forall x\in \Omega$  


$u(0,x)=2 \;\;\forall x \in [0,L/20]$  

$f(u(0,x_1),v(0,x))=g(u(0,x_1),v(0,x))=0 \;\;\forall x_1 \in [L/20,L] \;\; and\;\; \forall x \in [0,L]$

Where $f(u,v) = u - \frac{u^3}{3} - v$ and $g(u,v) = u + \beta -\gamma v$
and Neumann boundary conditions:

$\partial_x u(0)=\partial_x u(L)=0 \;\; and \;\; \partial_x u_e(0)=\partial_x u_e(L)=0$
Since pure Neumann boundary conditions are ill-posed, to avoid solving a singular system we set:
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

 storage, recation and flux terms in our system are not dependent on $\vec{x}$ and $t$, but in general they can be.
"""

# ╔═╡ 90479d25-4018-4074-8e31-fe2ef9d881eb
md"""
### Discretization over Finite Volumes $\omega_k$:
#### Constructing control volumes
We construct and formulate exactly as it was done in the lecture 
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
By use of Gauss's thereom, the integral of divergence of flux over volume become an integral of flux multiplied by the normal vector over the boundary:
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

After following same procedure as the discretization of Equation 1 and since there is no storage term in this equation:
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

# ╔═╡ 35973cf2-18b6-4e1d-b5a7-70bd5cb85569
md"""
### Observations

We successfully recreated the 1d problem without much difficulty, but had a lot of trouble recreating the spiral pattern observed in [1]. It turns out that the solution is very sensitive to spatial step-size-- less than roughly 100 cells per L=70 in each dimension would result in an incorrect solution that did not display conductivity in that direction. To avoid this, we ran our spiral solution overnight with N=120x120 and ``\Delta t = 10^{-2}``.
"""

# ╔═╡ ead7c30f-72c4-48bf-8e13-95b4caa0ccd1
md"""
# Bibliography
"""

# ╔═╡ 2374c5d0-bdad-11ec-1718-05131c0e0731
begin
	references=[
	DOI("10.1137/070680503")
]
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
ShortCodes = "f62ebe17-55c5-4640-972f-b59c0dd11ccf"

[compat]
LaTeXStrings = "~1.3.0"
ShortCodes = "~0.3.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON3]]
deps = ["Dates", "Mmap", "Parsers", "StructTypes", "UUIDs"]
git-tree-sha1 = "8c1f668b24d999fb47baf80436194fdccec65ad2"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.9.4"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Memoize]]
deps = ["MacroTools"]
git-tree-sha1 = "2b1dfcba103de714d31c033b5dacc2e4a12c7caa"
uuid = "c03570c3-d221-55d1-a50c-7939bbd78826"
version = "0.4.4"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "621f4f3b4977325b9128d5fae7a8b4829a0c2222"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.4"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[ShortCodes]]
deps = ["Base64", "CodecZlib", "HTTP", "JSON3", "Memoize", "UUIDs"]
git-tree-sha1 = "0fcc38215160e0a964e9b0f0c25dcca3b2112ad1"
uuid = "f62ebe17-55c5-4640-972f-b59c0dd11ccf"
version = "0.3.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "d24a825a95a6d98c385001212dc9020d609f2d4f"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.8.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
"""

# ╔═╡ Cell order:
# ╠═bb68c5f8-d906-474b-8492-bec212fb5def
# ╠═3626452a-23be-44de-be1e-a0ebefd7d170
# ╠═6c937839-6ac4-4dec-878b-8a18d0f5d6bd
# ╟─c7aaa1d6-e85a-4075-9191-a2cca74f163c
# ╟─90479d25-4018-4074-8e31-fe2ef9d881eb
# ╟─82660825-5862-49c7-be5f-d034b1e684f5
# ╟─eba8c349-7097-401d-bb86-5c5ec529f197
# ╟─b9d4522b-a589-4a96-befa-ea4ee5ab0fe2
# ╟─af8bc723-a46a-4621-af80-d3a319da173b
# ╠═198ec4a6-27dc-41ff-be30-d4e1fe6c35d3
# ╠═35973cf2-18b6-4e1d-b5a7-70bd5cb85569
# ╟─ead7c30f-72c4-48bf-8e13-95b4caa0ccd1
# ╠═2374c5d0-bdad-11ec-1718-05131c0e0731
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
