### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 559db84a-26b5-4ae3-883a-468c5e27a7d2
begin
	using PlutoUI;
	TableOfContents();
end

# ╔═╡ 2374c5d0-bdad-11ec-1718-05131c0e0731
begin
	using ShortCodes
	references=[
		DOI("10.1137/070680503"),
		DOI("10.5281/zenodo.3529808"),
		DOI("10.21105/joss.01848")
	]
end

# ╔═╡ 063322b5-fcb1-48c9-9fd2-8a92101bb44f
md"""

# Bidomain Model
## Introduction

This is the report of Implementation of Bidomain Model as Formulated in [1], Usng **VoronoiFVM.jl** [2], as the Final project of the course "Scienfic Computing" at Techinical University Of Belin in winter semster of 21/22, under Supervision of Dr. Juergen Fuhrmann.  

Presented and implemented by Alexander Quinlan and Erfan Baradarantohidi

## Bidomain Model 
The bidomain  is a system of partial differential equations used to model the propagation of electrical potential waves in myocardium. It is composed of coupled parabolic and eliptic PDEs, as well as at least one ordinary differential equation to model the ion activity through the cardic cell membrance [1].  

The electrical properties of the myocardium are generally described by the bidomain equations, a set of coupled parabolic and elliptic partial differential equations (PDEs) that represents the tissue as two separate, distinct continua - one intracellular and the other extracellular. The intracellular and the extracellular media are connected via the cell membrane, and thus the two PDEs are coupled at each point in space through a set of complex, non-linear ordinary differential equations (ODEs), which describe the ionic transport across the cell membrane. Certain modelling environments use the monodomain representation of cardiac activity, which involves solving a single parabolic PDE, by assuming either that the extracellular potentials are negligible, or that the anisotropy ratios are equal in the intracellular and the extracellular domains[3].
 
"""

# ╔═╡ 6c937839-6ac4-4dec-878b-8a18d0f5d6bd
md"""
# Biodomain Problem Modeled as a system of PDEs


$u= u_i-u_e, \ u_e \; and \; \ v \;on \: \ Q\times\Omega= [0,\ T=30] \times [0, \ L=70]$

$\frac{\partial u}{\partial t} -
	\nabla\cdot(\sigma_i(\nabla u+ \nabla u_e))
     - \frac{1}{\epsilon} f(u,v) = 0$

$\nabla\cdot(\sigma_i\nabla u+(\sigma_i +\sigma_e)\nabla u_e)=0$

$\frac{\partial v}{\partial t} - \epsilon g(u,v)=0$

We take the initial and boundary conditions $u$ and $v$ constant at the equilibrium value of the system $f(u,v)=g(u,v)=0$ and $u_e$ constant at $0$ , except on the interval $[0,L/20]$ where we fix $u$ at a supercritical value, $u(x)=2$  

So we write :  

$u_e(0,x)=0 \:\; \forall x\in \Omega$  


$u(0,x)=2 \;\;\forall x \in [0,L/20]$  

$f(u(0,x_1),v(0,x))=g(u(0,x_1),v(0,x))=0 \;\;\forall x_1 \in [L/20,L] \;\; and\;\; \forall x \in [0,L]$

Where $f(u,v) = u - \frac{u^3}{3} - v$ and $g(u,v) = u + \beta -\gamma v$
and Neumann boundary conditions:

$\partial_x u(0)=\partial_x u(L)=0 \;\; and \;\; \partial_x u_e(0)=\partial_x u_e(L)=0$
Since pure Neumann boundary conditions are ill-posed, we avoid solving a singular system with:
$u_e(0)=0.$
"""

# ╔═╡ b7926230-80ed-4a06-9a45-f5191ccf2aef
md"""
### What do these variables describe? 
At each point in the computational domain, two electrical potential:
- ``u_i`` the **interacellular** potentionial
- ``u_e`` the **extracellular** potentionial
are recovered, representing the average of the electrical potential over the extracellular and the interacellular space, recpectively, in the neighborhoods of that point. 
- ``u=u_i - u_e`` is the **transmembrance** potential
- ``\sigma_i`` and ``\sigma_e`` are second order tensor, representing **electrical conductivity** in each spatial direcation (in this implementation, it is assumed that conductivity of heart tissue is same in each direction, hence it is constant)
- ``v`` is a lumped **ionic variable** (the ODE as last equation is not properly part of the bidomain model but rather models the ionic activity across the cellular membrance which is responsible foe electrical activation of the tissue)
- ``\epsilon`` is a paarmeter linled to the **ratio between the repolarization rate and the tissue excitation rate**
- function ``f`` and ``g`` :

$f(u,v) = u - \frac{u^3}{3} - v$
$g(u,v) = u + \beta -\gamma v$

- ``\gamma `` control **ion transport** , ``\beta`` is linked to **cell excitabillity**

"""

# ╔═╡ c7aaa1d6-e85a-4075-9191-a2cca74f163c
md"""
### Bidomain as a reaction-diffusion system
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
# Discretization over Finite Volumes $\omega_k$:
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
### Second Equation  

After following same procedure as the discretization of Equation 1:
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
### Last Equation 

```math
\begin{align*}
\int_{\omega_k}\frac{\partial v}{\partial t} d\omega = \int_{\omega_k} \eps(u + \beta - \gamma v)d\omega 
\end{align*}
```
Approximation of Integrals:
```math
\begin{align*}
|\omega_k|\frac{\partial v_k}{\partial t} = \eps |\omega_k| (u_k + \beta - \gamma v_k) 
\end{align*}
```
Then forward discretizing in time :
```math
\begin{align*}
\frac{|\omega_k|}{\Delta t}(v_k^{n+1} - v_k^{n}) - \eps |\omega_k| (u_k^{n} + \beta - \gamma v_k^{n})
\end{align*}
```
"""

# ╔═╡ 5fb78a95-d6bb-45f8-800b-486713d4851a
md"""
### Observations

We successfully recreated the 1d problem without much difficulty, but had a lot of trouble recreating the spiral pattern observed in [1]. It turns out that the solution is very sensitive to spatial step-size-- less than roughly 100 cells per L=70 in each dimension would result in an incorrect solution that did not display conductivity in that direction. To avoid this, we ran our spiral solution overnight with N=120x120 and ``\Delta t = 10^{-2}``.
"""

# ╔═╡ 0a809ed0-9aec-4238-8b21-49eeb6d9c0fb
md"""
# Bibliography
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ShortCodes = "f62ebe17-55c5-4640-972f-b59c0dd11ccf"

[compat]
PlutoUI = "~0.7.38"
ShortCodes = "~0.3.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.6.3"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "StructTypes", "UUIDs"]
git-tree-sha1 = "8c1f668b24d999fb47baf80436194fdccec65ad2"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.9.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Memoize]]
deps = ["MacroTools"]
git-tree-sha1 = "2b1dfcba103de714d31c033b5dacc2e4a12c7caa"
uuid = "c03570c3-d221-55d1-a50c-7939bbd78826"
version = "0.4.4"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "621f4f3b4977325b9128d5fae7a8b4829a0c2222"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.ShortCodes]]
deps = ["Base64", "CodecZlib", "HTTP", "JSON3", "Memoize", "UUIDs"]
git-tree-sha1 = "0fcc38215160e0a964e9b0f0c25dcca3b2112ad1"
uuid = "f62ebe17-55c5-4640-972f-b59c0dd11ccf"
version = "0.3.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "d24a825a95a6d98c385001212dc9020d609f2d4f"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.8.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─559db84a-26b5-4ae3-883a-468c5e27a7d2
# ╟─063322b5-fcb1-48c9-9fd2-8a92101bb44f
# ╟─6c937839-6ac4-4dec-878b-8a18d0f5d6bd
# ╟─b7926230-80ed-4a06-9a45-f5191ccf2aef
# ╟─c7aaa1d6-e85a-4075-9191-a2cca74f163c
# ╟─90479d25-4018-4074-8e31-fe2ef9d881eb
# ╟─82660825-5862-49c7-be5f-d034b1e684f5
# ╟─eba8c349-7097-401d-bb86-5c5ec529f197
# ╟─b9d4522b-a589-4a96-befa-ea4ee5ab0fe2
# ╟─af8bc723-a46a-4621-af80-d3a319da173b
# ╟─198ec4a6-27dc-41ff-be30-d4e1fe6c35d3
# ╟─5fb78a95-d6bb-45f8-800b-486713d4851a
# ╟─0a809ed0-9aec-4238-8b21-49eeb6d9c0fb
# ╟─2374c5d0-bdad-11ec-1718-05131c0e0731
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
