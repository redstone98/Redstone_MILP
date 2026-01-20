# Redstone_MILP
Redstone's LEO satellite constellation optimization using Mixed-Integer Linear Programming

# Satellite Constellation Scheduling to MILP/MIQP Problem

This note summarizes how an unconstrained SAT–GS contact chart from constellation simulation can be transformed into standard **MILP / MIQP** formulations for revisit-time optimization.

---

## 1. Given: Contact Matrix (No SAT Operational Constraint)

- From satellite constellation simulation, we obtain an **unconstrained SAT–GS contact chart**.
- Each contact row is indexed as a constant entry.

Let the contact matrix be:

\[
A=
\begin{array}{c|c|c}
 \text{SAT} & \text{GS} & t_m \\
\hline
 a_{11} & a_{12} & a_{13} \\
 a_{21} & a_{22} & a_{23} \\
 a_{31} & a_{32} & a_{33} \\
\vdots & \vdots & \vdots  \\
 a_{N1} & a_{N2} & a_{N3} \\
\end{array}
\]

**Definitions**

- Total number of contacts: \(N\)
- Number of satellites: \(p\), where \(a_{i1} \in \{1,2,\ldots,p\}\)
- Number of ground stations: \(q\), where \(a_{i2} \in \{1,2,\ldots,q\}\)
- Time is sorted in ascending order:
\[
a_{13} \le a_{23} \le \cdots \le a_{N3}.
\]

---

## 2. Key Input Parameters for MILP/MIQP

- Unconstrained contact matrix: \(A\)
- Time vector: \(\bm{t} = A(:,3)\)
- Initial and final time: \(t_{\text{start}}, t_{\text{end}}\)
- Satellite cadence constraint: \(\tau\)
- Number of satellites: \(p\), number of GS: \(q\)
- Total contacts: \(N\)
- Contacts per satellite: \(|S_i|\)
- Contacts per GS: \(|G_j|\)
- Binary decision variable:
\[
\bm{x} \in \{0,1\}^N
\]

---

## 3. Selection Matrix Generation

### 3.1. \(E_{S_i}^{1}\): Map global contacts to each satellite’s contact sequence

\[
E_{S_i}^{1} : |S_i| \times N,\quad i=1,\ldots,p
\]

Construction outline:

1. Initialize \(E_{S_i}^{1}=\text{zeros}(|S_i|,N)\)
2. Collect indices where \(A(:,1)=i\Rightarrow a_{i1},a_{i2},\ldots,a_{i|S_i|}\)
3. Set one-hot rows:
\[
E_{S_i}^{1}(r, a_{ir}) = 1,\quad r=1,\ldots,|S_i|
\]

Then:

\[
\bm{x}_{S_i}=E_{S_i}^{1}\bm{x},\qquad \bm{t}_{S_i}=E_{S_i}^{1}\bm{t}
\]

---

### 3.2. \(E_{S_i,x}^{2}, E_{S_i,t}^{2}\): Pairwise-difference selection matrices for each satellite

Dimensions:

\[
E_{S_i,x}^{2},\, E_{S_i,t}^{2} \in \mathbb{R}^{\frac{|S_i|(|S_i|-1)}{2}\times |S_i|},\quad i=1,\ldots,p
\]

Construction (as given):

- For \(\alpha=1,\ldots,|S_i|-1\), create blocks \(E_{S_i,x,\alpha}^{2}\) and stack them.
- Similarly construct \(E_{S_i,t,\alpha}^{2}\) with \(-1\) at the \(\alpha\)-th column and \(+1\) at later contacts.

---

### 3.3. \(E_{G_j}^{1}\): Map global contacts to each GS revisit sequence (with boundaries)

\[
E_{G_j}^{1} : (|G_j|+2)\times (N+2),\quad j=1,\ldots,q
\]

Construction outline:

1. Initialize \(E_{G_j}^{1}=\text{zeros}(|G_j|+2,N+2)\)
2. Collect indices where \(A(:,2)=j\Rightarrow b_{j1},\ldots,b_{j|G_j|}\)
3. Set boundary rows:
\[
E_{G_j}^{1}(1,1)=1,\qquad E_{G_j}^{1}(|G_j|+2,N+2)=1
\]
4. Map contacts (with index shift):
\[
E_{G_j}^{1}(\alpha+1,b_{j\alpha}+1)=1,\quad \alpha=1,\ldots,|G_j|
\]

Affine forms:

\[
\bm{x}_{G_j}
=
E_{G_j}^{1}
\begin{bmatrix}
1\\ \bm{x}\\ 1
\end{bmatrix}
\quad\Rightarrow\quad \bm{x}_{G_j}\ \text{is affine w.r.t.}\ \bm{x}
\]

\[
\bm{t}_{G_j}
=
E_{G_j}^{1}
\begin{bmatrix}
t_{\text{start}}\\ \bm{t}\\ t_{\text{end}}
\end{bmatrix}
\quad\Rightarrow\quad \bm{t}_{G_j}\ \text{is affine w.r.t.}\ \bm{t}
\]

---

### 3.4. \(E_{G_j,x}^{2}, E_{G_j,t}^{2}\): Pairwise-difference selection matrices for each GS

Dimensions:

\[
E_{G_j,x}^{2},\, E_{G_j,t}^{2}\in \mathbb{R}^{\frac{(|G_j|+2)(|G_j|+1)}{2}\times (|G_j|+2)},\quad j=1,\ldots,q
\]

Construction follows the provided nested-loop stacking approach.

---

## 4. Derivation of \(L_1\), \(L_2\), \(L_\infty\) Revisit-Time Problems to MILP/MIQP

### 4.1. Cadence constraint \(\tau\) as linear inequality (Big-\(M\))

For satellite \(S_i\), the cadence constraint:

\[
\bm{t}_{S_i,\ell}-\bm{t}_{S_i,k}
\ge
\tau - \tau\left(2-\bm{x}_{S_i,\ell}-\bm{x}_{S_i,k}\right),
\quad \forall\,\ell>k,\ i=1,\ldots,p
\]

Using the difference–selection matrices:

\[
E_{S_i,x}^{2}\bm{x}_{S_i}\le
\mathbf{1}+\frac{1}{\tau}E_{S_i,t}^{2}\bm{t}_{S_i},
\qquad i=1,\ldots,p
\]

With \(\bm{x}_{S_i}=E_{S_i}^{1}\bm{x}\) and \(\bm{t}_{S_i}=E_{S_i}^{1}\bm{t}\), stacking all satellites yields:

\[
A\bm{x}\le \bm{b},
\quad
A=
\begin{bmatrix}
E^{2}_{S_{1,x}}E^{1}_{S_1}\\
E^{2}_{S_{2,x}}E^{1}_{S_2}\\
\vdots\\
E^{2}_{S_{p,x}}E^{1}_{S_p}
\end{bmatrix},
\quad
\bm{b}=\mathbf{1}+\frac{1}{\tau}
\begin{bmatrix}
E^{2}_{S_{1,t}}E^{1}_{S_1}\\
E^{2}_{S_{2,t}}E^{1}_{S_2}\\
\vdots\\
E^{2}_{S_{p,t}}E^{1}_{S_p}
\end{bmatrix}\bm{t}
\]

---

### 4.2. \(L_1\) optimization: maximize number of activated contacts

\[
\max\ \mathbf{1}^{\mathsf T}\bm{x}
\quad\text{s.t.}\quad
A\bm{x}\le \bm{b},\ \bm{x}\in\{0,1\}^N
\]

---

### 4.3–4.4. \(L_{\infty}\) optimization: minimize maximum revisit time

Introduce scalar \(R\) and enforce:

\[
\min\ R
\quad\text{s.t.}\quad
A\bm{x}\le \bm{b},\ \bm{x}\in\{0,1\}^N,\quad
R\mathbf{1}\ge C\bm{x}+\bm{d}
\]

where each GS contributes \((C_j,\bm{d}_j)\), and stacking \(j=1,\ldots,q\) yields \(C,\bm{d}\).

---

### 4.5–4.6. \(L_2\) optimization: MIQP minimizing squared revisit-time sum

Define \(\bm{y}\in\mathbb{R}^M\) where:

\[
M=\sum_{j=1}^{q}\frac{(|G_j|+2)(|G_j|+1)}{2},\qquad
\mathbf{1}\le \bm{y}\le 2\mathbf{1}
\]

Constraint:

\[
\bm{y}\ge E\bm{x}+\bm{f}
\]

Quadratic objective (square-sum of revisit time):

\[
\min\ (\bm{y}-\mathbf{1})^{\mathsf T}G(\bm{y}-\mathbf{1})
\]

Final MIQP:

\[
\min\ (\bm{y}-\mathbf{1})^{\mathsf T} G (\bm{y}-\mathbf{1})
\]
\[
\text{s.t.}\quad
A\bm{x}\le \bm{b},\ \bm{x}\in\{0,1\}^N
\]
\[
\bm{y}\ge E\bm{x}+\bm{f},\quad \mathbf{1}\le \bm{y}\le 2\mathbf{1},\quad \bm{y}\in\mathbb{R}^M
\]

with

\[
E=
\begin{bmatrix}
E^{2}_{G_1,x} E^{1}_{G_1} \\
E^{2}_{G_2,x} E^{1}_{G_2} \\
\vdots \\
E^{2}_{G_q,x} E^{1}_{G_q}
\end{bmatrix}
\begin{bmatrix}
\mathbf{0}_{1\times N}\\
\mathbf{I}_{N\times N}\\
\mathbf{0}_{1\times N}
\end{bmatrix},
\qquad
\bm{f}=
\begin{bmatrix}
E^{2}_{G_1,x} E^{1}_{G_1} \\
E^{2}_{G_2,x} E^{1}_{G_2} \\
\vdots \\
E^{2}_{G_q,x} E^{1}_{G_q}
\end{bmatrix}
\begin{bmatrix}
1\\
\mathbf{0}_{N\times 1}\\
1
\end{bmatrix}
\]

and

\[
G=
\begin{bmatrix}
t_{\text{start}}\\ \bm{t}\\ t_{\text{end}}
\end{bmatrix}^{\mathsf T}
\left(
\begin{bmatrix}
E^{2}_{G_1,t} E^{1}_{G_1}\\
E^{2}_{G_2,t} E^{1}_{G_2}\\
\vdots\\
E^{2}_{G_q,t} E^{1}_{G_q}
\end{bmatrix}^{\mathsf T}
\begin{bmatrix}
E^{2}_{G_1,t} E^{1}_{G_1}\\
E^{2}_{G_2,t} E^{1}_{G_2}\\
\vdots\\
E^{2}_{G_q,t} E^{1}_{G_q}
\end{bmatrix}
\right)
\begin{bmatrix}
t_{\text{start}}\\ \bm{t}\\ t_{\text{end}}
\end{bmatrix}
\]

---

## Notes

- This README preserves the derivation structure from the original LaTeX draft.
- If you want, I can also:
  - add **“How to reproduce (MATLAB scripts)”** section,
  - add **notation table** (symbols & dimensions),
  - or convert selection-matrix construction into **pseudo-code + MATLAB snippets**.

