# Redstone_MILP
Redstone's LEO satellite constellation optimization using Mixed-Integer Linear Programming

# Satellite Constellation Scheduling to MILP/MIQP Problem

This note summarizes how an unconstrained SAT–GS contact chart from constellation simulation can be transformed into standard **MILP / MIQP** formulations for revisit-time optimization.

---

## 1. Given: Contact Matrix (No SAT Operational Constraint)

- From satellite constellation simulation, we obtain an **unconstrained SAT–GS contact chart**.
- Each contact row is indexed as a constant entry.

Let the contact matrix be:

$$
A=
\begin{array}{c|c|c}
\text{SAT} & \text{GS} & t_m \\
\hline
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33} \\
\vdots & \vdots & \vdots \\
a_{N1} & a_{N2} & a_{N3}
\end{array}
$$

**Definitions**

- Total number of contacts: $N$
- Number of satellites: $p$, where $a_{i1} \in \{1,2,\ldots,p\}$
- Number of ground stations: $q$, where $a_{i2} \in \{1,2,\ldots,q\}$
- Time is sorted in ascending order:

$$
a_{13} \le a_{23} \le \cdots \le a_{N3}
$$

---

## 2. Key Input Parameters for MILP/MIQP

- Unconstrained contact matrix: $A$
- Time vector: $\mathbf{t} = A(:,3)$
- Initial and final time: $t_{\text{start}},\ t_{\text{end}}$
- Satellite cadence constraint: $\tau$
- Number of satellites: $p$, number of GS: $q$
- Total contacts: $N$
- Contacts per satellite: $|S_i|$
- Contacts per GS: $|G_j|$
- Binary decision variable:

$$
\mathbf{x} \in \{0,1\}^N
$$

---

## 3. Selection Matrix Generation

### 3.1. $E_{S_i}^{1}$: Map global contacts to each satellite’s contact sequence

$$
E_{S_i}^{1} : |S_i| \times N,\quad i=1,\ldots,p
$$

Construction outline:

1. Initialize $E_{S_i}^{1} = \text{zeros}(|S_i|,N)$  
2. Collect indices where $A(:,1)=i \Rightarrow a_{i1},\ldots,a_{i|S_i|}$  
3. Set one-hot rows:

$$
E_{S_i}^{1}(r,a_{ir}) = 1,\quad r=1,\ldots,|S_i|
$$

Then:

$$
\mathbf{x}_{S_i} = E_{S_i}^{1}\mathbf{x},\qquad
\mathbf{t}_{S_i} = E_{S_i}^{1}\mathbf{t}
$$

---

### 3.2. $E_{S_i,x}^{2},\ E_{S_i,t}^{2}$: Pairwise-difference selection matrices for each satellite

$$
E_{S_i,x}^{2},\ E_{S_i,t}^{2}
\in \mathbb{R}^{\frac{|S_i|(|S_i|-1)}{2}\times |S_i|},\quad i=1,\ldots,p
$$

- Constructed by stacking difference blocks over $\alpha=1,\ldots,|S_i|-1$  
- $E_{S_i,t}^{2}$ uses $-1/+1$ to encode time differences

---

### 3.3. $E_{G_j}^{1}$: Map global contacts to each GS revisit sequence

$$
E_{G_j}^{1} : (|G_j|+2)\times (N+2),\quad j=1,\ldots,q
$$


---

### 3.4. $E_{G_j,x}^{2},\ E_{G_j,t}^{2}$: Pairwise-difference selection matrices for each GS

$$
E_{G_j,x}^{2},\ E_{G_j,t}^{2}
\in \mathbb{R}^{\frac{(|G_j|+2)(|G_j|+1)}{2}\times(|G_j|+2)},\quad j=1,\ldots,q
$$

---

## 4. Revisit-Time Optimization Formulations

### 4.1. Cadence constraint (Big-$M$)

$$
\mathbf{t}_{S_i,\ell}-\mathbf{t}_{S_i,k}
\ge
\tau-\tau(2-\mathbf{x}_{S_i,\ell}-\mathbf{x}_{S_i,k})
$$

Equivalent linear inequality:

$$
E_{S_i,x}^{2}\mathbf{x}_{S_i}
\le
\mathbf{1}+\frac{1}{\tau}E_{S_i,t}^{2}\mathbf{t}_{S_i}
$$

Stacked form:

$$
A\mathbf{x}\le\mathbf{b}
$$

---

### 4.2. $L_1$ Optimization (MILP)

$$
\max\ \mathbf{1}^\mathsf{T}\mathbf{x}
\quad\text{s.t.}\quad
A\mathbf{x}\le\mathbf{b},\ \mathbf{x}\in\{0,1\}^N
$$

---

### 4.3–4.4. $L_\infty$ Optimization (MILP)

$$
\min\ R
\quad\text{s.t.}\quad
A\mathbf{x}\le\mathbf{b},\quad
R\mathbf{1}\ge C\mathbf{x}+\mathbf{d}
$$

---

### 4.5–4.6. $L_2$ Optimization (MIQP)

$$
\min\ (\mathbf{y}-\mathbf{1})^\mathsf{T}G(\mathbf{y}-\mathbf{1})
\quad\text{s.t.}\quad
\mathbf{y}\ge E\mathbf{x}+\mathbf{f},\quad
\mathbf{1}\le\mathbf{y}\le2\mathbf{1}
$$

---

## Notes

- Fully compatible with **GitHub MathJax**
- Preserves the full mathematical derivation
- Ready for **research-grade public repository**

