# The Conjugate Gradient Method for Solving Sparse Symmetric Positive Definite Linear Systems
**Report (PDF):** [[Zenodo DOI](https://doi.org/10.5281/zenodo.19674996)]  
**Code Repository:** [[GitHub link](https://github.com/jacekschmidt/Conjugate_Gradient_Method/edit/main/README.md)]
## Overview

This project implements the Conjugate Gradient Method (CGM) for solving large linear systems of the form:

A x = b

Such systems arise widely in applied mathematics, particularly in the numerical solution of partial differential equations (PDEs), optimisation problems, and data-driven applications.

While direct methods such as LU decomposition are robust, they become computationally infeasible for large-scale problems due to their cubic complexity. Iterative methods provide a scalable alternative, and the conjugate gradient method is especially effective for symmetric positive definite (SPD) systems.

This project focuses on the formulation, implementation, and performance of CGM, including comparisons with steepest descent and improvements through algorithmic optimisation.

---

## Mathematical Model

We consider linear systems:

A x = b

where:

- A ∈ R^(n × n) is symmetric positive definite (SPD)  
- x is the unknown solution vector  
- b is the known right-hand side  

The CGM solves this system by iteratively minimising the quadratic functional:

(1/2) xᵀ A x − bᵀ x

---

## Algorithms Implemented

### Conjugate Gradient Method (CGM)

- Constructs conjugate search directions  
- Guarantees convergence in at most n iterations  
- Typically converges much faster in practice  

### Fast Conjugate Gradient

- Optimised version avoiding repeated matrix-vector products  
- Improves computational efficiency  

### Steepest Descent (Baseline)

- Uses gradient direction only  
- Much slower convergence compared to CGM  
- Included for comparison  

---

## Model Problems

### 1D Poisson Equation

- Discretised using finite differences  
- Produces tridiagonal SPD matrix  

### 2D Poisson Equation

- Discretised on a square grid  
- Produces large sparse SPD matrix  
- Represents physical problems such as heat distribution or steady-state diffusion  

### Modified Systems

- Parameterised matrices used to study condition number effects  
- Demonstrates impact on convergence rate  

---

## Project Structure

```
.
├── src/
│   └── Conjugate_Gradient_Method.cpp
│
├── include/
│   ├── mvector.h
│   ├── mmatrix.h
│   ├── timing.h
│
├── results/
│   ├── solutionspace.eps
│   ├── contoursolutionspace.eps
│
├── Conjugate_Gradient_Method.pdf
└── README.md
```

---

## How to Run

### Requirements

- C++ compiler (e.g. g++)

### Compile

```bash
g++ -O2 main.cpp -o cg_solver
```

### Run

```bash
./cg_solver
```

The `main()` function includes multiple experiments, including:

- Small system validation  
- Comparison with steepest descent  
- Convergence behaviour analysis  
- Performance scaling with system size  
- 2D Poisson equation solver  

To run specific experiments, uncomment the relevant sections in `main()`.

---

## Experiments

### Convergence Behaviour

The CGM is tested on systems of increasing size. Results demonstrate:

- Rapid convergence compared to steepest descent  
- Iteration count significantly lower than system dimension  

### Effect of Condition Number

Modified systems illustrate:

- Higher condition numbers slow convergence  
- Algorithm performance strongly depends on matrix properties  

### Performance Scaling

Timing experiments show:

- Efficient scaling for large sparse systems  
- Fast CG implementation reduces computational cost  

### 2D Poisson Problem

The method is applied to a discretised 2D Poisson equation with fixed boundary conditions (e.g. steady-state heat distribution on a square plate).

---

## Results

### Numerical Solution of 2D Poisson Equation

- 3D surface plot shows smooth spatial variation of the solution  
- Contour plot highlights level sets and gradient structure  

These results demonstrate the effectiveness of CGM for solving large sparse PDE systems.

---

## Key Observations

- CGM significantly outperforms steepest descent  
- Convergence depends heavily on the condition number of A  
- Sparse structure is essential for computational efficiency  
- Fast CG implementation reduces redundant computations  
- PDE-derived systems are particularly well-suited to CGM  

---

## Notes

- Implementation uses custom vector and matrix classes  
- Designed for clarity and extensibility  
- Easily adaptable to include preconditioning techniques  

---

## Author

Jacek Schmidt 
jacek14schmidt@gmail.com
