## Schrodinger Equation Solver

This program is a 1D universal Schrodinger equation solver implemented in pure C language. This program is to find out the limit of the compuatational efficiency for numerical Schrodinger equation solver.

```
	Author: Zhang Chang-kai
	E-mail: phy.zhangck@gmail.com
```

This program is licensed under General Public License-3.0.

###  Usage

Makefile is provided. Use `make` to automatically compile, execute and write data into `data.log`. A visualizer written in Python is used to generate animation for the solution.

### Configuration

Basic configurations lie in the front of the program and the main function. Settings as spacial and temporal step size, number of spacial and temporal steps are located in the beginning of the program using macro definition. The initial condition is set in the main function.