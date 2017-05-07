## Schrodinger Equation Solver

This program is a universal Schrodinger equation solver implemented in Python and accelerated using numpy. Two versions of solving engine are provided for a comparison in efficiency.

```
	Author: Zhang Chang-kai
	E-mail: phy.zhangck@gmail.com
```

This program is licensed under General Public License-3.0. The author would like to acknowledge helpful discussions with Robert Brown.

### Usage

The program consists of mainly five files

- intfce.py  —  interface file, the entrance of the program
- config.py  —  stores all the configuration of the program
- engine_cplx.py  —  solving engine using numpy complex data type
- engine_real.py  —  solving engine using real data type
- visual.py  —  Python visualization of the solution

The program can be initiated by

```
  $ python intfce.py
```

Default settings will generate an animation for the solution. Settings can be modified so as to generate plot or any other types of visualization.