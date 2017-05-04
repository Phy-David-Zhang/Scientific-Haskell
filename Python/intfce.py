#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Interface for Numerical Schrodinger Equation Solver

    # Copyright (C) 2017 Zhang Chang-kai #
    # Contact via: phy.zhangck@gmail.com #
    # General Public License version 3.0 #

'''Numerical Schrodinger Equation Solve Interface'''

# import engine and Visualization
import engine_ori as engine
import visual

# designate action
control = "animate"
i = 9000

_genSolution = engine.Schdger()
import time
import cProfile
start = time.time()
cProfile.run("_genSolution.solveEq()")
end = time.time()
print(end - start)

if control == "animate":
    visual.set_sol(_genSolution)
    visual.genAnimate()
    
elif control == "plot":
    visual.set_sol(_genSolution)
    visual.plot(i)
    

# End of Interface for Schrodinger Equation Solver
