-- ============================================
-- Main Entrance of Schrodinger Equation Solver
-- ============================================
-- 
--   Copyright (C) 2017 Zhang Chang-kai
--   Contact via: phy.zhangck@gmail.com
--   General Public License version 3.0
--
--   This is the main entrance of Schrodinger
--   equation solver. This file will load the 
--   configuartion file and call the solve
--   engine. Finally write the solution into
--   a txt file (default data.txt).

-- import solve engine and configuration
import SchdgerEq
import Configure

-- main: write data into a txt file
main = writeFile destination (show result)
    where result = every gap (solve initial step_t rec_tot)

-- extract every nth element of a list
every n xs = case drop (n-1) xs of
              (y:ys) -> y : every n ys
              [] -> []

-- End of Main Entrance of Shrodinger Sovler