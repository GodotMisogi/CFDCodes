import AeroModules

main = do
    putStrLn "------------------------------"
    putStrLn "|       H a s k F O I L      |"
    putStrLn "------------------------------"
    putStrLn "Enter coordinates file path:"
    -- putStrLn "Enter NACA 4-digit series"
    filename <- getLine
    -- putStrLn "Enter number of panels"
    -- num <- getLine
    -- putStrLn "Enter chord length (in m)"
    -- chord <- getLine
    -- putStrLn "Closed trailing edge?"
    -- te <- getLine
    coords <- importCoords filename
    -- coords <- let [a, b, c, d] = read $ take 4 filename :: Double in importCoords $ naca4 ((a, b, c, d)) num (if head getLine == "y" then "closed" else "no u")
    putStrLn "Enter speed:"
    speed <- getLine
    putStrLn "Enter angle (in degrees):"
    angle <- getLine
    putStrLn "Pick panel type: dblt or dblt-src"
    panel_type <- getLine
    -- putStrLn $ show coords
    let airfoil = if panel_type == "dblt-src" then doubletSourcePanels coords else doubletPanels coords
    -- putStrLn $ show airfoil
    let uniform2D = Uniform2D (read speed) (read angle)
    let (cl, kcl, cps) = aeroCoefficients airfoil uniform2D
    putStrLn "Lift Coefficient:"
    putStrLn $ show (cl, kcl)
    -- putStrLn $ show (influenceMatrix airfoil uniform2D)
    -- putStrLn $ show (boundaryCondition airfoil uniform2D)
    -- putStrLn $ show cps

