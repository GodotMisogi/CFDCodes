getsFloats :: FilePath -> IO [Float]
getsFloats path = do
    contents <- readFile path
    let someFloats = map read  . lines $ contents
    return someFloats

fileData :: String -> IO [[Float]]
fileData str = readFile str >>=
                \file -> return $ map ((map $ \x -> read x::Float) . words $) $ lines file