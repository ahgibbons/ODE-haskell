module Main where

import System.IO
import Data.List
import Control.Applicative

import Solvers (euler_ivp, rk4_ivp)

h = 0.01
f = shm
t0 = 0
tf = 100
y0 = [10,0]

shm :: a -> [Double] -> [Double]
shm _ [x,v] = [-v,x]

result_to_csv :: String -> ([Double], [[Double]]) -> String
result_to_csv sep = unlines . map (intercalate sep) 
              . (map . map) show . liftA2 (zipWith (:)) fst snd 

main :: IO ()
main = do
    writeFile "euler.dat" $ result_to_csv " " $ euler_ivp h f t0 tf y0
    writeFile "rk4.dat" $ result_to_csv " " $ rk4_ivp h f t0 tf y0
    putStrLn "Finished!"
