module Solvers where

import Data.List
import Control.Applicative

infixl 6 .+.
infixl 7 *.
a .+. b = zipWith (+) a b
a *. b  =  (a*) <$> b


eulerstep :: Double -> (Double -> [Double] -> [Double]) -> (Double, [Double]) 
          -> (Double, [Double])
eulerstep h f (t,y) = (t+h, y .+. h *. (f t y) )


euler_ivp :: Double -> (Double -> [Double] -> [Double])
          -> Double -> Double -> [Double]
          -> ([Double], [[Double]])  
euler_ivp h f t0 tf y0 = unzip $ takeWhile (\(t,_) -> t*signum h' < tf*signum h') 
                       $ iterate (eulerstep h' f) (t0, y0)
    where
        h' = h * signum (tf-t0)

rk4_step :: Double -> (Double -> [Double] -> [Double])
         -> (Double, [Double]) -> (Double, [Double])
rk4_step h f (t,y) = (t+h, y .+. (h/6) *. (k1.+.2*.k2.+.2*.k3.+.k4))
    where
        k1 = f t y
        k2 = f (t + h/2) (y .+. ((h/2)*.k1))
        k3 = f (t + h/2) (y .+. ((h/2)*.k2))
        k4 = f (t + h)   (y .+. (h*.k3))

rk4_ivp h f t0 tf y0 = unzip $ takeWhile (\(t,_) -> t*signum h' < tf*signum h') 
                       $ iterate (rk4_step h' f) (t0, y0)
    where
        h' = h * signum (tf-t0)

{-rk4 :: Double -> (Double -> [Double] -> [Double])
          -> Double -> Double -> [Double]
          -> ([Double], [[Double]])
rk4 h f t0 tf y0 =
    where-}
