import Data.Complex as Complex
import Bits
import Data.WAVE as WAVE
import qualified Data.List as List (transpose, zip5, groupBy, sortBy)
import qualified System (getArgs)
import qualified Numeric.FFT (fft, ifft)



{------------------------------------------------------------------- 
  some settings for tweaking

    sample frequency: should match that of source audio file

    frame size: >= 2048 samples per frame sounds decent

    over sampling: >= 4 sounds decent
  
    frequency multiplier: 2**(x/12) increases pitch by x half-steps
 -------------------------------------------------------------------}
sampleFreq :: Double
sampleFreq = 44100

frameSize :: Int
frameSize = 1024*2 -- do make this a power of 2

overSamp :: Int
overSamp = 4 -- do make this a power of 2

freqMultiplier :: Double
freqMultiplier = 2**(( 7 )/12)



{------------------------------
  just converting some numbers
 ------------------------------}
frameSizeD :: Double
frameSizeD = fromIntegral frameSize

overSampD :: Double
overSampD = fromIntegral overSamp

frameLen = frameSizeD / sampleFreq
freqSpacing = sampleFreq / frameSizeD



{--------------------------------------
  several options for the window
  hann seemed to work best (it should)
 --------------------------------------}
window :: [Double]

--window = [1,1..]
--invWindow = window

window = map ((0.5+) . ((-0.5)*) . cos . ((2 * pi / (frameSizeD-1))*)) (take frameSize [0,1..]) -- hann
invWindow = window

--window = map ((0.5434782608695652+) . ((-0.45652173913043476)*) . cos . ((2 * pi / (frameSizeD-1))*)) (take frameSize [0,1..]) -- hamming 25/46 and 21/46
--invWindow = window

--blackman a0 a1 a2 x = a0 - a1 * cos (2 * pi * x / (frameSizeD - 1)) + a2 * cos (4 * pi * x / (frameSizeD - 1))
--window = map (blackman 0.42 0.5 0.08) (take frameSize [0,1..])
--invWindow = [1,1..]



{----------------------------------------
  converts wave to samples (and inverse)
 ----------------------------------------}
getSamples :: WAVE -> WAVESamples
getSamples (WAVE h ss) = ss

restoreSamples :: WAVESamples -> WAVE -> WAVE
restoreSamples newWAVESamples (WAVE oldWAVEHeader oldWAVESamples) = WAVE oldWAVEHeader newWAVESamples



{-------------------------------------------
  converts samples to doubles (and inverse)
 -------------------------------------------}
samplesToDoubles :: WAVESamples -> [Double]
samplesToDoubles = map (sampleToDouble . head)

doublesToSamples :: [Double] -> WAVESamples
doublesToSamples x = List.transpose [map doubleToSample x]



{---------------------------------------------------
  converts doubles to complex doubles (and inverse)
  imaginary component is zero
  just setting up for FFT input
 ---------------------------------------------------}
mkComplexList :: [Double] -> [Complex Double]
mkComplexList x = zipWith id (map (Complex.:+) x) [0, 0 ..]

mkDoubleList = map Complex.realPart



{--------------------------------------------------
  takes a list of doubles (which were samples) and
  splits it up into frames (and inverse)
 --------------------------------------------------}
framer :: Int -> Int -> [Double] -> [[Double]]
framer _ _ [] = []
framer size overlap x = take size x : framer size overlap (drop (div size overlap) x)

deframer' :: Int -> Int -> [[Double]] -> [Double]
deframer' _ _ [] = []
deframer' size overlap (x:xs) = zipWith (+) (x ++ [0,0..]) (take padding [0,0..] ++ deframer' size overlap xs)
	where
	padding = div size overlap

deframer :: Int -> Int -> [[Double]] -> [Double]
deframer size overlap x = map (/(fromIntegral overlap :: Double)) (deframer' size overlap x)



{-------------------------------------------------------------------------------------------------------
  this is where the FFT (and inverse) takes place
  takes a list of sample points and returns a pair of lists - amplitudes and phases
  notice that the inverse reverses the list's order
  if I remember correctly, the "R" versions double the number of bins by using only the real components
 -------------------------------------------------------------------------------------------------------}
transform :: [Double] -> ([Double], [Double])
transform x = (map Complex.magnitude rawFFT, map Complex.phase rawFFT)
	where
	rawFFT = (Numeric.FFT.fft . mkComplexList) x

transformR :: [Double] -> ([Double], [Double])
transformR x = (map (2*) {-to copy smb-} $ chop $ fst $ transform x, chop $ snd $ transform x)
	where
	chop = take (div (length x) 2 + 1)

invTransform :: ([Double], [Double]) -> [Double]
invTransform (x, y) = (mkDoubleList . reverse . Numeric.FFT.ifft) (zipWith Complex.mkPolar x y)

invTransformR :: ([Double], [Double]) -> [Double]
invTransformR (x, y) = invTransform (x ++ (drop 1 $ reverse $ drop 1 x), y ++ (map ((-1)*) $ drop 1 $ reverse $ drop 1 y))



{------------------------------------------------------------------------------------------------
  estimates the "true" bin frequencies by comparing change in phase between (oversampled) frames
  (and inverse)
  takes a list of [(mag, phase)] and returns a list of [(mag, trueFreq)]
 ------------------------------------------------------------------------------------------------}
trueBin :: Double -> Double -> Double -> Double
trueBin overSamp bin phaseDiff = overSamp * ((phaseDiff/2/pi) - fromIntegral qpd2/2)
	where
	qpd2 :: Int
	qpd2 = if qpd >= 0 then qpd + (qpd .&. 1) else qpd - (qpd .&. 1)
	qpd :: Int
	qpd = truncate((phaseDiff/pi) - (2*bin/overSamp)) -- changed from floor to ceiling because smb seems to do that

trueBinList' :: Double -> [[Double]] -> [[Double]]
trueBinList' _ [x] = []
trueBinList' overSamp (x0:x1:xs) = zipWith (trueBin overSamp) [0..] (zipWith (-) x1 x0) : trueBinList' overSamp (x1:xs)

trueBinList :: Double -> [[Double]] -> [[Double]]
trueBinList overSamp x = trueBinList' overSamp (take (length $ head x) [0,0..] : x)


invTrueBin :: Double -> Double -> Double -> Double
invTrueBin overSamp bin trueBin = trueBin {- - bin -} * 2 * pi / overSamp

invTrueBinList :: Double -> [[Double]] -> [[Double]]
invTrueBinList overSamp x = scanl1 (zipWith (+)) (map (zipWith ($) ((map $ invTrueBin overSamp) [0..])) x)



{-----------------------------
  some black magic I guess :\
 -----------------------------}
gComper :: (Int, Double, Double) -> (Int, Double, Double) -> Bool
gComper (roundedScaledBin0, _, _) (roundedScaledBin1, _, _) = roundedScaledBin0 == roundedScaledBin1

grouper :: [(Int, Double, Double)] -> [[(Int, Double, Double)]]
grouper = List.groupBy gComper


combinerOp :: (Int, Double, Double) -> (Int, Double, Double) -> (Int, Double, Double)
combinerOp (rSBin0, x0, _) (_, x1, y1) = (rSBin0, x0 + x1, y1) -- tweek this!!!!!

combiner :: [[(Int, Double, Double)]] -> [(Int, Double, Double)]
combiner = map (foldl1 combinerOp)


expander' :: Int -> [(Int, Double, Double)] -> [(Int, Double, Double)]
expander' _ [] = repeat (0, 0, 0)
expander' k ((rSBin, amp, trueBin):xs) = if k == rSBin
                        then (rSBin, amp, trueBin) : expander' (k+1) xs
                        else (0, 0, 0) : expander' (k+1) ((rSBin, amp, trueBin):xs)

expander :: [(Int, Double, Double)] -> [(Int, Double, Double)]
expander = expander' 0

addScaledIndexes' :: Double -> Double -> [(Double, Double)] -> [(Int, Double, Double)]
addScaledIndexes' _ _ [] = []
addScaledIndexes' scale k ((amp, trueBin):xs) = (round (scale*k), amp, trueBin) :
												addScaledIndexes' scale (k+1) xs

addScaledIndexes :: Double -> [(Double, Double)] -> [(Int, Double, Double)]
addScaledIndexes scale = addScaledIndexes' scale 0

rmScaledIndexes :: [(Int, Double, Double)] -> [(Double, Double)]
rmScaledIndexes [] = []
rmScaledIndexes ((rSBin, mag, trueBin):xs) = (mag, trueBin) : rmScaledIndexes xs



{----------------------------------------------------
  this is where we take the info in frequency domain
  and actually shift the pitches
 ----------------------------------------------------}
shiftPitch' :: Int -> Double -> ([Double], [Double]) -> ([Double], [Double])
shiftPitch' frameSize scale (mags, bins)
	= (unzip . take (div frameSize 2 + 1) . rmScaledIndexes . expander . combiner . grouper . addScaledIndexes scale) $
		zip mags (map (scale*) bins)

shiftPitch :: Int -> Double -> [([Double], [Double])] -> [([Double], [Double])]
shiftPitch frameSize scale = map (shiftPitch' frameSize scale)



{--------------------------------------------------------------------
  put it all together
  it looks weird because this was just the easiest way to experiment
 --------------------------------------------------------------------}
doublesToDoubles :: Int -> Double -> Double -> [Double] -> [Double]
doublesToDoubles frameSize overSamp pitchScale x =
	deframer frameSize (truncate overSamp) $ map (zipWith (*) invWindow . invTransformR) (zip c (invTrueBinList overSamp d))
	where
	(a, b) = unzip $ map (transformR . zipWith (*) window) (framer frameSize (truncate overSamp) x)
	(c, d) = unzip $ shiftPitch frameSize pitchScale (zip a (trueBinList overSamp b))

bigThing2 :: WAVE -> WAVE
bigThing2 audio = restoreSamples finishedSamples audio
	where
	doubs =  (samplesToDoubles . getSamples) audio
	finishedSamples = (doublesToSamples . doublesToDoubles frameSize overSampD freqMultiplier) doubs



{------------------------
  and the non-pure stuff
 ------------------------}
main = do
	args <- System.getArgs
	audio <- getWAVEFile (args !! 0)
	let b = bigThing2 audio
	putWAVEFile (args !! 1) b
