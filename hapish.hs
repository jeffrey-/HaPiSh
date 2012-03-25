import Data.Complex as Complex
import Bits
import Data.WAVE as WAVE
import qualified List (transpose, zip5, groupBy, sortBy)
import qualified System (getArgs)
import qualified Numeric.FFT (fft, ifft)


sampleFreq :: Double
sampleFreq = 44100

frameSize :: Int
frameSize = 1024*2 -- do make this a power of 2

overSamp :: Int
overSamp = 4 -- do make this a power of 2

freqMultiplier, ampMult :: Double
freqMultiplier = 2**((-12)/12)
ampMult = 1 -- 4 -- 0.001
--
frameSizeD :: Double
frameSizeD = fromIntegral frameSize

overSampD :: Double
overSampD = fromIntegral overSamp

frameLen = frameSizeD / sampleFreq
freqSpacing = sampleFreq / frameSizeD

--                window = -.5*cos(2.*M_PI*(double)k/(double)fftFrameSize)+.5;
blackman = (\a0 a1 a2 x -> a0 - a1 * cos (2 * pi * x / (frameSizeD - 1)) + a2 * cos (4 * pi * x / (frameSizeD - 1)))

window :: [Double]
--window = [1,1..]
window = map ((0.5+) . ((-0.5)*) . cos . ((2 * pi / (frameSizeD-1))*)) (take frameSize [0,1..]) -- hann
invWindow = window
--window = map ((0.5434782608695652+) . ((-0.45652173913043476)*) . cos . ((2 * pi / (frameSizeD-1))*)) (take frameSize [0,1..]) -- hamming 25/46 and 21/46
--invWindow = window
--window = map (blackman 0.42 0.5 0.08) (take frameSize [0,1..])
--invWindow = [1,1..]


getSamples :: WAVE -> WAVESamples
getSamples (WAVE h ss) = ss

restoreSamples :: WAVESamples -> WAVE -> WAVE
restoreSamples newWAVESamples (WAVE oldWAVEHeader oldWAVESamples) = WAVE oldWAVEHeader newWAVESamples
--
samplesToDoubles :: WAVESamples -> [Double]
samplesToDoubles [] = []
samplesToDoubles (x:xs) = sampleToDouble (x !! 0) : (samplesToDoubles xs)

doublesToSamples :: [Double] -> WAVESamples
doublesToSamples x = List.transpose [(map doubleToSample x)]
--
mkComplexList :: [Double] -> [Complex Double]
mkComplexList x = zipWith id (map (Complex.:+) x) [0, 0 ..]

mkDoubleList :: [Complex Double] -> [Double]
mkDoubleList x = map Complex.realPart x
--
framer :: Int -> Int -> [Double] -> [[Double]]
framer _ _ [] = []
framer size overlap x = (take size x) : framer size overlap (drop (div size overlap) x)

deframer' :: Int -> Int -> [[Double]] -> [Double]
deframer' _ _ [] = []
deframer' size overlap (x:xs) = zipWith (+) (x ++ [0,0..]) ((take padding [0,0..]) ++ deframer' size overlap xs)
	where
	padding = div size overlap

--deframer' :: Int -> Int -> [[Double]] -> [Double]
--deframer' _ _ [] = []
--deframer' size overlap (x:xs) = offset x (deframer' size overlap xs)
--	where
--	offset a b = zipWith (+) (a ++ padding) (padding ++ b)
--	padding = take (div size overlap) [0,0..]

deframer :: Int -> Int -> [[Double]] -> [Double]
deframer size overlap x = map (/((fromIntegral overlap) :: Double)) (deframer' size overlap x)

--deframer :: Int -> Int -> [[Double]] -> [Double]
--deframer size overlap x = foldl (++) [] (map ((take takeN) . (drop 0)) x)
--	where
--	takeN = div size overlap
--	dropN = div (size - takeN) 2

--
reorder :: [a] -> [a] 
reorder (x:xs) = x : (reverse xs)
--

-- takes a list of sample points and returns a pair of lists - amplitudes and phases
transform :: [Double] -> ([Double], [Double])
transform x = (map Complex.magnitude rawFFT, map Complex.phase rawFFT)
	where
	rawFFT = (Numeric.FFT.fft . mkComplexList) x

transformR :: [Double] -> ([Double], [Double])
transformR x = ((map (2*) {-to copy smb-}) $ chop $ fst $ transform x, chop $ snd $ transform x)
	where
	chop = (take ((div (length x) 2) + 1))

invTransform :: ([Double], [Double]) -> [Double]
invTransform (x, y) = (mkDoubleList . reorder . Numeric.FFT.ifft) (zipWith (Complex.mkPolar) x y)

-- this will fuck up the amplitude!
invTransformR :: ([Double], [Double]) -> [Double]
invTransformR (x, y) = invTransform (x ++ (drop 1 $ reverse $ drop 1 x), y ++ (map ((-1)*) $ drop 1 $ reverse $ drop 1 y))

-- THIS IS WHERE SOME REAL NUMERICAL DIFFERENTIATION COULD HELP
-- takes a list of [(mag, phase)] and returns a list of [(mag, trueFreq)]
trueBinOld :: Double -> Double -> Double -> Double
trueBinOld overSamp bin phaseDiff = bin + (phaseDiff * overSamp / (2 * pi))

trueBin :: Double -> Double -> Double -> Double
trueBin overSamp bin phaseDiff = overSamp * ((phaseDiff/2/pi) - (fromIntegral qpd2)/2)
	where
	qpd2 :: Int
	qpd2 = if (qpd >= 0) then qpd + (qpd .&. 1) else qpd - (qpd .&. 1)
	qpd :: Int
	qpd = truncate((phaseDiff/pi) - (2*bin/overSamp)) -- changed from floor to ceiling because smb seems to do that

trueBinList' :: Double -> [[Double]] -> [[Double]]
trueBinList' _ [x] = []
trueBinList' overSamp (x0:x1:xs) = map (uncurry $ trueBin overSamp) (zip [0..] (zipWith (-) x1 x0))
									: trueBinList' overSamp (x1:xs)

trueBinList :: Double -> [[Double]] -> [[Double]]
trueBinList overSamp x = trueBinList' overSamp ((take (length $ x !! 0) [0,0..]) : x)

invTrueBin :: Double -> Double -> Double -> Double
invTrueBin overSamp bin trueBin = (trueBin {- - bin -}) * 2 * pi / overSamp

--invTrueBinList :: Double -> [[Double]] -> [[Double]]
--invTrueBinList overSamp x = List.transpose $ map phaseDiffs $ List.transpose x
--	where
--	phaseDiffs y = scanl1 (+) (map (uncurry $ invTrueBin overSamp) (zip [0..] y))

invTrueBinList :: Double -> [[Double]] -> [[Double]]
invTrueBinList overSamp x = scanl1 (zipWith (+)) (map (zipWith ($) ((map $ invTrueBin overSamp) [0..])) x)


dPhase5 :: Double -> (Double, Double, Double, Double, Double) -> Double
dPhase5 overSamp (p0, p1, p2, p3, p4) = overSamp * (p0 - p4 + (8 * (p3 - p1))) / 12

trueBinPlus :: Double -> Double -> Double -> Double
trueBinPlus overSamp bin dPhase = dPhase / (2 * pi) - bigMap (dPhase / (2 * pi) - bin)
	where
	map :: Int -> Int
	map x = if (x >= 0) then x + (x .&. 1) else x - (x .&. 1)
	bigMap :: Double -> Double
	bigMap x = overSamp / 2 * (fromIntegral $ map $ truncate $ 2 / overSamp * x)

trueBinList5' :: Double -> [[Double]] -> [[Double]]
trueBinList5' _ [x] = []
trueBinList5' overSamp (x0:x1:x2:x3:x4:xs)
	= map (uncurry $ trueBinPlus overSamp) (zip [0..] (map (dPhase5 overSamp) (List.zip5 x0 x1 x2 x3 x4)))
		: trueBinList5' overSamp (x1:x2:x3:x4:xs)

trueBinList5 :: Double -> [[Double]] -> [[Double]]
trueBinList5 overSamp x = trueBinList5' overSamp ([zeros, zeros] ++ x ++ [zeros, zeros])
	where
	zeros = take (length $ x !! 0) [0,0..]


dPhase3 :: Double -> (Double, Double, Double) -> Double
dPhase3 overSamp (p0, p1, p2) = overSamp * (p2 - p0) / 2

trueBinList3' :: Double -> [[Double]] -> [[Double]]
trueBinList3' _ [x] = []
trueBinList3' overSamp (x0:x1:x2:xs)
	= map (uncurry $ trueBinPlus overSamp) (zip [0..] (map (dPhase3 overSamp) (zip3 x0 x1 x2)))
		: trueBinList3' overSamp (x1:x2:xs)

trueBinList3 :: Double -> [[Double]] -> [[Double]]
trueBinList3 overSamp x = trueBinList3' overSamp ([zeros] ++ x ++ [zeros])
	where
	zeros = take (length $ x !! 0) [0,0..]


-- wtf was I thinking?
--trueBinList' :: Double -> Double -> [Double] -> [Double]
--trueBinList' _ _ [x] = []
--trueBinList' overSamp 0 (phase0:phases) = trueBin overSamp 0 phase0 : trueBinList' overSamp 1 (phase0:phases)
--trueBinList' overSamp bin (phase0:phase1:phases) = trueBin overSamp bin (phase1 - phase0)
--	: trueBinList' overSamp (bin+1) (phase1:phases)
--
--trueBinList :: Double -> [Double] -> [Double]
--trueBinList overSamp x = trueBinList' overSamp 0 x
--
--invTrueBinList :: Double -> [Double] -> [Double]
--invTrueBinList overSamp trueBins = scanl1 (+) (zipWith ($) (map (invTrueBin overSamp) [0..]) trueBins)



--shiftPitch' :: [(Double, Double)] -> [(Double, Double)]
--shiftPitch' [] = []
--shiftPitch' (old:olds) = case ((snd old) - ((fromIntegral (length (shiftPitch' olds)) :: Double) < (-0.5)) of
--	True	-> (((fst ((shiftPitch' olds news) !! 0)) + (fst old)), (snd old)) : drop 1 (shiftPitch' olds news)
--	False	-> shiftPitch' (old:olds) ((0,0):news)

-- try this with the second (maximum [0, snd old]) as just (snd old)
{-
shiftPitch'' :: ([(Double, Double)], [(Double, Double)]) -> ([(Double, Double)], [(Double, Double)])
shiftPitch'' ([], new) = ([], new)
shiftPitch'' ((old:olds), news) = case ((maximum [0, snd old]) - ((fromIntegral (length news)) :: Double) < (-0.5)) of
	True	-> (olds, ((fst old) + (fst (news !! 0)), (snd old{-maximum [0, snd old]-})) : (drop 1 news))
	False	-> (old:olds, (0, (fromIntegral (length news))):news)

shiftPitch' :: [(Double, Double)] -> [(Double, Double)]
shiftPitch' x = snd $ until ((==[]) . fst) shiftPitch'' (x, [])

shiftPitch :: Double -> ([Double], [Double]) -> ([Double], [Double])
shiftPitch scale (mag, freq) = unzip $ reverse $ shiftPitch' $ zip mag (map (scale*) freq)

makeRightSize :: Int -> ([Double], [Double]) -> ([Double], [Double])
makeRightSize size (x, y) 
	| (length x == size)	= (x, y)
	| (length x < size)		= (take size (x ++ [0,0..]), take size (y ++ (drop (length x) [0..])))
	| (length x > size)		= (take size x, take size y)

shiftPitchSize :: Double -> [([Double], [Double])] -> [([Double], [Double])]
shiftPitchSize scale x = map (makeRightSize (length $ fst (x !! 0))) $ map (shiftPitch scale) x
-}

insertify :: Double -> (Double, Double) -> Int -> [(Double, Double)] -> [(Double, Double)]
insertify scale (y, z) k stuff = (take k stuff) ++ ((fst (stuff !! k) + y, z*scale) : (drop (k+1) stuff))

slowShiftPitch'' :: Int -> Double -> Int -> ([(Double, Double)], [(Double, Double)]) -> ([(Double, Double)], [(Double, Double)])
slowShiftPitch'' _ _ _ ([], news) = ([], news)
slowShiftPitch'' frameSize scale k ((old:olds), news) = case (floor ((fromIntegral (k*2)) * scale) <= frameSize) of
	True	-> slowShiftPitch'' frameSize scale (k+1) (olds, insertify scale old (floor((fromIntegral k)*scale)) news)
	False	-> ((old:olds), news)
--	where
--	insertify :: (Double, Double) -> Int -> [(Double, Double)] -> [(Double, Double)]
--	insertify (y, z) k stuff = (take k stuff) ++ ((fst (stuff !! k) + y, z*scale) : (drop (k+1) stuff))

slowShiftPitch' :: Int -> Double -> ([Double], [Double]) -> ([Double], [Double])
slowShiftPitch' frameSize scale (x, y) = unzip $ snd $ slowShiftPitch'' frameSize scale 0 (zip x y, replicate (div frameSize  2 + 1) (0, 0))

slowShiftPitch :: Int -> Double -> [([Double], [Double])] -> [([Double], [Double])]
slowShiftPitch frameSize scale x = map (slowShiftPitch' frameSize scale) x


--sComper :: (Double, Double) -> (Double, Double) -> Ordering
--sComper (_, x) (_, y) = if (x > y) then GT else LT
--
--sorter :: [(Double, Double)] -> [(Double, Double)]
--sorter x = List.sortBy sComper x

--gComper :: (Double, Double) -> (Double, Double) -> Bool
--gComper (_, x) (_, y) = round x == round y
--
--grouper :: [(Double, Double)] -> [[(Double, Double)]]
--grouper x = List.groupBy gComper x

gComper :: (Int, Double, Double) -> (Int, Double, Double) -> Bool
gComper (roundedScaledBin0, _, _) (roundedScaledBin1, _, _) = roundedScaledBin0 == roundedScaledBin1

grouper :: [(Int, Double, Double)] -> [[(Int, Double, Double)]]
grouper x = List.groupBy gComper x

--combinerOp :: (Double, Double, Double) -> (Double, Double) -> (Double, Double)
--combinerOp (x0, x1) (y0, y1) = (x0 + y0, x1) -- tweek this!!!!!!
--
--combiner :: [[(Double, Double)]] -> [(Double, Double)]
--combiner x = map (foldl1 combinerOp) x

combinerOp :: (Int, Double, Double) -> (Int, Double, Double) -> (Int, Double, Double)
combinerOp (rSBin0, x0, _) (_, x1, y1) = (rSBin0, x0 + x1, y1) -- tweek this!!!!!

combiner :: [[(Int, Double, Double)]] -> [(Int, Double, Double)]
combiner x = map (foldl1 combinerOp) x

--expander' :: Int -> [(Double, Double)] -> [(Double, Double)]
--expander' _ [] = repeat (0, 0)
--expander' k (x:xs) = if (k == round (snd x))
--						then (x : expander' (k+1) xs)
--						else ((0, 0) : expander' (k+1) (x:xs))
--
--expander :: [(Double, Double)] -> [(Double, Double)]
--expander x = expander' 0 x

expander' :: Int -> [(Int, Double, Double)] -> [(Int, Double, Double)]
expander' _ [] = repeat (0, 0, 0)
expander' k ((rSBin, amp, trueBin):xs) = if (k == rSBin)
                        then ((rSBin, amp, trueBin) : expander' (k+1) xs)
                        else ((0, 0, 0) : expander' (k+1) ((rSBin, amp, trueBin):xs))

expander :: [(Int, Double, Double)] -> [(Int, Double, Double)]
expander x = expander' 0 x

addScaledIndexes' :: Double -> Double -> [(Double, Double)] -> [(Int, Double, Double)]
addScaledIndexes' _ _ [] = []
addScaledIndexes' scale k ((amp, trueBin):xs) = (round (scale*k), amp, trueBin) :
												addScaledIndexes' scale (k+1) xs

addScaledIndexes :: Double -> [(Double, Double)] -> [(Int, Double, Double)]
addScaledIndexes scale x = addScaledIndexes' scale 0 x

rmScaledIndexes :: [(Int, Double, Double)] -> [(Double, Double)]
rmScaledIndexes [] = []
rmScaledIndexes ((rSBin, mag, trueBin):xs) = (mag, trueBin) : rmScaledIndexes xs

--shiftPitch' :: Int -> Double -> ([Double], [Double]) -> ([Double], [Double])
--shiftPitch' frameSize scale (mags, freqs)
--	= (unzip . (take (div frameSize 2 + 1)) . expander . combiner . grouper . sorter) $ 
--		zip mags (map (\x -> maximum [0, x]) (map (scale*) freqs))

shiftPitch' :: Int -> Double -> ([Double], [Double]) -> ([Double], [Double])
shiftPitch' frameSize scale (mags, bins)
	= (unzip . (take (div frameSize 2 + 1)) . rmScaledIndexes . expander . combiner . grouper . (addScaledIndexes scale)) $
		(zip mags (map (scale*) bins))

shiftPitch :: Int -> Double -> [([Double], [Double])] -> [([Double], [Double])]
shiftPitch frameSize scale x = map (shiftPitch' frameSize scale) x




doublesToDoubles :: Int -> Double -> Double -> [Double] -> [Double]
doublesToDoubles frameSize overSamp pitchScale x =
	map (ampMult*) $ deframer frameSize (truncate overSamp) $ map (zipWith (*) invWindow) $ map invTransformR (zip c (invTrueBinList overSamp d))
	where
	(a, b) = unzip $ map transformR $ map (zipWith (*) window) $ framer frameSize (truncate overSamp) x
	(c, d) = unzip $ (shiftPitch frameSize pitchScale) (zip a (trueBinList3 overSamp b))

bigThing2 :: WAVE -> WAVE
bigThing2 audio = restoreSamples finishedSamples audio
	where
	doubs =  (samplesToDoubles . getSamples) audio
	finishedSamples = (doublesToSamples . (doublesToDoubles frameSize overSampD freqMultiplier)) doubs

{-
framesToFrames :: Int -> Double -> Double -> [[Double]] -> [[Double]]
framesToFrames frameSize overSamp pitchScale x = map invTransformR (zip c (invTrueBinList overSamp d))
	where
	(a, b) = unzip $ (map transformR x)
	(c, d) = unzip $ (map ((makeRightSize (div frameSize 2 + 1)) . (shiftPitch pitchScale))) (zip a (trueBinList overSamp b))

bigThing :: WAVE -> WAVE
bigThing audio = restoreSamples finishedSamples audio
	where
	frameList =  ((framer frameSize overSamp) . samplesToDoubles . getSamples) audio
	finishedSamples = (doublesToSamples . (deframer frameSize overSamp) . (framesToFrames frameSize overSampD freqMultiplier) . (map (zipWith (*) window))) frameList
-}




main = do
	args <- System.getArgs
	audio <- getWAVEFile (args !! 0)
	let b = bigThing2 audio
	putWAVEFile (args !! 1) b
