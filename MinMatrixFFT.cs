using System;
using System.Diagnostics.Contracts;
using System.Numerics;
using MathNet.Numerics.IntegralTransforms; // Using MathNet.Numerics for FFT operations
using MathNet.Numerics;
using System.Linq;
using MathNet.Numerics.Distributions;
using System.Reflection.Metadata.Ecma335;
using static System.Formats.Asn1.AsnWriter;

namespace MinMatrixFFT
{
    class MinMatrixMethod
    {
        ///summary
        ///Declare upperBound and lowerBound vectors
        ///</summary>

        static void Main(string[] args)
        {
            /// <summary>
            /// declare time array
            /// </summary> 
            try
            {
                // Declare timeInterval vector
                List<DateTime> timeInterval = new List<DateTime>();

                // Input initial date
                Console.Write("Enter the initial date (yyyy-mm-dd): ");
                DateTime initialDate = DateTime.Parse(Console.ReadLine());
                timeInterval.Add(initialDate);

                // Input maturity T
                Console.Write("Enter the maturity date (yyyy-mm-dd): ");
                DateTime maturity = DateTime.Parse(Console.ReadLine());

                if (maturity <= DateTime.MinValue)
                {
                    Console.WriteLine("Maturity must be a future date");
                    return;
                }
                // Validate maturity date
                if (maturity <= initialDate)
                {
                    Console.WriteLine("Maturity date must be after the initial date.");
                    return;
                }

                // Input valuation dates (T1, T2, ..., Tn)
                Console.Write("Enter the valuation dates separated by commas (yyyy-mm-dd): ");
                string input = Console.ReadLine();
                DateTime[] valuationDates = input.Split(',').Select(date => DateTime.Parse(date.Trim())).ToArray();

                // Check valuation dates are before maturity and after the initial date
                if (valuationDates.Any(date => date < initialDate || date > maturity))
                {
                    Console.WriteLine("All valuation dates must be valid and before the maturity date.");
                    return;
                }
                // Add valuation dates to timeInterval
                timeInterval.AddRange(valuationDates);

                // Add maturity date to timeInterval
                timeInterval.Add(maturity);

                // Sort the timeInterval vector
                timeInterval = timeInterval.OrderBy(date => date).ToList();


                // declare the tempTimeVar d_i
                int n = timeInterval.Count - 1;
                double[] tempTimeVar = new double[n];

                for (int i = 1; i <= n; i++)
                {
                    // Calculate time difference between consecutive dates in days
                    double timeDifference = (timeInterval[i] - timeInterval[i - 1]).TotalDays;

                    // Calculate d_i as the reciprocal of the time difference
                    tempTimeVar[i] = 1.0 / timeDifference;
                }

                //declare (1xn) vector variable y_i
                double[] Variable = new double[n];

                //define (1x(n-1)) vector Scale 
                double[] Scale = new double[n - 1];
                for (int i = 1; i < n; i++)                         // Start from i = 1 because delta_i starts at i = 2
                {
                    // Calculate (1xn) vector delta_i 
                    Scale[i - 1] = Math.Sqrt(tempTimeVar[n - i] / tempTimeVar[n - i - 1]);
                }

                double[] upperBound = new double[n];
                double[] lowerBound = new double[n];

                //input the inputUpperBound vector
                Console.Write("Enter the values for inputUpperBound separated by commas (u1, u2, ..., un): ");
                string uInput = Console.ReadLine();
                double[] inputUpperBound = uInput.Split(',').Select(val => double.Parse(val.Trim())).ToArray();

                // Check if the length of inputUpperBound matches the number of intervals
                if (inputUpperBound.Length != n)
                {
                    Console.WriteLine("The length of inputUpperBound must be equal to the number of intervals (n)");
                    return;
                }

                /// Calculate the upperBound vector (b_i)
                for (int i = 0; i < n; i++)
                {
                    // Calculate b_i 
                    upperBound[i] = inputUpperBound[n - i - 1] * Math.Sqrt(tempTimeVar[n - i - 1]);
                }

                ///input the l vector (inputLowerBound)
                Console.Write("Enter the values for inputLowerBound separated by commas (l1, l2, ..., ln): ");
                string lInput = Console.ReadLine();
                double[] inputLowerBound = lInput.Split(',').Select(val => double.Parse(val.Trim())).ToArray();

                /// Check if the length of inputLowerBound matches the number of intervals
                if (inputLowerBound.Length != n)
                {
                    Console.WriteLine("The length of inputLowerBound must be equal to the number of intervals.");
                    return;
                }
                /// Calculate the lowerBound vector (a_i) using inputLowerBound
                for (int i = 0; i < n; i++)
                {
                    // Calculate a_i based on the formula a_i = l_{n+1-i} * sqrt(d_{n+1-i})
                    lowerBound[i] = inputLowerBound[n - i - 1] * Math.Sqrt(tempTimeVar[n - i - 1]);
                }

                // Input the unitInterval value (length of each divided unit of the integrating interval)
                Console.Write("Enter the unit interval length: ");
                double unitInterval = double.Parse(Console.ReadLine());

                // Validate that the unitInterval is a positive number
                if (unitInterval <= 0)
                {
                    Console.WriteLine("Unit interval must be a positive number.");
                    return;
                }


                ///summary
                ///Declare the tempVar z
                ///</summary>
                double[] tempVar = new double[n];

                for (int i = 0; i < n; i++)
                {
                    tempVar[i] = Variable[i] * Scale[i];
                }

            }
            catch (Exception ex)
            {
                Console.WriteLine($"An error occurred: {ex.Message}");
            }
        }

        ///summary
        ///Declare the originalFunction with initial value f_1(y_1) = 1
        ///</summary>

        //public static List<List<double>> originalFunction(int i, int n)
        //{
        //        originalFunction[1].Add(1.0);


        //    originalFunction.Add(originalFunction[i]);


        //    return originalFunction;
        //}
        /// <summary>
        /// RegularGrid function to generate the grid for tempVar
        /// </summary>
        /// <param name="lowerBound"></param>
        /// <param name="upperBound"></param>
        /// <param name="unitInterval"></param>
        /// <returns></returns>
        // RegularGrid function to generate the grid for tempVar
        public static List<double> RegularGrid(double[] lowerBound, double[] upperBound, double unitInterval, double[] Scale)
        {
            List<double> regularGrid = new List<double>();
            for (int i = 0; i < lowerBound.Length; i++)
            {
                double gridLowerBound = Math.Floor(lowerBound[i] * Scale[i]);
                double gridUpperBound = Math.Ceiling((upperBound[i] + 1) * Scale[i]);
                int numGridPoints = (int)((gridUpperBound - gridLowerBound) / unitInterval) + 1;

                for (int j = 0; j < numGridPoints; j++)
                {
                    regularGrid.Add(gridLowerBound + j * unitInterval);
                }
            }
            return regularGrid;
        }



        // Normal density function
        public static double NormalDensity(double x)                                // LN: NormalDensity
        {
            return (1.0 / Math.Sqrt(2 * Math.PI)) * Math.Exp(-0.5 * (x * x));
        }

        /// <summary>
        /// Compute tempFunction at all grid points by using the sum of all convolutions, each convolution is calculated by FFT method 
        /// </summary>
        /// <param name="regularGrid"></param>
        /// <param name="originalFunction"></param>
        /// <param name="unitInterval"></param>
        /// <param name="lowerBound"></param>
        /// <param name="upperBound"></param>
        /// <returns></returns>
        public static List<double> tempFunctionMethod(List<double> regularGrid, int i, int n, int p, int q, double unitInterval, int gridLowerBound, int gridUpperBound, double[] upperBound, double[] lowerBound)
        {
            List<List<double>> tempFunction = new List<List<double>>(n);
            List<List<double>> originalFunction = new List<List<double>>(n);

            originalFunction.Add(new List<double> { 1.0 });
            for (i = 0; i < n; i++)
            {
                tempFunction.Add(new List<double>(gridUpperBound - gridLowerBound + 1));
            }
            for (i = 1; i < n; i++)
            {
                // Define sum upperbound ceil(upperBound_(i-1) 
                int[] sumUpperBound = new int[upperBound.Length];
                for (i = 1; i < upperBound.Length; i++)
                {
                    sumUpperBound[i] = (int)Math.Ceiling(upperBound[i - 1]);
                }
                // Define sum lowerbound floor((lowerBound_(i-1) + 1) 
                int[] sumLowerBound = new int[lowerBound.Length];
                for (i = 1; i < upperBound.Length; i++)
                {
                    sumLowerBound[i] = (int)Math.Floor(lowerBound[i - 1]);
                }
                ///summary
                /// loop through the grid to add tempFunction values at each p
                ///</summary>
                double tempFunctionSum = 0.0;

                for (p = gridLowerBound; p <= gridUpperBound; p++)
                {
                    for (q = sumLowerBound[i - 1]; q <= sumUpperBound[i - 1]; q++)
                    {
                        double normalDensityValue = NormalDensity((p - q) * unitInterval);
                        ///summary
                        ///Do FFT for values of tempFunction_i
                        ///</summary>
                        int fftSize = (int)Math.Pow(2, Math.Ceiling(Math.Log2(regularGrid.Count)));
                        Complex[] normalDensityComplex = new Complex[fftSize];
                        Complex[] originalFunctionComplex = new Complex[fftSize];
                        normalDensityComplex[q] = new Complex(NormalDensity((p - q) * unitInterval), 0);
                        originalFunctionComplex[i - 1] = new Complex(originalFunction[i - 1][q], 0);
                        // Do FFT

                        Fourier.Forward(normalDensityComplex, FourierOptions.Matlab);
                        Fourier.Forward(originalFunctionComplex, FourierOptions.Matlab);

                        // Calculate the convolution in the frequency domain
                        Complex[] convolutionResult = new Complex[fftSize];
                        convolutionResult[q] = normalDensityComplex[q] * originalFunctionComplex[i - 1];

                        // Inverse of the FFT
                        Fourier.Inverse(convolutionResult, FourierOptions.Matlab);

                        // Sum of all the convolution[q] values to obtain the tempFunction values
                        tempFunctionSum += convolutionResult[q].Real;
                    }
                    tempFunction[i].Add(tempFunctionSum);
                }
            }
            return tempFunction[i];
        }

        /// <summary>
        /// Perform linear interpolation of all tempFunction values in the regulargrid to obtain originalFunction values in the targetGrid
        ///     Linear Interpolation. Method: Connect two consecutive points with a straight line 
        ///     in the given points to interpolate the function. 
        /// </summary>
        ///     Points in the x-axis to intepolate.
        /// </param>
        // Create targetGrid for the originalFunction
        public static List<double> originalGrid(double lowerBound, double upperBound, double unitInterval, double scale)
        {
            // Adjust the lower and upper bounds using floor and ceil
            double originalGridLowerBound = Math.Floor(lowerBound + 1);
            double originalGridUpperBound = Math.Ceiling(upperBound);

            // Declare the number of points in the grid
            int orignialNumGridPoints = (int)((originalGridUpperBound - originalGridLowerBound) / unitInterval) + 1;

            // Create a list to hold the grid points
            List<double> originalGrid = new List<double>();

            // Generate the grid points
            for (int i = 0; i < orignialNumGridPoints; i++)
            {
                originalGrid.Add(originalGridLowerBound + i * unitInterval);
            }
            return originalGrid;
        }
        /// <summary>
        /// Interpolate f_i(y_i) from h_i(z) using linear interpolation
        /// </summary>
        /// <param name="originalGrid"></param>
        /// <param name="regularGrid"></param>
        /// <param name="tempFunction"></param>
        /// <param name="originalFunction"></param>
        /// <param name="unitInterval"></param>
        /// <param name="gridLowerBound"></param>
        /// <param name="gridUpperBound"></param>
        /// <param name="originalGridLowerBound"></param>
        /// <param name="originalGridUpperBound"></param>
        /// <param name="i"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static List<double> originalFunctionInterpolation(List<double> originalGrid, List<List<double>> tempFunction, List<List<double>> originalFunction, double unitInterval, double[] Scale, int originalGridLowerBound, int originalGridUpperBound, int r, int i, int n)
        {
            for (i = 1; i <= n; i++)
            {
                for (r = originalGridLowerBound; r < originalGridUpperBound; r++)
                {
                    originalFunction[i][r] = tempFunction[i][(int)Math.Floor(r * Scale[i])] + (((Scale[i] * r) - Math.Floor(r * Scale[i])) * ((tempFunction[i][(int)Math.Ceiling(Scale[i] * r)] - tempFunction[i][(int)Math.Ceiling(Scale[i] * r)])));
                }

            }
            return originalFunction[i];

        }
        public static double CalculateProbability(int n, double[] lowerBound, double[] upperBound, double unitInterval, List<DateTime> timeInterval, DateTime maturity, DateTime[] valuationDates)
        {
            try
            {
                // Validate that the timeInterval contains at least n+1 dates
                if (timeInterval.Count != n + 1)
                {
                    throw new ArgumentException("The timeInterval must contain exactly n+1 dates.");
                }

                // Compute tempTimeVar (d_i) 
                double[] tempTimeVar = new double[n];
                for (int i = 1; i <= n; i++)
                {
                    tempTimeVar[i - 1] = 1.0 / (timeInterval[i] - timeInterval[i - 1]).TotalDays;
                }

                // Compute Scale array 
                double[] Scale = new double[n - 1];
                for (int i = 1; i < n; i++)
                {
                    Scale[i - 1] = Math.Sqrt(tempTimeVar[n - i] / tempTimeVar[n - i - 1]);
                }

                List<double> regularGrid = RegularGrid(lowerBound, upperBound, unitInterval, Scale);
                List<List<double>> originalFunction = new List<List<double>>(n);
                originalFunction.Add(new List<double> { 1.0 });
                List<List<double>> tempFunction = new List<List<double>>(n);

                for (int i = 1; i < n; i++)
                {
                    tempFunction.Add(tempFunctionMethod(regularGrid, i, n, 1, 1, unitInterval, (int)Math.Floor(lowerBound[i]), (int)Math.Ceiling(upperBound[i]), upperBound, lowerBound));
                    originalFunction.Add(originalFunctionInterpolation(originalGrid(lowerBound[i], upperBound[i], unitInterval, Scale[i - 1]), tempFunction, originalFunction, unitInterval, Scale, (int)Math.Floor(lowerBound[i] + 1), (int)Math.Ceiling(upperBound[i]), 1, i, n));
                }

                /// Calculate the final probability by summing over the grid

                double probability = 0.0;

                for (int p = 1; p < regularGrid.Count; p++)
                {
                    double gridPoint = regularGrid[p];
                    double normalDensityValue = NormalDensity(gridPoint);
                    double originalFunctionValue = originalFunction[n - 1][p];
                    probability += normalDensityValue * originalFunctionValue;
                }

                // Multiply by unitInterval to calculate the final probability
                probability *= unitInterval;

                return probability;
            }
            catch (Exception ex)
            {
                Console.WriteLine($"An error occurred: {ex.Message}");
                return 0.0;
            }
        }

        //public static double CumulativeDistributionFunction(double[] upperBound)
        //{
        //    const double NEGATIVE_INFTY = double.NegativeInfinity;
        //    double[] negativeInfty = new double[upperBound.Length];
        //    for (int i = 0; i < upperBound.Length; i++)
        //    {
        //        negativeInfty[i] = NEGATIVE_INFTY;
        //    }
        //    return CumulativeDistributionFunction(negativeInfty, upperBound);
        //}
        //public static double CumulativeDistributionFunction(double[] upperBound, List<double> tempFunction, double[] lowerBound, double unitInterval, List<double> regularGrid)
        //{
        //    throw new NotImplementedException();
        //    {
        //        tempFunction.ToArray() = originalFunction 
        //    }
        //}
        //    public static double CalculateProbability(int i, int n, double[] lowerBound, double[] upperBound, double[] Scale, double unitInterval, List<DateTime> timeInterval, DateTime maturity, DateTime[] valuationDates)
        //    {
        //        List<List<double>> originalFunction = new List<List<double>>();
        //        List<List<double>> tempFunction = new List<List<double>>();
        //        List<double> regularGrid = RegularGrid(lowerBound, upperBound, unitInterval, Scale);

        //        // Initialize the first element of originalFunction with 1.0
        //        originalFunction.Add(new List<double> { 1.0 });


        //        originalFunction.Add(new List<double> { 1.0 });

        //        for (i = 1; i < n; i++)
        //        {
        //            tempFunction.Add(tempFunctionMethod(regularGrid, originalFunction, unitInterval, lowerBound, upperBound));
        //            originalFunction.Add(originalFunctionInterpolation(regularGrid, tempFunction, originalGrid, originalFunction, unitInterval, Scale));
        //        }

        //        double probability = 0.0;

        //        // Generate the grid for summing

        //        for (p = 0; p < regularGrid.Count; p++)
        //        {
        //            double gridPoint = regularGrid[p] * unitInterval;
        //            double normalDensityValue = NormalDensity(gridPoint);
        //            double originalFunctionValue = originalFunction[n][p];

        //            // Add to probability
        //            probability += normalDensityValue * originalFunctionValue;
        //        }

        //        // Step 4: Multiply by unitInterval to calculate final probability
        //        probability *= unitInterval;

        //        return probability;
        //    }
        //}

        //public static double CumulativeDistributionFunction(double[] lowerBound, double[] higherBound)
        //{
        //    //LN: to implement
        //    throw new NotImplementedException();
        //}
        //public double CalculateMinMatrixFFT()
        //{
        //    // Create grid
        //    double[] grid = CreateCalculationGrid();

        //    // Initialize h_0(z) as Gaussian function on the grid
        //    double[] h0 = new double[grid.Length];
        //    for (int i = 0; i < grid.Length; i++)
        //    {
        //        h0[i] = Gaussian(grid[i]);
        //    }

        //    // Initial value for h_i(z) and f_i(y_i)
        //    double[] hi = h0;
        //    double[] fi = null;

        //    // Initial value for delta array 
        //    double[] delta = new double[bounder];
        //    for (int i = 0; i < delta.Length; i++)
        //    {
        //        delta[i] = 1.0;
        //    }
        //    // Iteratively compute h_i(z) and f_i(y_i) for each i
        //    for (int i = 0; i < delta.Length; i++)
        //    {
        //        // Compute h_i(z) from h_{i-1}(z) using FFT convolution
        //        hi = ComputeHi(grid, hi, delta);

        //        // Interpolate f_i(y_i) from h_i(z)
        //        fi = InterpolateFi(grid, hi, delta[i], unitInterval, grid.Length);
        //    }
        //    // Compute the probability (sum over fi with integration)
        //    double probability = ComputeProbability(fi, grid);
        //    // Return the final probability value
        //    return probability;
        //}
        //    }
        //double[] result = new double[regularGrid.Length];
        //for (int q = 0; q < regularGrid.Length; q++)
        //{
        //    if (regularGrid[q] >= lowerBound && regularGrid[q] <= upperBound)
        //    {
        //    }
        //}

        //return result;
    }
}


//// Compute y_i from r and unitInterval
//double yi = r * unitInterval;

//// Turn y_i into z with z = delta[i] * y_i
//double z = delta * yi;

//// find grid index that correspond to z
//int lowerIndex = (int)Math.Floor(z / unitInterval);
//int upperIndex = (int)Math.Ceiling(z / unitInterval);

//// Ensure that index is on the grid
//if (lowerIndex < 0 || upperIndex >= grid.Length)
//{
//    fi[r] = 0.0; // if outside the bounder, return 0
//    continue;
//}

//// get h_i(z) values at lowerIndex and upperIndex
//double hLower = hi[lowerIndex];
//double hUpper = hi[upperIndex];

//// Compute interpolation parameter
//double weight = (z - grid[lowerIndex]) / (grid[upperIndex] - grid[lowerIndex]);

//// Interpolation
//fi[r] = hLower + weight * (hUpper - hLower);

/// <summary>
/// Declare the tempFunction h_i(z)
/// </summary>
/// <param name="t"></param>

/// <summary>
/// Covariance matrix computing function
/// </summary>
/// <param name="t"></param>
/// <returns></returns>
//        static double[,] ComputeCovarianceMatrix(double[] t)
//        {
//            int n = t.Length;
//            double[,] covarianceMatrix = new double[n, n];
//            for (int i = 0; i < n; i++)
//            {
//                for (int j = 0; j < n; j++)
//                {
//                    covarianceMatrix[i, j] = Math.Min(t[i], t[j]);
//                }
//            }
//            return covarianceMatrix;
//        }

//        ///<summary>
//        /// Create the covariance matrix (Minmatrix)
//        /// </summary>
//        double[,] covarianceMarix = ComputeCovarianceMatrix(x);
//                foreach (var time in x)
//                {
//                    Console.Write($"{time:F2} ");
//                }
//    Console.WriteLine();

//                ///<summary>
//                /// Print the Wiener process at each time ti
//                /// </summary>
//                Console.WriteLine("The Wiener process at each time t_i is: ");
//                foreach (var time in x)
//                {
//                    Console.Write("{time:F2} ");
//                }
//Console.WriteLine();

////Print the covariance matrix
//Console.WriteLine("The covariance matrix is: ");
//for (int i = 0; i < n; i++)
//{
//    for (int j = 0; j < n; j++)
//    {
//        Console.Write("{ covarianceMatrix[i, j]:F2} ");
//    }
//    Console.WriteLine();
//}

//// Compute d_k
//double[] d = new double[n - 1];
//for (int k = 1; k < n; k++)                                           // d[1] is Wt2 - Wt1
//{
//    d[k - 1] = 1 / (x[k] - x[k - 1]);
//}

//// Create the tridiagonal inverse matrix
//double[,] inverseMatrix = new double[n, n];

//for (int i = 0; i < n; i++)
//{
//    inverseMatrix[i, i] = i == n - 1 ? d[i] : d[i] + d[i + 1];
//    if (i > 0)
//        inverseMatrix[i, i - 1] = -d[i];                    // a_{k, k+1}
//    if (i < n - 1)
//        inverseMatrix[i, i + 1] = -d[i];                    // a_{k+1, k}            
//}
/////<summary>
///// Print the inverse matrix
///// </summary>
//Console.WriteLine("Inverse matrix is:");
//for (int i = 0; i < n; i++)
//{
//    for (int j = 0; j < n; j++)
//    {
//        Console.Write($"{inverseMatrix[i, j]:F2} ");
//    }
//    Console.WriteLine();
//}

//// Input the vector l (lower bound multipliers)
//double[] l = new double[n];
//Console.WriteLine("Enter the elements of vector l (lower bound multipliers, separated by spaces):");
//string[] lInput = Console.ReadLine().Split();
//for (int i = 0; i < n; i++)
//{
//    l[i] = double.Parse(lInput[i]);
//}

//// Input the vector u (upper bound multipliers)
//double[] u = new double[n];
//Console.WriteLine("Enter the elements of vector u (upper bound multipliers, separated by spaces):");
//string[] uInput = Console.ReadLine().Split();
//for (int i = 0; i < n; i++)
//{
//    u[i] = double.Parse(uInput[i]);
//}

//// Create and compute the y vector
//double[] y = new double[n];
//for (int i = 0; i < n; i++)
//{
//    y[i] = Math.Sqrt(d[n - 1 - i]) * x[n - 1 - i]; // Change from x to y based on formula
//}

//// Print the y vector
//Console.WriteLine("The y vector is:");
//foreach (var yi in y)
//{
//    Console.Write($"{yi:F2} ");
//}
//Console.WriteLine();

//// Compute delta(i) for i = 1 to n-1
//double[] delta = new double[n - 1];
//for (int i = 0; i < n - 1; i++)
//{
//    delta[i] = Math.Sqrt(d[n - 2 - i] / d[n - 1 - i]);
//}
//            }



//        /// <summary>
//        /// Declare the unit interval and bounder for grid construction
//        /// </summary>
//        public double unitInterval { get; set; }
//public int bounder { get; set; }

//public MinMatrixMethod(double unitinterval_ = 0.005, int bounder_ = 25)
//{
//    unitInterval = unitinterval_;
//    bounder = bounder_;
//}

///// <summary>
///// create calculation grid
///// </summary>
///// <returns></returns>

//public double[] CreateCalculationGrid()// LN: RegularGrid
//{
//    int gridSize = (int)(2 * bounder / unitInterval) + 1;
//    double[] grid = new double[gridSize];
//    for (int i = 0; i < gridSize; i++)
//    {
//        grid[i] = -bounder + i * unitInterval;
//    }
//    return grid;

//}

//// Define Gaussian function
//public static double Gaussian(double x)// LN: NormalDensity
//{
//    return Math.Exp(-0.5 * x * x) / Math.Sqrt(2 * Math.PI);
//}

//// Compute h_i(z) using FFT convolution
//public static double[] ComputeHi(double[] grid, double[] hiMinus1, double[] delta) //LN: tempFunction
//{
//    int n = grid.Length;
//    Complex[] gaussian = new Complex[n]; // LN: normalDensity
//    Complex[] hiPrev = new Complex[n]; // LN: hPrev

//    // Prepare Gaussian function and h_{i-1}(y_{i-1}) for FFT
//    for (int i = 0; i < n; i++)
//    {
//        gaussian[i] = new Complex(Gaussian(grid[i] / delta[i]), 0);
//        hiPrev[i] = new Complex(hiMinus1[i], 0); // Previous h_{i-1}
//    }

//    // Apply FFT to both arrays
//    Fourier.Forward(gaussian, FourierOptions.Matlab);
//    Fourier.Forward(hiPrev, FourierOptions.Matlab);

//    // Pointwise multiplication in frequency domain
//    for (int i = 0; i < n; i++)
//    {
//        gaussian[i] *= hiPrev[i];
//    }

//    // Apply inverse FFT to get convolution result
//    Fourier.Inverse(gaussian, FourierOptions.Matlab);

//    // Normalize and take the real part as result
//    double[] hi = new double[n];
//    for (int i = 0; i < n; i++)
//    {
//        hi[i] = gaussian[i].Real / n;
//    }

//    return hi;
//}


////LN: New method
//public

//    
//public double ComputeProbability(double[] fi, double[] grid)
//{
//    double sum = 0;
//    for (int i = 0; i < fi.Length; i++)
//    {
//        sum += fi[i] * unitInterval;
//    }
//    return sum;
//}

//LN: Computation of Multi-dimensional Gaussian Repartition
//public static double CumulativeDistributionFunction(double[] higherBound)
//{
//    const double NEGATIVE_INFTY = -10defcef6e6edvtyss4dvy6e4dvyde.0;
//    double[] negativeInfty = new double[higherBound.Length];
//    for (int i = 0; i < higherBound.Length; i++)
//    {
//        negativeInfty[i] = NEGATIVE_INFTY;
//    }
//    return CumulativeDistributionFunction(negativeInfty, higherBound);
//}

//public static double CumulativeDistributionFunction(double[] lowerBound, double[] higherBound)
//{
//    //LN: to implement
//    throw new NotImplementedException();
//}
//public double CalculateMinMatrixFFT()
//{
//    // Create grid
//    double[] grid = CreateCalculationGrid();

//    // Initialize h_0(z) as Gaussian function on the grid
//    double[] h0 = new double[grid.Length];
//    for (int i = 0; i < grid.Length; i++)
//    {
//        h0[i] = Gaussian(grid[i]);
//    }

//    // Initial value for h_i(z) and f_i(y_i)
//    double[] hi = h0;
//    double[] fi = null;

//    // Initial value for delta array 
//    double[] delta = new double[bounder];
//    for (int i = 0; i < delta.Length; i++)
//    {
//        delta[i] = 1.0;
//    }
//    // Iteratively compute h_i(z) and f_i(y_i) for each i
//    for (int i = 0; i < delta.Length; i++)
//    {
//        // Compute h_i(z) from h_{i-1}(z) using FFT convolution
//        hi = ComputeHi(grid, hi, delta);

//        // Interpolate f_i(y_i) from h_i(z)
//        fi = InterpolateFi(grid, hi, delta[i], unitInterval, grid.Length);
//    }
//    // Compute the probability (sum over fi with integration)
//    double probability = ComputeProbability(fi, grid);
//    // Return the final probability value
//    return probability;
//}
//    }

//class Program
//{
//    static void Main(string[] args)
//    {
//        MinMatrixFFT minMatrixFFT = new MinMatrixFFT(0.1, 25);



//double[] input_1 = { 0.4, 0.5, 0.9, 0.2 };
//double[] input_2 = { 0.9, 0.2, 0.5, 1.8 };
//Complex[] result = minMatrixFFT.ConvolveWithFFT(input_1, input_2);


//for (int i = 0; i < result.Length; i++)
//    Console.WriteLine(result[i].Real);
//for (int i = 0; i<result.Length; i++)
//        //    Console.WriteLine(result[i]);

//        double probability = minMatrixFFT.CalculateMinMatrixFFT();
//        Console.WriteLine("Computed probability: " + probability);
//    }
//}
