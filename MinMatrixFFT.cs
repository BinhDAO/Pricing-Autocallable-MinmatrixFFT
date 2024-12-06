using System;
using System.Diagnostics.Contracts;
using System.Numerics;
using MathNet.Numerics.IntegralTransforms; // Using MathNet.Numerics for FFT operations
using MathNet.Numerics;
using System.Linq;
using MathNet.Numerics.Distributions;


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
                DateTime[] valuationDates = input
                    .Split(',')
                    .Select(date => DateTime.Parse(date.Trim()))
                    .ToArray();

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

                // Convert DateTime to double based on the time difference from initialDate
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
                double[] inputUpperBound = uInput
                    .Split(',')
                    .Select(val => double.Parse(val.Trim()))
                    .ToArray();

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

        public static List<double> originalFunction(double[] Variable, int i, int n)
        {
            if (i == 1)
            {
                return new List<double> { 1.0 };
            }
            else if (i > 1 && i <= n)
            {
                throw new InvalidOperationException("we are estimating f_i(y_i)");
            }
            else
            {
                throw new ArgumentOutOfRangeException("Index i is out of range");
            }
        }
        /// <summary>
        /// RegularGrid function to generate the grid for tempVar
        /// </summary>
        /// <param name="lowerBound"></param>
        /// <param name="upperBound"></param>
        /// <param name="unitInterval"></param>
        /// <returns></returns>
        // RegularGrid function to generate the grid for tempVar
        public static List<double> RegularGrid(double[] lowerBound, double[] upperBound, double unitInterval, double scale)
        {
            List<double> regularGrid = new List<double>();
            for (int i = 0; i < lowerBound.Length; i++)
            {
                double gridLowerBound = Math.Floor(lowerBound[i] * scale);
                double gridUpperBound = Math.Ceiling((upperBound[i] + 1) * scale);
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
        public List<double> tempFunction(double[] regularGrid, int i, int n, int p, int q, double unitInterval, int gridLowerBound, int gridUpperBound, double[] upperBound, double[] lowerBound)
        {
            List<double> tempFunction = new List<double>(regularGrid.Length);

            for (p = gridLowerBound; p <= gridUpperBound; p++)
            {
                tempFunction.Add(tempFunction[p]);
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
                ///Do the sum of all values convolution[q] from sumLowerBound to sumUpperBound
                ///</summary>
                for (q = sumLowerBound[i - 1]; q <= sumUpperBound[i - 1]; q++)
                {
                    ///summary
                    ///Do FFT for values of tempFunction at each point q in the sumInterval
                    ///</summary>
                    int fftSize = (int)Math.Pow(2, Math.Ceiling(Math.Log2(regularGrid.Length)));
                    Complex[] normalDensityValues = new Complex[fftSize];
                    Complex[] originalFunctionValues = new Complex[fftSize];

                    normalDensityValues[q] = new Complex(NormalDensity((p - q) * unitInterval), 0);
                    originalFunctionValues[i - 1] = new Complex(originalFunction(new double[] { q * unitInterval }, i - 1, n)[0], 0);

                    // Do FFT

                    Fourier.Forward(normalDensityValues, FourierOptions.Matlab);
                    Fourier.Forward(originalFunctionValues, FourierOptions.Matlab);

                    // Calculate the convolution in the frequency domain
                    Complex[] convolutionResult = new Complex[fftSize];
                    convolutionResult[q] = normalDensityValues[q] * originalFunctionValues[i - 1];

                    // Inverse of the FFT
                    Fourier.Inverse(convolutionResult, FourierOptions.Matlab);

                    // Sum of all the convolution[q] values to obtain tempFunction[p]
                    tempFunction[p] += convolutionResult[q].Real;
                }

            }
            return tempFunction;
        }

        /// <summary>
        /// Perform linear interpolation of all tempFunction values in the regulargrid to obtain originalFunction values in the targetGrid
        ///     Linear Interpolation. Method: Connect two consecutive points with a straight line 
        ///     in the given points to interpolate the function. 
        /// </summary>
        ///     Points in the x-axis to intepolate.
        /// </param>
        // Create targetGrid for the originalFunction
        public static List<double> TargetGrid(double lowerBound, double upperBound, double unitInterval, double scale)
        {
            // Adjust the lower and upper bounds using floor and ceil
            double originalGridLowerBound = Math.Floor(lowerBound + 1);
            double originalGridUpperBound = Math.Ceiling(upperBound);

            // Declare the number of points in the grid
            int orignialNumGridPoints = (int)((originalGridUpperBound - originalGridLowerBound) / unitInterval) + 1;

            // Create a list to hold the grid points
            List<double> targetGrid = new List<double>();

            // Generate the grid points
            for (int i = 0; i < orignialNumGridPoints; i++)
            {
                targetGrid.Add(originalGridLowerBound + i * unitInterval);
            }
            return targetGrid;
        }
        /// <summary>
        /// Interpolate f_i(y_i) from h_i(z) using linear interpolation
        /// </summary>
        /// <param name="targetGrid"></param>
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
        public static double[] originalFunctionInterpolation(List<double> targetGrid, List<double> tempFunction, List<double> originalFunction, double unitInterval, double[] Scale, int originalGridLowerBound, int originalGridUpperBound, int r, int i, int n)
        {
            for (r = originalGridLowerBound; r < originalGridUpperBound; r++)
            {
                originalFunction[r] = tempFunction[(int)Math.Floor(r * Scale[i])] + (((Scale[i] * r) - Math.Floor(r * Scale[i])) * ((tempFunction[(int)Math.Ceiling(Scale[i] * r)] - tempFunction[(int)Math.Ceiling(Scale[i] * r)])));
            }
            return originalFunction.ToArray();

        }
        public static double CumulativeDistributionFunction(double[] upperBound)
        {
            const double NEGATIVE_INFTY = double.NegativeInfinity;
            double[] negativeInfty = new double[upperBound.Length];
            for (int i = 0; i < upperBound.Length; i++)
            {
                negativeInfty[i] = NEGATIVE_INFTY;
            }
            return CumulativeDistributionFunction(negativeInfty, upperBound);
        }
        public static double CumulativeDistributionFunction(double[] upperBound, double[] lowerBound, double unitInterval)
        {
            throw new NotImplementedException();
        }
        public double CalculateMinMatrixFFT()
        // Create grid
        {
            // Create grid
            double[] grid = RegularGrid();

            // Initialize h_0(z) as Gaussian function on the grid
            double[] tempFunction_0 = new double[grid.Length];
            for (int i = 0; i < grid.Length; i++)
            {
                tempFunction_0[i] = NormalDensity(grid[i]);
            }
            // Initial value for h_i(z) and f_i(y_i)

            double[] tempFunction_i = tempFunction_0;
            double[] originalFunction_i = null;
            

            
        }
