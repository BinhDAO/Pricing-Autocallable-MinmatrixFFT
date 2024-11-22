using System;
using System.Diagnostics.Contracts;
using System.Numerics;
using MathNet.Numerics.IntegralTransforms; // Using MathNet.Numerics for FFT operations
using MathNet.Numerics;


public class MinMatrixFFT
{   /// <summary>
/// declare time array
/// </summary>
    static void Main(string[] args)
    {
        // Input the number of elements n
        Console.Write("Enter the number of elements n: ");
        int n = int.Parse(Console.ReadLine());

        // Input the vector x (1xN), representing W[ti] at different times i
        double[] x = new double[n];
        Console.WriteLine("Enter the elements of vector x (W[ti] values, separated by spaces): ");
        string[] xInput = Console.ReadLine().Split();
        for (int i = 0; i < n; i++)
        {
            x[i] = double.Parse(xInput[i]);
        }
        /// <summary>
        /// Covariance matrix computing function
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        static double[,] ComputeCovarianceMatrix(double[] t)
        {
            int n = t.Length;
            double[,] covarianceMatrix = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    covarianceMatrix[i, j] = Math.Min(t[i], t[j]);
                }
            }
            return covarianceMatrix;
        }

        ///<summary>
        /// Create the covariance matrix (Minmatrix)
        /// </summary>
        double[,] covarianceMarix = ComputeCovarianceMatrix(x);
        foreach (var time in x)
        {
            Console.Write("${time:F2} ");
        }
        Console.WriteLine();

        ///<summary>
        /// Print the Wiener process at each time ti
        /// </summary>
        Console.WriteLine("The Wiener process at each time t_i is: ");
        foreach (var time in x)
        {
            Console.Write("{time:F2} ");
        }
        Console.WriteLine();

        //Print the covariance matrix
        Console.WriteLine("The covariance matrix is: ");
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Console.Write("{ covarianceMatrix[i, j]:F2} ");
            }
            Console.WriteLine();
        }

        // Compute d_k
        double[] d = new double[n - 1];
        for (int k = 1; k < n; k++)                                           // d[1] is Wt2 - Wt1
        {
            d[k - 1] = 1 / (x[k] - x[k - 1]);
        }

        // Create the tridiagonal inverse matrix
        double[,] inverseMatrix = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)                                             // main diagonal
                {
                    if (i == n)
                        inverseMatrix[i, j] = d[n];                     // last element of the main diagonal is d_n + d_{n+1} = d_n
                    else
                        inverseMatrix[i, j] = d[i] + d[i + 1];          // d_k + d_{k+1}
                }
                else if (Math.Abs(i - j) == 1)                          // elements that are above and below the main diagonal
                {
                    if (i < j)
                        inverseMatrix[i, j] = -d[i];                    // a_{k, k+1}
                    else
                        inverseMatrix[i, j] = -d[j];                    // a_{k+1, k}
                }
                else                                                    // other elements outside the 3 diagonals are all 0
                {
                    inverseMatrix[i, j] = 0;
                }
            }
        }
        ///<summary>
        /// Print the inverse matrix
        /// </summary>
        Console.WriteLine("Inverse matrix is:");
        for (int i = 0; i < n; i++)
        {
            {
                for (int j = 0; j < n; j++)
                {
                    {
                        Console.Write($"{inverseMatrix[i, j]:F2} ");
                    }
                    Console.WriteLine();
                }
            }
        }

        // Input the vector l (lower bound multipliers)
        double[] l = new double[n];
        Console.WriteLine("Enter the elements of vector l (lower bound multipliers, separated by spaces):");
        string[] lInput = Console.ReadLine().Split();
        for (int i = 0; i < n; i++)
        {
            l[i] = double.Parse(lInput[i]);
        }

        // Input the vector u (upper bound multipliers)
        double[] u = new double[n];
        Console.WriteLine("Enter the elements of vector u (upper bound multipliers, separated by spaces):");
        string[] uInput = Console.ReadLine().Split();
        for (int i = 0; i < n; i++)
        {
            u[i] = double.Parse(uInput[i]);
        }

        // Create and compute the y vector
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            y[i] = Math.Sqrt(d[n - 1 - i]) * x[n - 1 - i]; // Change from x to y based on formula
        }

        // Print the y vector
        Console.WriteLine("The y vector is:");
        foreach (var yi in y)
        {
            Console.Write($"{yi:F2} ");
        }
        Console.WriteLine();

        // Compute delta(i) for i = 1 to n-1
        double[] delta = new double[n - 1];
        for (int i = 0; i < n - 1; i++)
        {
            delta[i] = Math.Sqrt(d[n - 2 - i] / d[n - 1 - i]);
        }
    }
          


/// <summary>
/// Declare the unit interval and bounder for grid construction
/// </summary>
public double unitInterval { get; set; }
    public int bounder { get; set; }

    public MinMatrixFFT(double unitinterval_ = 0.005, int bounder_ = 25) {
        unitInterval = unitinterval_;
        bounder = bounder_;
    }

    /// <summary>
    /// create calculation grid
    /// </summary>
    /// <returns></returns>

    public double[] CreateCalculationGrid()
    {
            int gridSize = (int) (2 * bounder / unitInterval) + 1;
            double[] grid = new double[gridSize];
            for (int i = 0; i < gridSize; i++)
            {
                grid[i] = -bounder + i * unitInterval;
            }
            return grid;
    
    }

    // Define Gaussian function
    public static double Gaussian(double x)
    {
        return Math.Exp(-0.5 * x * x) / Math.Sqrt(2 * Math.PI);
    }

    // Compute h_i(z) using FFT convolution
    public static double[] ComputeHi(double[] grid, double[] hiMinus1, double[] delta)
    {
        int n = grid.Length;
        Complex[] gaussian = new Complex[n];
        Complex[] hiPrev = new Complex[n];

        // Prepare Gaussian function and h_{i-1}(y_{i-1}) for FFT
        for (int i = 0; i < n; i++)
        {
            gaussian[i] = new Complex(Gaussian(grid[i] / delta[i]), 0); 
            hiPrev[i] = new Complex(hiMinus1[i], 0); // Previous h_{i-1}
        }

        // Apply FFT to both arrays
        Fourier.Forward(gaussian, FourierOptions.Matlab);
        Fourier.Forward(hiPrev, FourierOptions.Matlab);

        // Pointwise multiplication in frequency domain
        for (int i = 0; i < n; i++)
        {
            gaussian[i] *= hiPrev[i];
        }

        // Apply inverse FFT to get convolution result
        Fourier.Inverse(gaussian, FourierOptions.Matlab);

        // Normalize and take the real part as result
        double[] hi = new double[n];
        for (int i = 0; i < n; i++)
        {
            hi[i] = gaussian[i].Real / n;
        }

        return hi;
    }
    
     
        /// <summary>
        ///     Linear Interpolation. Method: Connect two consecutive points with a straight line 
        ///     in the given points to interpolate the function. 
        /// </summary>
        ///     Points in the x-axis to intepolate.
        /// </param>
        // Interpolate f_i(y_i) from h_i(z) using linear interpolation
        public static double[] InterpolateFi(double[] grid, double[] hi, double delta, double unitInterval, int n)
        {
            // number of y_i elements
            double[] fi = new double[n];

            for (int r = 0; r < n; r++)
            {
                // Compute y_i from r and unitInterval
                double yi = r * unitInterval;

                // Turn y_i into z with z = delta[i] * y_i
                double z = delta * yi;

                // find grid index that correspond to z
                int lowerIndex = (int)Math.Floor(z / unitInterval);
                int upperIndex = (int)Math.Ceiling(z / unitInterval);

                // Ensure that index is on the grid
                if (lowerIndex < 0 || upperIndex >= grid.Length)
                {
                    fi[r] = 0.0; // if outside the bounder, return 0
                    continue;
                }

                // get h_i(z) values at lowerIndex and upperIndex
                double hLower = hi[lowerIndex];
                double hUpper = hi[upperIndex];

                // Compute interpolation parameter
                double weight = (z - grid[lowerIndex]) / (grid[upperIndex] - grid[lowerIndex]);

                // Interpolation
                fi[r] = hLower + weight * (hUpper - hLower);
            }

            return fi;
        }

        public double ComputeProbability(double[] fi, double[] grid)
    {
        double sum = 0;
        for (int i = 0; i < fi.Length; i++)
        {
            sum += fi[i] * unitInterval;
        }
        return sum;
    }

    public double CalculateMinMatrixFFT()
    {
        // Create grid
        double[] grid = CreateCalculationGrid();

        // Initialize h_0(z) as Gaussian function on the grid
        double[] h0 = new double[grid.Length];
        for (int i = 0; i < grid.Length; i++)
        {
            h0[i] = Gaussian(grid[i]);
        }

        // Initial value for h_i(z) and f_i(y_i)
        double[] hi = h0;
        double[] fi = null;

        // Initial value for delta array 
        double[] delta = new double[bounder];
        for (int i = 0; i < delta.Length; i++)
        {
            delta[i] = 1.0;                             
        }
        // Iteratively compute h_i(z) and f_i(y_i) for each i
        for (int i = 0; i < delta.Length; i++)
        {
            // Compute h_i(z) from h_{i-1}(z) using FFT convolution
            hi = ComputeHi(grid, hi, delta);

            // Interpolate f_i(y_i) from h_i(z)
            fi = InterpolateFi(grid, hi, delta[i], unitInterval, grid.Length);
        }
        // Compute the probability (sum over fi with integration)
        double probability = ComputeProbability(fi, grid);
        // Return the final probability value
        return probability;
    }
}

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
