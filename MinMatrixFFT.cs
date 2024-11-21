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

        // Input the vector x (1xN), representing Wti at different times i
        double[] x = new double[n];
        Console.WriteLine("Enter the elements of vector x (Wti values, separated by spaces):");
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
            Console.Write("{time:F2} ");
        }
        Console.WriteLine();

        ///<summary>
        /// Print the Wiener process at each time ti
        /// </summary>
        Console.WriteLine("The Wiener process at each time ti is: ");
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
        Console.WriteLine("\nMa trận inverse:");
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

    public MinMatrixFFT(double unitinterval_ = 0.1, int bounder_ = 25) {
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

    public double[] ComputeGaussianFunction(double[] grid)
    {
        double sqrt2Pi = Math.Sqrt(2.0 * Math.PI);
        double[] gaussian = new double[grid.Length];
        for (int i = 0; i < grid.Length; i++)
        {
            gaussian[i] = Math.Exp(-0.5 * grid[i] * grid[i]) / sqrt2Pi;
        }
        return gaussian;
    }

    public Complex[] ConvolveWithFFT(double[] f1, double[] f2)
    {
        Contract.Requires(f1.Length == f2.Length);
        int nbElements = f1.Length;
        Complex[] fft1 = new Complex[nbElements];
        Complex[] fft2 = new Complex[nbElements];
        for (int i = 0; i < nbElements; i++)
        {
            fft1[i] = new Complex(f1[i], 0);
            fft2[i] = new Complex(f2[i], 0);
        }

        Fourier.Forward(fft1, FourierOptions.Matlab);
        Fourier.Forward(fft2, FourierOptions.Matlab);

        Complex[] h = new Complex[nbElements];
        for (int i = 0; i < nbElements; i++)
        {
            h[i] = fft1[i] * fft2[i];
        }

        Fourier.Inverse(h, FourierOptions.Matlab);
        return h;
    }
    public double[] LinearInterpolation(double[] grid, double[] h, double[] xValues)
    {
        /// <summary>
        ///     Linear Interpolation. Method: Connect two consecutive points with a straight line 
        ///     in the given points to interpolate the function. 
        /// </summary>
        /// <param name="grid">
        ///     Points in the x-axis of the given points.
        /// </param>
        /// <param name="h">
        ///     Points in the y-axis of the given points.
        /// </param>
        /// <param name="xValues">
        ///     Points in the x-axis to intepolate.
        /// </param>
        double[] fi = new double[xValues.Length];

        for (int j = 0; j < xValues.Length; j++)
        {
            double x = xValues[j];
            int i = 0;

            // Find the 2 nearest grid points grid[i] <= x <= grid[i+1]
            while (i < grid.Length - 1 && grid[i + 1] < x)
            {
                i++;
            }

            // If x is outside the grid length, return the boundary
            if (i == grid.Length - 1)
            {
                fi[j] = h[i];
            }
            else
            {
                // Do Linear interpolation
                fi[j] = h[i] + (x - grid[i]) * (grid[i + 1] - grid[i]) * (1 / (h[i + 1] - h[i]));
            }
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
        // Step 1: Create grid
        double[] grid = CreateCalculationGrid();

        // Step 2: Compute Gaussian function
        double[] gaussian = ComputeGaussianFunction(grid);

        // Step 3: Perform convolution
        Complex[] hi = ConvolveWithFFT(gaussian, gaussian);

        // Step 4: Interpolate function fi
        double[] hi_real = new double[hi.Length];
        for (int i = 0; i < hi.Length; i++) 
        {
            hi_real[i] = hi[i].Real;
        }         
        double[] fi = LinearInterpolation(grid, hi_real, grid);

        // Step 5: Compute the probability
        return ComputeProbability(fi, grid);
    }
}

class Program
{
    static void Main(string[] args)
    {
        MinMatrixFFT minMatrixFFT = new MinMatrixFFT(0.1, 25);



        //double[] input_1 = { 0.4, 0.5, 0.9, 0.2 };
        //double[] input_2 = { 0.9, 0.2, 0.5, 1.8 };
        //Complex[] result = minMatrixFFT.ConvolveWithFFT(input_1, input_2);


        //for (int i = 0; i < result.Length; i++)
        //    Console.WriteLine(result[i].Real);
        //for (int i = 0; i<result.Length; i++)
        //    Console.WriteLine(result[i]);

        double probability = minMatrixFFT.CalculateMinMatrixFFT();
        Console.WriteLine("Computed probability: " + probability);
    }
}

