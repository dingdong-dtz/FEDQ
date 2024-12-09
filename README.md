The main functions of running the program are as follows:

'FEDQ_theory'  folder is used to calculate the balance of parsing and verify the accuracy of the program

FDEQ_theory_sample.m:Calculation verification of analytical equilibrium examples

convergence.m:Calculation of program convergence

Other functional functions:

FDEQ_theory.m	FDEQ_theory_sample In the form of a function

evalJs.m:Interpolate the current voltage gradient onto the grid points

Gsolver_SOL_steed. m: Solving difference equations

Metric_SOL_steady. m: Solve the metric

Nabla2d'all_steed. m: Constructing difference equations

UpRZ_SOL_old. m: Update magnetic surface coordinate system grid

plot~~~:to plot figures

Theory_stample-m: The magnetic flux function used in the example

Diff_steady_all. m: A stable differential format

FdcoeffF.m: Ordinary differential format

The explanation of other small functions used in the local area can refer to the explanation of the specific code location used



'FEDQ'  folder is used to calculate the magnetic flux results on a rough rectangular grid and obtain accurate results on the magnetic surface coordinate system grid

FDEQ_sample.m:The main function of the example

change_saperatrix:Application examples of changing attack points

other functins are same to the functions in 'FEDQ_theory' folder

