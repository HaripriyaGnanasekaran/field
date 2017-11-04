//Here should the input data come

// define system
// The size of the box in x y amd z direction.
int Mx = 60;
int My = 60;
int Mz = 60;		
int N_g = 25; 		// length of the polymer chain or the number of beads in a segment of the dendrimer in cluding the barnch point.
double f = 1; 		// the number of molecules in the box
double chi = 0.0;
const int dendrimer = 1;	// If the molecule simulated is a dendrimer it should be one
int generations = 4; 	// the number of generations incase the molecule is a dendrimer (it should be 1 in case of simple polymers)
const int dendfunc = 3;	// the number of branches in each node incase the molecule is a dendrimer
int nfixed = 1; 		// the first nfixed nodes are not moved by the MC algorithm.
const int endnode = 0; 	// If the endpoints of the dendrimer should also have 
double partrad = 1.0; 	// The size of a particle all lattice sites within this radius will be set in the mask
const int Npart =22;		// the number of nodes 
bool allowoverlap = false; // determines wether nodes are allowed to overlap
				
double H_poslist[Npart*3] = {30.5,30.5,30.5,40.5,30.5,30.5,20.5,30.5,30.5,30.5,40.5,30.5,40.5,20.5,30.5,40.5,30.5,40.5,20.5,30.5,20.5,10.5,30.5,30.5,30.5,50.5,30.5,30.5,40.5,40.5,40.5,20.5,20.5,50.5,20.5,30.5,30.5,30.5,40.5,40.5,40.5,40.5,20.5,20.5,20.5,20.5,30.5,10.5,10.5,40.5,30.5,10.5,20.5,30.5,30.5,50.5,40.5,30.5,50.5,20.5,20.5,40.5,40.5,30.5,40.5,50.5}; //Vector containing all particle positions.
double H_oldposlist[Npart*3] = {30.5,30.5,30.5,40.5,30.5,30.5,20.5,30.5,30.5,30.5,40.5,30.5,40.5,20.5,30.5,40.5,30.5,40.5,20.5,30.5,20.5,10.5,30.5,30.5,30.5,50.5,30.5,30.5,40.5,40.5,40.5,20.5,20.5,50.5,20.5,30.5,30.5,30.5,40.5,40.5,40.5,40.5,20.5,20.5,20.5,20.5,30.5,10.5,10.5,40.5,30.5,10.5,20.5,30.5,30.5,50.5,40.5,30.5,50.5,20.5,20.5,40.5,40.5,30.5,40.5,50.5}; //Vector containing all particle positions.

// calculation parameters
int Nequisteps = 2000; 	// the number of equilibration steps during which the stepsize may be adjusted.
double desiredaccept = 0.25;// target acceptance for monte carlo moves
int NMCsteps = 42000; 	// the number of MC steps
int n_moves = 1; 		// the number of particles moved with monte carlo per step
int underflowinterval=20;  // after how many segments should the propagator of the polymer chain be adjusted to prevent an underflow
int underflowprotect =1;	// should the calculation use underflow protection.
int maxiter = 1000;		// the number of maximum number iteration the program may tajke to solve a problem



//convergence parameters
double alphaweight =1;
int m_user = 8;		// number of steps kept in memeory by the diis method 8 seems optimal the higher the value tha faster the convergence is but also the less stable.
double eta = 0.1; 		// how quickly the alpha parameters should be adjusted. the larger the quicker it will change but it may start to oscillate. a safe value =is 0.1
double tolerance = 1e-7; 	// How accurately should the the result be calculated. The smaller the more accurate.  Typically between 1e-5 to 1e-7
int mmin = 6;			// when the calculation fails to find the correct result with _user it tries to solve it remebering fewer steps. 

// output
int vtkstore= 100; 		// How many monte carlo steps have to take place before a new density profile is stored.




//double H_poslist[Npart*3] = {10.5,10.5,10.5,10.5,10.5,0.5,10.5,0.5,10.5,0.5,10.5,10}; //Vector containing all particle positions.
//double H_oldposlist[Npart*3] ={10.5,10.5,10.5,10.5,10.5,0.5,10.5,0.5,10.5,0.5,10.5,10}; //Vector to store all old particle positions incase our Monte Carlo move is not accepted.
//double H_poslist[Npart*3] = {10.5,10.5,10.5,10.5,10.5,0.5,10.5,0.5,10.5,0.5,10.5,10.5,0.5,10.5,0.5,0.5,0.5,10.5,0.5,20.5,10.5,10.5,0.5,0.5,20.5,0.5,10.5,0.5,10.5,20}; //Vector containing all particle positions.
//double H_oldposlist[Npart*3] ={10.5,10.5,10.5,10.5,10.5,0.5,10.5,0.5,10.5,0.5,10.5,10.5,0.5,10.5,0.5,0.5,0.5,10.5,0.5,20.5,10.5,10.5,0.5,0.5,20.5,0.5,10.5,0.5,10.5,20}; //Vector to store all old particle positions incase our Monte Carlo move is not accepted.


//double H_poslist[Npart*3] = {10.5,10.5,10}; //Vector containing all particle positions.
//double H_oldposlist[Npart*3] = {10.5,10.5,10}; //Vector to store all old particle positions incase our Monte Carlo move is not accepted.