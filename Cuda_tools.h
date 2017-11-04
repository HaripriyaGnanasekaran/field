#ifndef CUDA_TOOLSxH
#define CUDA_TOOLSxH
  
void GPU_present(void);
double* AllocateMemoryOnDevice(int);
void FreeMemoryOnDevice(double*);
void TransferDataToHost(int, double*, double* );
void TransferDataToDevice(int, double*, double*);
void InitializeForward(int,int,int,double*,double*);
void InitializeForward(int,int,int,int,double*,double*);
void InitializeBackward(int ,int , double* ,double* );
void InitializeBackward(int ,int , int, double* ,double* );
void Forward(int,int,int, double*,double*,double,int,int); //1D cubic flat
void Forward(int,int,int,double*,double*,double, int*); //first order
void Forward(int,int,int,double*,double*,double,double, int*); //second order
void Backward(int,int, double*, double*,double*,double,int,int);
void Backward(int,int,double*,double*,double*,double, int*);//first order
void Backward(int,int,double*,double*,double*,double,double, int*); //second order
void Propagate(int, int, int, int, double* , double*, double, int*); //homefully general 
void Composition(int, int, int,double*, double*, double* );
void Composition(int, int, int, int, double*, double*, double* );
void CorrectDoubleCounting(int, int, double*, double*);
void CorrectDoubleCounting(int, int, int, double*, double*);
void NormPhi(int,int,double*,double);
void Zero(int,double*);
void Add(int, int, int, double* , double*);
double Norm2(int, int, double*);
#endif 
