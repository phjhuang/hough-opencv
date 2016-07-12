#ifndef MY_KALMAN_
#define MY_KALMAN_

typedef struct kalman_s_t
{   
	//The number of dimensions in signal (ns);
	//The number of dimensions in measurement (nm)
	//The relation between one signal and the next (ns*ns)
	double** A;//4*4

	//The relation between measurement and signal (nm*ns)
	double** H;//4*4

	//the state vector-actually a column matrix (ns*1)
	double** X;//4*1

	//The prediction vector-actually a column matrix (ns*1)
	//matrix x_pred;

	//The measurement vector-again a column matrix (nm*1)
	double** Z;//4*1

	//The error covariance matrix (ns*ns)
	double** P;//4*4

	//The Kalman gain matrix (ns*nm)
	double** K;//4*4

	//The covariance matrix of system (ns*ns)
	double** Q;//4*4

	//The covariance matrix of measurement (nm*nm)
	double** R;//4*4

	//temp var
	double **temp_1;    //（4*1）
	double **temp_2;    //(4*4)
	double **temp_2_1;   //(4*4)
	double **temp_3_1;  //(4*1)
	double **temp_3_2;  //(4*1)
	double **temp_4_1;  //(4*4)
	double **temp_4_2;  //(4*4)
	double **temp_4_3;  //（4*4）
	double **temp_4_4;  //（4*4）

	double **E;         //（4*4）
	double **temp_5;   //(4*4)
	double **temp_5_1;  //（4*4）
	double **temp_5_2;  //(4*4)
}kalman_s;

kalman_s*  kalmanInit();
void kalmanProcess(kalman_s *mykalman);

void kalmanFree(kalman_s *mykalman);

void freeMyMatrix(double **matrix,int rows,int column);

//矩阵的一些重载操作
bool multiplyMatrix(double** matrix1, int row1, int column1, double** matrix2, int row2, int column2, double** matrix3);

bool transpositionMatrix(double** matrix, double** transpositionMatirx, int row, int column);

bool addMatrix(double** matrix1, int row1, int column1, double** matrix2, int row2, int column2, double** matrix3);

bool subMatrix(double** matrix1, int row1, int column1, double** matrix2, int row2, int column2, double** matrix3);

double DeterminantOfMatrix(double** matrix, int row, int column);

bool InverseMatrix(double** matrix, double** inverseMatirx, int row, int column);

void printMatrix(double** matrix, int row, int column);

#endif  //MY_KALMAN_
