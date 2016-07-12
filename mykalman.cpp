 #include "stdafx.h"
#include "mykalman.h"

#define NS  4         //信号的维度，即状态个数
#define NM  4         //测量值的维数

kalman_s* kalmanInit()
{
	int i,j;
	kalman_s *myKalman = (kalman_s *)malloc(sizeof(kalman_s));

	//*********************kalman滤波矩阵分配内存*************/
	myKalman->A = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
	{
		myKalman->A[i] = (double*)malloc(sizeof(double)*NS);
		memset(myKalman->A[i],0,sizeof(double)*NS);
	}

	myKalman->H = (double**)malloc(sizeof(double*)*NM);         
	for(i = 0; i < NM; i++)
	{
		myKalman->H[i] = (double*)malloc(sizeof(double)*NS);
		memset(myKalman->H[i],0,sizeof(double)*NS);
	}

	myKalman->K = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
	{
		myKalman->K[i] = (double*)malloc(sizeof(double)*NM);
		memset(myKalman->K[i],0,sizeof(double)*NM);
	}

	myKalman->Z = (double**)malloc(sizeof(double*)*NM);         
	for(i = 0; i < NM; i++)
	{
		myKalman->Z[i] = (double*)malloc(sizeof(double)*1);
		memset(myKalman->Z[i],0,sizeof(double)*1);
	}

	myKalman->P = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->P[i] = (double*)malloc(sizeof(double)*NS);

	myKalman->Q = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->Q[i] = (double*)malloc(sizeof(double)*NS);

	myKalman->R = (double**)malloc(sizeof(double*)*NM);         
	for(i = 0; i < NM; i++)
		myKalman->R[i] = (double*)malloc(sizeof(double)*NM);

	myKalman->X = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
	{
		myKalman->X[i] = (double*)malloc(sizeof(double)*1);
	}
	//一些矩阵计算用到的临时变量
	myKalman->temp_1 = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->temp_1[i] = (double*)malloc(sizeof(double)*1);

	myKalman->temp_2 = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->temp_2[i] = (double*)malloc(sizeof(double)*NS);

	myKalman->temp_2_1 = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->temp_2_1[i] = (double*)malloc(sizeof(double)*NS);

	myKalman->temp_3_1 = (double**)malloc(sizeof(double*)*NM);         
	for(i = 0; i < NM; i++)
		myKalman->temp_3_1[i] = (double*)malloc(sizeof(double)*1);

	myKalman->temp_3_2 = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->temp_3_2[i] = (double*)malloc(sizeof(double)*1);

	myKalman->temp_4_1 = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->temp_4_1[i] = (double*)malloc(sizeof(double)*NM);

	myKalman->temp_4_2= (double**)malloc(sizeof(double*)*NM);         
	for(i = 0; i < NM; i++)
		myKalman->temp_4_2[i] = (double*)malloc(sizeof(double)*NM);

	myKalman->temp_4_3 = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->temp_4_3[i] = (double*)malloc(sizeof(double)*NM);

	myKalman->temp_4_4 = (double**)malloc(sizeof(double*)*NM);         
	for(i = 0; i < NM; i++)
		myKalman->temp_4_4[i] = (double*)malloc(sizeof(double)*NM);

	myKalman->temp_5= (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->temp_5[i] = (double*)malloc(sizeof(double)*NS);

	myKalman->temp_5_1 = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->temp_5_1[i] = (double*)malloc(sizeof(double)*NS);

	myKalman->temp_5_2 = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->temp_5_2[i] = (double*)malloc(sizeof(double)*NS);

	myKalman->E = (double**)malloc(sizeof(double*)*NS);         
	for(i = 0; i < NS; i++)
		myKalman->E[i] = (double*)malloc(sizeof(double)*NS);

	//***一些矩阵的初始化，注意这里默认的是NS=4,NM=4的初始化,若修改，则下面的赋值维度也需要修改*****/
	myKalman->A[0][0] = 1.0;
	myKalman->A[1][1] = 1.0;
	myKalman->A[2][2] = 1.0;
	myKalman->A[3][3] = 1.0;       //状态转移矩阵的初始化
	myKalman->A[0][2] = 1.0;
	myKalman->A[1][3] = 1.0;

	for (i=0;i<NS;i++)
		for (j=0;j<NS;j++)
		{		
			if(i==j)
			{
				myKalman->P[i][j]=1.0;
				myKalman->Q[i][j]=0.001;  //过程噪声，值越小响应越慢
				myKalman->E[i][j]=1.0;
			}
			else
			{
				myKalman->P[i][j]=0;
				myKalman->Q[i][j]=0;
				myKalman->E[i][j]=0;
			}
		}

	for (i=0;i<NM;i++)
		for (j=0;j<NS;j++)
		{
			if(i==j)
				myKalman->H[i][j]=1.0;
			else
				myKalman->H[i][j]=0.0;
		}

	for (i=0;i<NM;i++)
		for (j=0;j<NM;j++)
		{
			if(i==j)
              myKalman->R[i][j]=0.1;
			else
			  myKalman->R[i][j]=0.0;
		}

	for(int k=0;k<NS;k++)
	myKalman->X[k][0]=0.0;       //状态变量的初始化

	return myKalman;
}
void kalmanProcess(kalman_s *mykalman)
{
	int i,j;
	for(i=0;i<NM;i++)
		memset(mykalman->temp_4_4[i],0,sizeof(double)*NM);

	//第一个公式：x(k|k-1) = Ax(k-1|k-1)
	multiplyMatrix(mykalman->A,NS,NS,mykalman->X,NS,1,mykalman->temp_1);                //X=A*X
	for (i=0;i<NS;i++)
		mykalman->X[i][0]=mykalman->temp_1[i][0];

	//第二个公式： P = A*P*A'+Q
	multiplyMatrix(mykalman->A,NS,NS,mykalman->P,NS,NS,mykalman->temp_2_1);	             //temp_2_1 = A*P
	transpositionMatrix(mykalman->A, mykalman->temp_2, NS, NS);                        //temp_2 = A'
	multiplyMatrix(mykalman->temp_2_1,NS,NS,mykalman->temp_2,NS,NS,mykalman->P);         //P = A*P*A’ 
	addMatrix(mykalman->P,NS,NS,mykalman->Q,NS,NS,mykalman->P);                         //P = A*P*A’+Q

	//第三个公式： X = X+K*[Z-H*X]
    multiplyMatrix(mykalman->H,NM,NS,mykalman->X,NS,1,mykalman->temp_3_1);               //temp_3_1=H*X
	subMatrix(mykalman->Z,NM,1,mykalman->temp_3_1,NM,1,mykalman->temp_3_1);             //temp_3_1=Z-H*X    
	multiplyMatrix(mykalman->K,NS,NM,mykalman->temp_3_1,NM,1,mykalman->temp_3_2);       //temp_3_2 = K*(Z-H*X)
	addMatrix(mykalman->X,NS,1,mykalman->temp_3_2,NS,1,mykalman->X);                 //X = X+ K*(Z-H*X)

	//第四个公式：K = P*H'/[H*P*H'+R]
	transpositionMatrix(mykalman->H,mykalman->temp_4_3,NM,NS);                      //temp_4_3 = H'
	multiplyMatrix(mykalman->P,NS,NS,mykalman->temp_4_3,NS,NM,mykalman->temp_4_1);    //temp_4_1 = P*H'
	multiplyMatrix(mykalman->H,NM,NS,mykalman->temp_4_1,NS,NM,mykalman->temp_4_2);    //temp_4_2 =H*P*H'
	addMatrix(mykalman->temp_4_2,NM,NM,mykalman->R,NM,NM,mykalman->temp_4_2);         //temp_4_2 =H*P*H'+R
	InverseMatrix(mykalman->temp_4_2, mykalman->temp_4_4, NM,NM);                  //temp_4_4=~(H*P*H'+R)
	multiplyMatrix(mykalman->temp_4_1,NS,NM,mykalman->temp_4_4,NM,NM,mykalman->K);   //K = P*H'*~(H*P*H'+R)

	//第五个公式：P = [I-K*H]*P
	multiplyMatrix(mykalman->K,NS,NM,mykalman->H,NM,NS,mykalman->temp_5);            //temp_5 = K*H
	subMatrix(mykalman->E,NS,NS,mykalman->temp_5,NS,NS,mykalman->temp_5_1);          //temp_5_1 = E-K*H
	multiplyMatrix(mykalman->temp_5_1,NS,NS,mykalman->P,NS,NS,mykalman->temp_5_2);  //temp_5_2 = (E-K*H)*P

	for (i=0;i<NS;i++)
		for (j=0;j<NS;j++)
          mykalman->P[i][j] = mykalman->temp_5_2[i][j];
}

void kalmanFree(kalman_s *mykalman)
{
	freeMyMatrix(mykalman->A,NS,NS);
	freeMyMatrix(mykalman->K,NS,NM);
	freeMyMatrix(mykalman->P,NS,NS);
	freeMyMatrix(mykalman->Q,NS,NS);
	freeMyMatrix(mykalman->R,NM,NM);
	freeMyMatrix(mykalman->X,NS,1);
	freeMyMatrix(mykalman->Z,NM,1);
	freeMyMatrix(mykalman->temp_1,NS,1);
	freeMyMatrix(mykalman->temp_2,NS,NS);
	freeMyMatrix(mykalman->temp_2_1,NS,NS);
	freeMyMatrix(mykalman->temp_3_1,NM,1);
	freeMyMatrix(mykalman->temp_3_2,NS,1);
	freeMyMatrix(mykalman->temp_4_1,NS,NM);
	freeMyMatrix(mykalman->temp_4_2,NM,NM);
	freeMyMatrix(mykalman->temp_4_3,NS,NM);
	freeMyMatrix(mykalman->temp_4_4,NM,NM);
	freeMyMatrix(mykalman->temp_5,NS,NS);
	freeMyMatrix(mykalman->temp_5_1,NS,NS);
	freeMyMatrix(mykalman->temp_5_2,NS,NS);
	freeMyMatrix(mykalman->E,NS,NS);

	free(mykalman);
}

void freeMyMatrix(double **matrix,int rows,int cloumn)
{
	int i;
	for (i=0;i<rows;i++)
	{
		free(matrix[i]);
		matrix[i]=NULL;
	}
	free(matrix);
	matrix = NULL;
}


/// 矩阵相乘
 bool multiplyMatrix(double** matrix1, int row1, int column1, double** matrix2, int row2, int column2, double** matrix3)
 {
	int i,j,k;
	if(column1 != row2)
	{
		printf("matrix cannot be multiply !!!\n");
		return false;
	}

	for(i=0;i<row1;i++)
		memset(matrix3[i],0,sizeof(double)*column2);

	for(i = 0; i < row1; i++)
	{
		for(j = 0; j < column2; j++)
		{
			for(k = 0; k < column1; k++)
			{
				matrix3[i][j] += matrix1[i][k]*matrix2[k][j];
			}
		}
	}
	return true;
}

 /// 求矩阵转置矩阵
 bool transpositionMatrix(double** matrix, double** transpositionMatirx, int row, int column)
 {
	 int i;
	 int j;

	 if(matrix == NULL || transpositionMatirx == NULL || row <= 0 || column <= 0)
	 {
		 printf("TranspositionMatrix error !!!\n");
		 return false;
	 }

	 for(i = 0; i < row; i++)
	 {
		 for(j = 0; j < column; j++)
		 {
			 transpositionMatirx[j][i] = matrix[i][j];
// 			 transpositionMatirx[j][i] = matrix[i][j];
// 			 transpositionMatirx[j][i] = matrix[i][j];
		 }
	 }
	 return true;
 }

 //矩阵相加
 bool addMatrix(double** matrix1, int row1, int column1, double** matrix2, int row2, int column2, double** matrix3)
 {
	 int i;
	 int j;
	 if (row1 != row2 || column1 !=column2 )
	 {
		 printf("add matrix error !!!\n");
		 return false;
	 }

	 for(i=0;i<row1;i++)
	 {
		 for (j=0;j<column1;j++)
		 {
			 matrix3[i][j] = matrix1[i][j]+matrix2[i][j];
		 }
	 }

	 return true;
 }

 //矩阵相减
 bool subMatrix(double** matrix1, int row1, int column1, double** matrix2, int row2, int column2, double** matrix3)
 {
	 int i;
	 int j;
	 if (row1 != row2 || column1 !=column2 )
	 {
		 printf("add matrix error !!!\n");
		 return false;
	 }

	 for(i=0;i<row1;i++)
	 {
		 for (j=0;j<column1;j++)
		 {
			 matrix3[i][j] = matrix1[i][j]-matrix2[i][j];
		 }
	 }
	 return true;
 }

 /// 求矩阵的行列式
 double DeterminantOfMatrix(double** matrix, int row, int column)
 {
	 int i,j,k,p,r; 
	 double X, temp = 1, temp1 = 1, s = 0, s1 = 0;
	 if(column == 2)
	 {
		 for(i=0;i < row;i++)
			 for(j=0;j <column;j++)
				 if((i+j)%2) 
					 temp1 *= matrix[i][j];
				 else 
					 temp *= matrix[i][j];
		 X = temp - temp1;
	 }
	 else
	 {
		 for(k = 0; k < column; k++)
		 {

			 for(i=0,j=k; i<row && j<column; i++,j++)
				 temp *= matrix[i][j];
			 if(row - i)
			 {
				 for(p = row-i, r = row-1; p > 0; p--,r--)
					 temp *= matrix[r][p-1];
			 }
			 s += temp;
			 temp = 1;
		 }
		 for(k = column-1; k >= 0; k--)
		 {
			 for(i=0, j=k; i<row && j>=0; i++,j--)
				 temp1 *= matrix[i][j];
			 if(row - i)
			 {
				 for(p = row-1, r=i; r <row; p--,r++)
					 temp1 *= matrix[r][p];
			 }
			 s1 += temp1;
			 temp1 = 1;
		 }
		 X = s-s1;
	 }
	 return X;
 }

 /// 求n阶矩阵的逆矩阵
 bool InverseMatrix(double** matrix, double** inverseMatirx, int row, int column)
 {
	 int i, j, x, y, k, l;
	 double **SP = NULL, **AB = NULL, **TempMatrix = NULL, X;
	 if(matrix == NULL || inverseMatirx == NULL || row <= 0 || column != row)
	 {
		 printf("matrix cannot be inversed !!!!\n");
		 return false;
	 }

	 SP = (double **)malloc(row*sizeof(double*));
	 AB = (double **)malloc(row*sizeof(double*));
	 TempMatrix = (double **)malloc(row*sizeof(double*));
	 X = DeterminantOfMatrix(matrix, row, column);
	 if(X == 0)
	 {
		 printf("det==0,matrix cannot be inversed !!!!\n");
		 return false;
	 }
	 X = 1/X;

	 for(i = 0; i < row; i++) 
	 {
		 SP[i] = (double *)malloc(column*sizeof(double));
		 AB[i] = (double *)malloc(column*sizeof(double));
		 TempMatrix[i] = (double *)malloc(column*sizeof(double));
	 }

	 //求矩阵伴随矩阵
	 for(i = 0; i < row; i++)
	 {
		 for(j = 0; j < column; j++)
		 {
			 for(k = 0; k < row; k++)
				 for(l = 0; l <column; l++)
					 TempMatrix[k][l] = matrix[k][l];	//TempMatrix变成matrix;
			 {
				 for(x = 0; x < column; x++)
					 TempMatrix[(i*column+x)/row][(i*column+x)%row] = 0;
				 for(y = 0; y < row; y++)
					 TempMatrix[y][j]=0;
				 TempMatrix[(int)((i*column+j)/row)][(i*column+j)%row] = 1;
				 SP[(int)((i*column+j)/row)][(i*column+j)%row] = DeterminantOfMatrix(TempMatrix, row, column);
				 AB[(int)((i*column+j)/row)][(i*column+j)%row] = X * SP[(int)((i*column+j)/row)][(i*column+j)%row];
			 }
		 }
	 }
	 transpositionMatrix(AB, inverseMatirx, row, column);

	 freeMyMatrix(AB,row,column);
	 freeMyMatrix(SP,row,column);
	 freeMyMatrix(TempMatrix,row,column);

	 return true;
 }

 void printMatrix(double** matrix, int row, int column)
 {
	 int m, n;
	 for(m = 0; m < row; m++)
	 {
		 for(n = 0; n < column; n++)
		 {
			 printf("%lf ",matrix[m][n]);
		 }
		 printf("\n");
	 }
 }