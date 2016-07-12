#include "stdafx.h"
#include "flow_linux.h"
#include "lineProcess.h"
#include <opencv2/opencv.hpp>
#include"opencv2/core/core_c.h"
#include"opencv2/highgui/highgui_c.h"
#include"opencv2/imgproc/imgproc_c.h"
#include"opencv2/core/core.hpp"
#include"opencv2/highgui/highgui.hpp"
#include"opencv2/imgproc/imgproc.hpp"
#include "mykalman.h"
#include <iostream>
#include <strstream>

using namespace cv;
using namespace std;

#define  focal_len 2.5    //相机焦距,单位mm
#define  cell_size 0.0034  //相机靶元大小,单位mm
#define  img_scale  4     // 原始分辨率640->使用的分辨率160
#define  time_betwen_image 0.04  //每帧处理时间
#define  mheight 1.0       //传感器测得到的高度
#define sign(x) (( x > 0 ) - ( x < 0 ))
#define CHAR_BIT 8
#define BITMASK(b) (1 << ((b) % CHAR_BIT))  //把1左移(b) % CHAR_BIT)位
#define BITSLOT(b) ((b) / CHAR_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + CHAR_BIT - 1) / CHAR_BIT)
//从图像中解imu
void DecodeIMU(Mat img,float* pitch,float* roll)
{
	Mat imagecode(img,Rect(0,0,64,4)) ;
	for (int y = 0; y <4; y += 2)
	{
		uchar bitarray[BITNSLOTS(32)];
		int aaa=0;
		for (int x = 0; x<64; x += 2)
		{
			int side[4];
			side[0]=imagecode.at<uchar>((y)*64+x);
			side[1]=imagecode.at<uchar>((y+1)*64+x);
			side[2]=imagecode.at<uchar>((y)*64+x+1);
			side[3]=imagecode.at<uchar>((y+1)*64+x+1);
			int sum=side[0]+side[1]+side[2]+side[3];
			if(sum>400)
			{
				BITSET(bitarray,x/2);
			}
			else
			{
				BITCLEAR(bitarray,x/2);
			}
		}
		if(y<2)
		{
			*roll= *((float *)(bitarray));
		}
		else
		{
			*pitch = *((float *)(bitarray));
		}
	}
}
void DecodeIMU1(Mat img,float* pitch,float* roll)
{
	Mat imagecode(img,Rect(80,0,64,4)) ;
	for (int y = 0; y <4; y += 2)
	{
		uchar bitarray[BITNSLOTS(32)];
		int aaa=0;
		for (int x = 0; x<64; x += 2)
		{
			int side[4];
			side[0]=imagecode.at<uchar>((y)*64+x);
			side[1]=imagecode.at<uchar>((y+1)*64+x);
			side[2]=imagecode.at<uchar>((y)*64+x+1);
			side[3]=imagecode.at<uchar>((y+1)*64+x+1);
			int sum=side[0]+side[1]+side[2]+side[3];
			if(sum>400)
			{
				BITSET(bitarray,x/2);
			}
			else
			{
				BITCLEAR(bitarray,x/2);
			}
		}
		if(y<2)
		{
			*roll= *((float *)(bitarray));
		}
		else
		{
			*pitch = *((float *)(bitarray));
		}
	}
}

int main()
{
	float pixel_flow_x=0.0,pixel_flow_y=0.0;//图像光流计算结果
	float sumx =0.0,sumy =0.0;              //图像光流累加结果
	int nFrmNum=0;                         //图像帧数统计
	int nFrmNum_tp=0;                      //用于去除开机时若干帧图像的控制变量
	float pitch=0.0,roll=0.0,last_pitch=0.0,last_roll=0.0,ax=0.0,ay=0.0;  //IMU数据
	
    Mat pFrame_col;   //输入RGB图像
	Mat pFrame_cur;	  //输入gray图像

	float x,y;  //用于显示的矩形框左上角坐标点
	CvRect rect; //用于显示的矩形框
	rect.x = 70;
	rect.y = 50;
	x = rect.x;
	y = rect.y;
	rect.width = 20;
	rect.height = 20;
	
	/********打开输入视频文件*******************/
    //char *imgPath = "E:\\test video\\flow_test\\indoor\\20151211\\result_2_160_120.avi";

	//char *imgPath = "E:\\test video\\flow_test\\indoor\\20160108_slow.avi";
    char *imgPath = "video_PV_320X240-3.avi";
	VideoCapture inputVideo(imgPath);  
	if ( !inputVideo.isOpened())
	{
		return 0;
	}
	inputVideo.set(CV_CAP_PROP_POS_FRAMES,10);
    //写视频文件
	VideoWriter output;
	output.open("result20151217_0.avi",CV_FOURCC('M','J','P','G'),25,Size(160,120));
	//写文档
	FILE *myfile_w = fopen("write_opticalflow.txt","w");
   

	Flow_info_t *flow_struct_test=new Flow_info_t; //块匹配输出结构体
	Line_info_t *line_struct_test=new Line_info_t; //直线检测输出结构体

	kalman_s* myKalman = kalmanInit();  //kalman滤波初始化;
	
	int flow_lineDetect_flag = 3;//使用块匹配还是直线检测的标志
	int sad_method_staus;        //块匹配返回的标志
	int lineDetect_method_staus; //直线检测返回的标志

	while(true)
	{
		inputVideo >> pFrame_col;
		if(pFrame_col.empty()) break;
		if(pFrame_col.rows ==240)
			resize(pFrame_col,pFrame_col,Size(160,120));

		cvtColor(pFrame_col,pFrame_cur,CV_RGB2GRAY);
		/***************************处理开始********************************************/

		nFrmNum_tp++;     //防止刚开启摄像头的时候，图像是暗的，进行纹理分析时不准确，所以舍去开始的5帧
		if(nFrmNum_tp<5)
			continue;

	    //第一步：对图像进行纹理分析
		if(nFrmNum%100 == 0)
		{
	    //计算图像梯度均值
		Mat row_der,col_der,grad,angles;
		Sobel(pFrame_cur, row_der, CV_32F, 0, 1); 
		Sobel(pFrame_cur, col_der, CV_32F, 1, 0); 
		cartToPolar(col_der, row_der, grad, angles, true);
		Scalar mean,dev;
		meanStdDev(grad,mean,dev);
		float my_mean = mean[0];
		cout<<"image grad ave:"<<my_mean<<endl;
		if(my_mean>5)
			flow_lineDetect_flag = 0;
		else
            flow_lineDetect_flag = 1;
		}

        //第二步分支一：若纹理丰富进入sad方法
		if(flow_lineDetect_flag == 0)
		{
			sad_method_staus = flow_run(pFrame_cur.data,flow_struct_test);  //块匹配方法主入口
			if(sad_method_staus==0)
			{
				cout<<"sad match method unstable!!!"<<endl;  //报警：SAD匹配算法不稳定。

// 				if(lineDetect_method_staus == 0)
// 					cout<<"sad match and line detect method both unstable!!!!"<<endl;
			}
			pixel_flow_x = flow_struct_test->pixel_flow_x;
			pixel_flow_y = flow_struct_test->pixel_flow_y;
		}
		//第二步分支二：若纹理稀疏进入直线检测方法
		else if(flow_lineDetect_flag == 1)
		{
			lineDetect_method_staus = lineDetect_run(pFrame_col,line_struct_test); //直线检测方法主入口
			pixel_flow_x = line_struct_test->pixel_flow_x;
			pixel_flow_y = line_struct_test->pixel_flow_y;
			if(lineDetect_method_staus == 0)
			{
//				flow_lineDetect_flag = 0;
			    cout<<"line detetect method unstable!!!"<<endl;  //报警：直线检测算法不稳定。
			}
		}
		//第二步分支三：上述两种算法都失效，目前暂未进入该分支。
		else
		{
			pixel_flow_x = 0.0;
			pixel_flow_y = 0.0;
		}

		//第三步：利用IMU进行旋转补偿(目前暂未进行偏航方向补偿)
		DecodeIMU(pFrame_cur,&pitch,&roll); //从图像中解出pitch，roll
		if(nFrmNum>0)
		{
			float offsetx = (-focal_len*(pitch-last_pitch)*0.1*3.14159)/(180*cell_size*img_scale);
			if(offsetx>8)
				offsetx = 8;
			else if(offsetx<-8)
				offsetx = -8;

			offsetx = 0.0;    //没有IMU数据时，有的话就去掉该语句
			pixel_flow_x -=offsetx;

			float offsety = (-focal_len*(roll-last_roll)*0.1*3.14159)/(180*cell_size*img_scale);
			if(offsety>8)
				offsety = 8;
			else if(offsety<-8)
				offsety = -8;

			offsety = 0.0;    //没有IMU数据时，有的话就去掉该语句
			pixel_flow_y -=offsety;
		}
		last_pitch=pitch;
		last_roll=roll;

		//光流累加偏移量
		sumx +=pixel_flow_x;
		sumy +=pixel_flow_y;
		fprintf(myfile_w,"%d  %f   %f \n",nFrmNum,sumx,sumy);

		//第四步：转换为物理量
		float vd_x = pixel_flow_x*mheight*img_scale*cell_size/(focal_len*time_betwen_image);
		float vd_y = -pixel_flow_y*mheight*img_scale*cell_size/(focal_len*time_betwen_image);  //注意这里有个坐标系转换的关系，使得图像的输出与无人机坐标系一致

		//第五步;光流物理量与加速度计kalman融合
		DecodeIMU1(pFrame_cur,&ay,&ax);   //从图像中解出pitch，roll
		ax = 0.0;ay = 0.0;    //如果没有IMU数据时，有的话就去掉该语句

		myKalman->Z[0][0] = vd_x;
		myKalman->Z[1][0] = vd_y;        //测量值输入
		myKalman->Z[2][0] = 0.001*ax*time_betwen_image;  
		myKalman->Z[3][0] = 0.001*ay*time_betwen_image;//ax,ay的单位为mm/s2，所以乘0.001

		kalmanProcess(myKalman);  //kalman处理函数

		float mvd_x = myKalman->X[0][0];           //最终输出的最优结果
		float mvd_y = myKalman->X[1][0];

//		cout<<"mvd_x: "<<mvd_x<<" mvd_y: "<<mvd_y<<endl;

        /***************************处理结束********************************************/

		//用于显示的矩形框
		x += pixel_flow_x;
		y += pixel_flow_y;
		rect.x = x;
		rect.y = y;
		if(rect.x<0||rect.x>pFrame_col.cols||rect.y<0||rect.y>pFrame_col.rows)
		{
			x = (pFrame_col.cols-rect.width)/2;
			y = (pFrame_col.rows-rect.height)/2;
		}

		//图像上显示光流数值
		char showflow[256];
		sprintf(showflow,"x_flow is:%0.4f",pixel_flow_x);
		putText(pFrame_col,showflow,Point(10,30),1,1,Scalar(0,255,0));
		sprintf(showflow,"y_flow is:%0.4f",pixel_flow_y);
		putText(pFrame_col,showflow,Point(10,50),1,1,Scalar(0,255,0));

		if(flow_lineDetect_flag == 0)
		{
			sprintf(showflow,"use optical flow !!!");
			putText(pFrame_col,showflow,Point(10,80),1,1,Scalar(0,0,255));
		}
		else if(flow_lineDetect_flag == 1)
		{
			sprintf(showflow,"use line detection!!!");
			putText(pFrame_col,showflow,Point(10,80),1,1,Scalar(0,0,255));
		}
		else
		{
			sprintf(showflow,"not any method used!!!");
			putText(pFrame_col,showflow,Point(10,80),1,1,Scalar(0,0,255));
		}
		//测试用矩形显示
		rectangle(pFrame_col,cvPoint(rect.x,rect.y),cvPoint(rect.x+rect.width,rect.y+rect.height),Scalar(255,255,255));

		//显示图像
 		imshow("1",pFrame_col);

  		char key = waitKey(2);
		if(key == 27)
			break;

		output<<pFrame_col;
		nFrmNum++;
	}
	fclose(myfile_w);

}


