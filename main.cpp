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

#define  focal_len 2.5    //�������,��λmm
#define  cell_size 0.0034  //�����Ԫ��С,��λmm
#define  img_scale  4     // ԭʼ�ֱ���640->ʹ�õķֱ���160
#define  time_betwen_image 0.04  //ÿ֡����ʱ��
#define  mheight 1.0       //��������õ��ĸ߶�
#define sign(x) (( x > 0 ) - ( x < 0 ))
#define CHAR_BIT 8
#define BITMASK(b) (1 << ((b) % CHAR_BIT))  //��1����(b) % CHAR_BIT)λ
#define BITSLOT(b) ((b) / CHAR_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + CHAR_BIT - 1) / CHAR_BIT)
//��ͼ���н�imu
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
	float pixel_flow_x=0.0,pixel_flow_y=0.0;//ͼ�����������
	float sumx =0.0,sumy =0.0;              //ͼ������ۼӽ��
	int nFrmNum=0;                         //ͼ��֡��ͳ��
	int nFrmNum_tp=0;                      //����ȥ������ʱ����֡ͼ��Ŀ��Ʊ���
	float pitch=0.0,roll=0.0,last_pitch=0.0,last_roll=0.0,ax=0.0,ay=0.0;  //IMU����
	
    Mat pFrame_col;   //����RGBͼ��
	Mat pFrame_cur;	  //����grayͼ��

	float x,y;  //������ʾ�ľ��ο����Ͻ������
	CvRect rect; //������ʾ�ľ��ο�
	rect.x = 70;
	rect.y = 50;
	x = rect.x;
	y = rect.y;
	rect.width = 20;
	rect.height = 20;
	
	/********��������Ƶ�ļ�*******************/
    //char *imgPath = "E:\\test video\\flow_test\\indoor\\20151211\\result_2_160_120.avi";

	//char *imgPath = "E:\\test video\\flow_test\\indoor\\20160108_slow.avi";
    char *imgPath = "video_PV_320X240-3.avi";
	VideoCapture inputVideo(imgPath);  
	if ( !inputVideo.isOpened())
	{
		return 0;
	}
	inputVideo.set(CV_CAP_PROP_POS_FRAMES,10);
    //д��Ƶ�ļ�
	VideoWriter output;
	output.open("result20151217_0.avi",CV_FOURCC('M','J','P','G'),25,Size(160,120));
	//д�ĵ�
	FILE *myfile_w = fopen("write_opticalflow.txt","w");
   

	Flow_info_t *flow_struct_test=new Flow_info_t; //��ƥ������ṹ��
	Line_info_t *line_struct_test=new Line_info_t; //ֱ�߼������ṹ��

	kalman_s* myKalman = kalmanInit();  //kalman�˲���ʼ��;
	
	int flow_lineDetect_flag = 3;//ʹ�ÿ�ƥ�仹��ֱ�߼��ı�־
	int sad_method_staus;        //��ƥ�䷵�صı�־
	int lineDetect_method_staus; //ֱ�߼�ⷵ�صı�־

	while(true)
	{
		inputVideo >> pFrame_col;
		if(pFrame_col.empty()) break;
		if(pFrame_col.rows ==240)
			resize(pFrame_col,pFrame_col,Size(160,120));

		cvtColor(pFrame_col,pFrame_cur,CV_RGB2GRAY);
		/***************************����ʼ********************************************/

		nFrmNum_tp++;     //��ֹ�տ�������ͷ��ʱ��ͼ���ǰ��ģ������������ʱ��׼ȷ��������ȥ��ʼ��5֡
		if(nFrmNum_tp<5)
			continue;

	    //��һ������ͼ������������
		if(nFrmNum%100 == 0)
		{
	    //����ͼ���ݶȾ�ֵ
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

        //�ڶ�����֧һ��������ḻ����sad����
		if(flow_lineDetect_flag == 0)
		{
			sad_method_staus = flow_run(pFrame_cur.data,flow_struct_test);  //��ƥ�䷽�������
			if(sad_method_staus==0)
			{
				cout<<"sad match method unstable!!!"<<endl;  //������SADƥ���㷨���ȶ���

// 				if(lineDetect_method_staus == 0)
// 					cout<<"sad match and line detect method both unstable!!!!"<<endl;
			}
			pixel_flow_x = flow_struct_test->pixel_flow_x;
			pixel_flow_y = flow_struct_test->pixel_flow_y;
		}
		//�ڶ�����֧����������ϡ�����ֱ�߼�ⷽ��
		else if(flow_lineDetect_flag == 1)
		{
			lineDetect_method_staus = lineDetect_run(pFrame_col,line_struct_test); //ֱ�߼�ⷽ�������
			pixel_flow_x = line_struct_test->pixel_flow_x;
			pixel_flow_y = line_struct_test->pixel_flow_y;
			if(lineDetect_method_staus == 0)
			{
//				flow_lineDetect_flag = 0;
			    cout<<"line detetect method unstable!!!"<<endl;  //������ֱ�߼���㷨���ȶ���
			}
		}
		//�ڶ�����֧�������������㷨��ʧЧ��Ŀǰ��δ����÷�֧��
		else
		{
			pixel_flow_x = 0.0;
			pixel_flow_y = 0.0;
		}

		//������������IMU������ת����(Ŀǰ��δ����ƫ�����򲹳�)
		DecodeIMU(pFrame_cur,&pitch,&roll); //��ͼ���н��pitch��roll
		if(nFrmNum>0)
		{
			float offsetx = (-focal_len*(pitch-last_pitch)*0.1*3.14159)/(180*cell_size*img_scale);
			if(offsetx>8)
				offsetx = 8;
			else if(offsetx<-8)
				offsetx = -8;

			offsetx = 0.0;    //û��IMU����ʱ���еĻ���ȥ�������
			pixel_flow_x -=offsetx;

			float offsety = (-focal_len*(roll-last_roll)*0.1*3.14159)/(180*cell_size*img_scale);
			if(offsety>8)
				offsety = 8;
			else if(offsety<-8)
				offsety = -8;

			offsety = 0.0;    //û��IMU����ʱ���еĻ���ȥ�������
			pixel_flow_y -=offsety;
		}
		last_pitch=pitch;
		last_roll=roll;

		//�����ۼ�ƫ����
		sumx +=pixel_flow_x;
		sumy +=pixel_flow_y;
		fprintf(myfile_w,"%d  %f   %f \n",nFrmNum,sumx,sumy);

		//���Ĳ���ת��Ϊ������
		float vd_x = pixel_flow_x*mheight*img_scale*cell_size/(focal_len*time_betwen_image);
		float vd_y = -pixel_flow_y*mheight*img_scale*cell_size/(focal_len*time_betwen_image);  //ע�������и�����ϵת���Ĺ�ϵ��ʹ��ͼ�����������˻�����ϵһ��

		//���岽;��������������ٶȼ�kalman�ں�
		DecodeIMU1(pFrame_cur,&ay,&ax);   //��ͼ���н��pitch��roll
		ax = 0.0;ay = 0.0;    //���û��IMU����ʱ���еĻ���ȥ�������

		myKalman->Z[0][0] = vd_x;
		myKalman->Z[1][0] = vd_y;        //����ֵ����
		myKalman->Z[2][0] = 0.001*ax*time_betwen_image;  
		myKalman->Z[3][0] = 0.001*ay*time_betwen_image;//ax,ay�ĵ�λΪmm/s2�����Գ�0.001

		kalmanProcess(myKalman);  //kalman������

		float mvd_x = myKalman->X[0][0];           //������������Ž��
		float mvd_y = myKalman->X[1][0];

//		cout<<"mvd_x: "<<mvd_x<<" mvd_y: "<<mvd_y<<endl;

        /***************************�������********************************************/

		//������ʾ�ľ��ο�
		x += pixel_flow_x;
		y += pixel_flow_y;
		rect.x = x;
		rect.y = y;
		if(rect.x<0||rect.x>pFrame_col.cols||rect.y<0||rect.y>pFrame_col.rows)
		{
			x = (pFrame_col.cols-rect.width)/2;
			y = (pFrame_col.rows-rect.height)/2;
		}

		//ͼ������ʾ������ֵ
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
		//�����þ�����ʾ
		rectangle(pFrame_col,cvPoint(rect.x,rect.y),cvPoint(rect.x+rect.width,rect.y+rect.height),Scalar(255,255,255));

		//��ʾͼ��
 		imshow("1",pFrame_col);

  		char key = waitKey(2);
		if(key == 27)
			break;

		output<<pFrame_col;
		nFrmNum++;
	}
	fclose(myfile_w);

}


