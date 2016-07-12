#include "stdafx.h"
#include <opencv2/opencv.hpp> 
#include "lineProcess.h"

using namespace std;
using namespace cv;

#define myDebug

void getMyClusterLine(vector<Vec4i> lines,vector<vector<Vec4i>> &out_line_cluster)
{
	vector<Vec4i> templines = lines;
	vector<vector<Vec4i>> final_line_cluster(0);

	while (templines.size()>1)
	{
		double line_len_temp = 0.0;
		int maxIndex = 0;
		vector<double>line_length_vec;
		vector<double> line_theta_vec;       //直线斜率 
		vector<Point2f> centre_point_vec;
		vector<Vec4i> temp_clustered;
		vector<Vec4i> temp_un_clustered;

		//求剩下直线中的最长直线
		for (int i=0;i<templines.size();i++)
		{
			Point pt1(templines[i][0],templines[i][1]);
			Point pt2(templines[i][2],templines[i][3]);

			double line_length = sqrt((pt1.x-pt2.x)*(pt1.x-pt2.x)+(pt1.y-pt2.y)*(pt1.y-pt2.y));
			line_length_vec.push_back(line_length);

			//保存线段的中心点，用于后续求点到直线的距离
			Point2f t = Point2f((pt2.x+pt1.x)/2.0,(pt2.y+pt1.y)/2.0);
			centre_point_vec.push_back(t);

			//求角度，弧度表示
			double k=asin((double)(pt2.y - pt1.y)/line_length) ;
			line_theta_vec.push_back(k);

			if(line_length>line_len_temp)
			{
				line_len_temp = line_length;
				maxIndex = i;
			}
		}

		//求最长的那条直线方程
		Point pt_m1(templines[maxIndex][0],templines[maxIndex][1]);
		Point pt_m2(templines[maxIndex][2],templines[maxIndex][3]);

		Vec4i keepLongestLine;
		//提取与最长直线平行且距离较近的点
		for (int i=0;i<line_theta_vec.size();i++)
		{
			double dis = 1000.0;

			if(i==maxIndex)
			{
				keepLongestLine = templines[i];
			}
			else
			{
				//求两直线的夹角		
				double intersection = intersection_theta(templines[i],templines[maxIndex]);
				//角度很接近的线段
				if(intersection<0.3) 
				{

					if(fabs(line_theta_vec[maxIndex]-1.5708)<0.01)
					{
						dis = fabs(double(templines[maxIndex][0]-templines[i][0])); 
					}
					else
					{
						dis = fabs((pt_m1.x-pt_m2.x)*(centre_point_vec[i].y-pt_m2.y)-(pt_m1.y-pt_m2.y)*(centre_point_vec[i].x-pt_m2.x))
							/sqrt((pt_m1.x-pt_m2.x)*(pt_m1.x-pt_m2.x)+(pt_m1.y-pt_m2.y)*(pt_m1.y-pt_m2.y));  //点到直线的距离
					}

					//距离很接近的直线
					if(dis<5.0)      //160*120
						temp_clustered.push_back(templines[i]);
					else
						temp_un_clustered.push_back(templines[i]);       
				}
				else
				{
					temp_un_clustered.push_back(templines[i]);
				}
			}
		}

		temp_clustered.push_back(keepLongestLine);

		for (int i=0;i<temp_un_clustered.size();i++)
			templines[i] = temp_un_clustered[i];

		templines.resize(temp_un_clustered.size());

		final_line_cluster.push_back(temp_clustered);
	}

	out_line_cluster = final_line_cluster;
}

void mycvGetPoint(Vec4i line1, Vec4i line2,Point2f &out_inter_point)   //求两直线交点
{
	float x1 = line1[0];
	float x2 = line1[2];
	float y1 = line1[1];
	float y2 = line1[3];

	float x3 = line2[0];
	float x4 = line2[2];
	float y3 = line2[1];
	float y4 = line2[3];
  
	if(x2!=x1&&x3!=x4)
	{
		double k1 = (y2-y1)/(x2-x1);
		double k2 = (y3-y4)/(x3-x4);

		if(abs(k1-k2)<0.0001)
		{
			out_inter_point.x = 10000.0;
			out_inter_point.y = 10000.0;
		}
		else
		{
			out_inter_point.x = (k1*x1-k2*x3+y3-y1)/(k1-k2);
			out_inter_point.y = y1+(out_inter_point.x-x1)*k1;
		}
	}
	else if(x1 == x2 && x3 != x4)
	{
		out_inter_point.x  = x1;
		double k2 = (y3-y4)/(x3-x4);
		out_inter_point.y = k2*(out_inter_point.x-x3)+y3;
	}
	else if(x1!=x2&&x3==x4)
	{
		out_inter_point.x  = x3;
		double k1 = (y1-y2)/(x1-x2);
		out_inter_point.y = k1*(out_inter_point.x-x1)+y1;
	}
	else
	{
		out_inter_point.x = 100000.0;
		out_inter_point.y = 100000.0;
	}
}

//求两直线的夹角(锐角)
float intersection_theta(Vec4i line1, Vec4i line2)
{
	float vec1_x = line1[2] - line1[0];
	float vec1_y = line1[3] - line1[1];

	float vec2_x = line2[2] - line2[0];
	float vec2_y = line2[3] - line2[1];

	double s = (vec1_x*vec2_x+vec1_y*vec2_y)/sqrt((vec1_x*vec1_x+vec1_y*vec1_y)*(vec2_x*vec2_x+vec2_y*vec2_y));

	double theta = acos(fabs(s));

	return theta;
}

void filterLine(Mat img,vector<vector<Vec4i>> in_line_clustered,vector<intersectionPoint_s> &my_intersectionPoint_vec_tmp)
{
#ifdef myDebug
	Mat showImage;
	img.copyTo(showImage);

	//*******聚类后的直线测试**************************
		for (int i=0;i<in_line_clustered.size();i++)
		{
			Scalar s; 
			if(i==0)
				s = Scalar(0,255,0);
			else if(i==1)
				s = Scalar(255,0,0);
			else if(i==2)
				s = Scalar(0,0,255);
			else
				s= Scalar(0,0,0);

			vector<Vec4i> temp = in_line_clustered[i];
			for (int j=0;j<temp.size();j++)
			{
				Point p1,p2;
				p1.x = temp[j][0];  //某类中最后面的那条直线存储最长直线
				p1.y = temp[j][1]; 
				p2.x = temp[j][2];
				p2.y = temp[j][3];
				line(showImage,p1,p2,s,1,CV_AA);
			}
		}
#endif // DEBUG

	    vector<Vec4i> intersection_point_use(0);
		int numof_clustered = in_line_clustered.size();   //聚类数
		for (int i=0;i<numof_clustered;i++)
		{
#ifdef myDebug
			Scalar s; 
			if(i==0)
				s = Scalar(0,255,0);
			else if(i==1)
				s = Scalar(255,0,0);
			else if(i==2)
				s = Scalar(0,0,255);
			else
				s= Scalar(0,0,0);
#endif
			vector<Vec4i> temp = in_line_clustered[i];
			int numof_lines =temp.size();     //每类的直线数
			double length_all_lines = 0.0;

			for (int j=0;j<numof_lines;j++)
			{
				Point p1,p2;
				p1.x = temp[j][0];  
				p1.y = temp[j][1]; 
				p2.x = temp[j][2];
				p2.y = temp[j][3];

				double dis = sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
				length_all_lines +=dis;
				//				 cout<<dis<<endl;
			}

			Point p1,p2,p3,p4;
			if(length_all_lines>50.0)          //160*120
			{
				intersection_point_use.push_back(temp[numof_lines-1]);
#ifdef myDebug
				p1.x = temp[numof_lines-1][0];  //某类中最后面的那条直线存储最长直线
				p1.y = temp[numof_lines-1][1]; 
				p2.x = temp[numof_lines-1][2];
				p2.y = temp[numof_lines-1][3];

				line(img,p1,p2,s,3,CV_AA);

				double leng = sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
				double k=asin((double)(p2.y - p1.y)/leng) ;

				double a = sin(k);double b = cos(k);

				if(k<=0)
				{
					p3.x = p1.x + 1000*b;
					p3.y = p1.y + 1000*a;
					p4.x = p1.x-1000*b;
					p4.y = p1.y - 1000*a;
				}
				else
				{
					p3.x = p1.x-1000*b;
					p3.y = p1.y - 1000*a;
					p4.x = p1.x+1000*b;
					p4.y = p1.y + 1000*a;
				} 

				line(img,p3,p4,s,1,CV_AA);	
#endif // myDebug
			}
		}

	//测试垂直直线的交点
	if (intersection_point_use.size()>=2)
	{
		Point2f intersection = Point(-1.0,-1.0);
		double intertheta=0.0;
		intersectionPoint_s my_intersectionPoint_s;

		for (int i=0;i<intersection_point_use.size()-1;i++)
		{
			for (int j=i+1;j<intersection_point_use.size();j++)
			{
				mycvGetPoint(intersection_point_use[i], intersection_point_use[j],intersection);
				intertheta = intersection_theta(intersection_point_use[i], intersection_point_use[j]);

				if(fabs(intertheta-1.5708)<0.4&&intersection.x>0&&intersection.x<img.cols&&intersection.y>0&&intersection.y<img.rows)  //有时候点在图像外其实也是可以的，
				{
					my_intersectionPoint_s.intersectionPoint = intersection;
					my_intersectionPoint_s.Line1 = intersection_point_use[i];
					my_intersectionPoint_s.Line2 = intersection_point_use[j];
					my_intersectionPoint_s.intersectionPoint2LineDis = twoLineDis(intersection_point_use[i],intersection_point_use[j],intersection);
					my_intersectionPoint_vec_tmp.push_back(my_intersectionPoint_s);
#ifdef myDebug					
					circle(img,intersection,3,Scalar(0,255,0),3);
#endif // myDebug
				}
			}
		}
	}
#ifdef myDebug
 	  	imshow("show_image",showImage);	
 	  	waitKey(2);
#endif // myDebug
}

//当前帧交点的滤波，主要是去除相隔很近的点
void intersectionPointFilter(vector<intersectionPoint_s> in_my_intersectionPoint_vec,vector<intersectionPoint_s>&out_my_intersectionPoint_vec)
{
	int size = in_my_intersectionPoint_vec.size();
	vector<int>staus(size,1);

	if(size>1)
	for (int i=0;i<size-1;i++)
	{
		for (int j=i+1;j<size;j++)
		{
			Point2f p1 = in_my_intersectionPoint_vec[i].intersectionPoint;
			Point2f p2 = in_my_intersectionPoint_vec[j].intersectionPoint;
			float d1 = in_my_intersectionPoint_vec[i].intersectionPoint2LineDis;
			float d2 = in_my_intersectionPoint_vec[j].intersectionPoint2LineDis;
			float dis = sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
			if (dis<20)   //将距离很小的交点滤除掉160*120
			{
				if (d1<d2)
					staus[j] = 0;
				else
					staus[i] = 0;
			}
		}
	}
	for (int i=0;i<size;i++)
	{
		if(staus[i])
			out_my_intersectionPoint_vec.push_back(in_my_intersectionPoint_vec[i]);
	}
}


Point2f trackedIntersectionPointFlow(Mat img,vector<intersectionPoint_s> preIntersectionPointVec,vector<intersectionPoint_s> curIntersectionPointVec,int &tracknum)
{
	int pre_size = preIntersectionPointVec.size();
	int cur_size = curIntersectionPointVec.size();
	Point2f out_temppoint(0.0,0.0);
	trackedPoint tracked;
	vector<trackedPoint>tracked_vec(0);

	for (int i=0;i<cur_size;i++)
	{
		Point2f p1 = curIntersectionPointVec[i].intersectionPoint;
		
		for (int j=0;j<pre_size;j++)
		{
			Point2f p2 = preIntersectionPointVec[j].intersectionPoint;
			float dis = sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));

			if(dis<8)    //160*120
			{
				tracked.cur_dis = curIntersectionPointVec[i].intersectionPoint2LineDis;
				tracked.pre_dis = preIntersectionPointVec[j].intersectionPoint2LineDis;
				tracked.cur_intersection = p1;
				tracked.pre_intersection = p2;
				tracked_vec.push_back(tracked);
			}
		}
	}

	tracknum = tracked_vec.size();

	if(tracknum>0)
	{
		float mindis = 100000000.0;
		int index = 0;
		for (int k=0;k<tracknum;k++)
		{
			Point2f p=Point2f(tracked_vec[k].cur_intersection.x,tracked_vec[k].cur_intersection.y);
			float dis2Centre = sqrt((p.x-img.cols/2)*(p.x-img.cols/2)+ (p.y-img.rows/2)*(p.y-img.rows/2));
			float dis = 0.3*abs(tracked_vec[k].cur_dis+tracked_vec[k].pre_dis)+0.7*dis2Centre;
			if(dis<mindis)
			{
				mindis = dis;
				index = k;
			}
		}
		circle(img, tracked_vec[index].cur_intersection,3,Scalar(0,0,255),5);  //最终用于定位的点（红色点）
		out_temppoint.x = tracked_vec[index].cur_intersection.x - tracked_vec[index].pre_intersection.x;
		out_temppoint.y = tracked_vec[index].cur_intersection.y - tracked_vec[index].pre_intersection.y;
	}
	return out_temppoint;
}


float twoLineDis(Vec4i line1,Vec4i line2,Point2f intersection)   //交点到两线段最近端点的距离和,若点在线段内，则到该线段的距离为零
{
	float x1 = line1[0];
	float x2 = line1[2];
	float y1 = line1[1];
	float y2 = line1[3];

	float x3 = line2[0];
	float x4 = line2[2];
	float y3 = line2[1];
	float y4 = line2[3];

	float x = intersection.x;
	float y = intersection.y;
	float l1,l2;
	if(x>x1&&x<x2)
	{
		l1 = 0.0;		
	}
	else
	{
		float dis1 = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
		float dis2 = sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2));
		l1 =  min(dis1,dis2);
	}

	if(x>x3&&x<x4)
	{
		l2 = 0.0;
	}
	else
	{
		float dis3 = sqrt((x-x3)*(x-x3)+(y-y3)*(y-y3));
		float dis4 = sqrt((x-x4)*(x-x4)+(y-y4)*(y-y4));
		l2 =  min(dis3,dis4);
	}

	return l1+l2;
}

void myLineMerge(vector<vector<Vec4i>> in_line_clustered,vector<Vec4i> &out_merged_line)
{	
	for (int i=0;i<in_line_clustered.size();i++)
	{
		vector<Vec4i>tmp = in_line_clustered[i];
		double line_length = 0.0;
		double centrex = 0.0;
		double centrey = 0.0;
		double theta = 0.0;
		float vec_x = 0.0;
		float vec_y = 0.0;

		float vec_x_pre = 0.0;
		float vec_y_pre = 0.0;

		float vec_x_res = 0.0;
		float vec_y_res = 0.0;
		
		for (int j=0;j<tmp.size();j++)
		{
			Point pt1(tmp[j][0],tmp[j][1]);
			Point pt2(tmp[j][2],tmp[j][3]);
			double theta_tmp;
			static double tp=0.0;

			line_length += sqrt((pt1.x-pt2.x)*(pt1.x-pt2.x)+(pt1.y-pt2.y)*(pt1.y-pt2.y));
			centrex += (pt1.x+pt2.x)/2.0;
			centrey += (pt1.y+pt2.y)/2.0;

			vec_x = (pt1.x-pt2.x);
			vec_y = (pt1.y-pt2.y);
			if(j>0)
			{
				double len = sqrt((vec_x*vec_x+vec_y*vec_y)*(vec_x_pre*vec_x_pre+vec_y_pre*vec_y_pre));
				double s = (vec_x*vec_x_pre+vec_y*vec_y_pre)/len;
				double s1 = (-vec_x*vec_x_pre-vec_y*vec_y_pre)/len;
				double theta = acos(s);
				double theta1 = acos(s1);
				if(theta1<theta)
				{
					vec_x = -vec_x;
					vec_y = -vec_y;
				}
			}

			if(j==0)
			{
				vec_x_pre = vec_x;
				vec_y_pre = vec_y;
			}

			vec_x_res += vec_x;
			vec_y_res += vec_y;
		}

		if(line_length>50)
		{
			Vec4i tmp1;
			centrex /= tmp.size();
			centrey /= tmp.size();

			double a = double(vec_x_res)/double(sqrt(vec_x_res*vec_x_res+vec_y_res*vec_y_res));
			double b = double(vec_y_res)/double(sqrt(vec_x_res*vec_x_res+vec_y_res*vec_y_res));

			tmp1[0] = cvRound(centrex+1000*(a));
			tmp1[1] = cvRound(centrey+1000*(b));
			tmp1[2] = cvRound(centrex-1000*(a));
			tmp1[3] = cvRound(centrey-1000*(b));

			out_merged_line.push_back(tmp1);
		}
	}
}

Point2f motionEstimate(Mat src,vector<historyPointVec> intersectionPointVecVec,vector<Point2f>historyFlow,int frame_count,int &tracknum)
{
	Point2f temp(0.0,0.0);
	tracknum=0;
	int history_num = historyFlow.size();
	vector<intersectionPoint_s> current_intersectionPoint_vec = intersectionPointVecVec.back().filter_intersectionPoint_vec;

	for (int i=intersectionPointVecVec.size()-1;i>-1;i--)
	{
		int number = frame_count - intersectionPointVecVec[i].index;
		switch (number)
		{
		case 1:  //当前帧数与index的差为1，意思也就是上一帧，如果跟踪到点，给i一个很小值，使for循环不再执行了，输出结果
			temp = trackedIntersectionPointFlow(src,intersectionPointVecVec[i].filter_intersectionPoint_vec,current_intersectionPoint_vec,tracknum);
			if(tracknum>0)
				i=-10000;
			break;
		case 2: //当前帧数与index的差为2，意思也就是上上一帧，如果跟踪到点，给i一个很小值，使for循环不再执行了，输出结果，但此时的结果应该减去上一帧的光流量
			temp = trackedIntersectionPointFlow(src,intersectionPointVecVec[i].filter_intersectionPoint_vec,current_intersectionPoint_vec,tracknum);
			if(tracknum>0)
			{
				temp.x = temp.x - historyFlow[history_num-1].x;
				temp.y = temp.y - historyFlow[history_num-1].y;
				i=-10000;
			}
			break;
		case 3: //当前帧数与index的差为3，意思也就是上上上一帧
			temp = trackedIntersectionPointFlow(src,intersectionPointVecVec[i].filter_intersectionPoint_vec,current_intersectionPoint_vec,tracknum);
			if(tracknum>0)
			{
				temp.x = temp.x - historyFlow[history_num-1].x-historyFlow[history_num-2].x;
				temp.y = temp.y - historyFlow[history_num-1].y-historyFlow[history_num-2].x;
				i=-10000;
			}
			break;
		case 4: //当前帧数与index的差为4，意思也就是上上上上一帧
			temp = trackedIntersectionPointFlow(src,intersectionPointVecVec[i].filter_intersectionPoint_vec,current_intersectionPoint_vec,tracknum);
			if(tracknum>0)
			{
				temp.x = temp.x - historyFlow[history_num-1].x-historyFlow[history_num-2].x-historyFlow[history_num-3].x;
				temp.y = temp.y - historyFlow[history_num-1].y-historyFlow[history_num-2].x-historyFlow[history_num-3].x;
				i=-10000;
			}
			else     //一直跟踪到前四帧，还是没有跟踪到交点，就认为当前帧的光流为前三帧的均值
			{
				// 						temp.x = (historyFlow[4].x+historyFlow[3].x+historyFlow[2].x)/3;
				// 						temp.y = (historyFlow[4].y+historyFlow[3].y+historyFlow[2].y)/3;
			}
			break;

		default:
			break;
		}

	}
	return temp;
}

int lineDetect_run(Mat src,Line_info_t *line_info_tmp)
{
	static int frame_count=0;
	static vector<historyPointVec> intersectionPointVecVec(0);
	static vector<Point2f>historyFlow(0);
	static int line_count=0,track_count = 0;
	int detect_staus = 1;
	static int tracknum = 0;

	//第一步：图像的预处理及二值化
	Mat src_gray,equal_dst;
	cvtColor(src,src_gray,CV_BGR2GRAY);
	equalizeHist(src_gray,equal_dst);
	GaussianBlur(equal_dst,equal_dst,Size(5,5),0,0); 
	adaptiveThreshold(equal_dst, equal_dst, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY_INV,7,7);

	//第二步：基于二值图像的直线检测
	vector<Vec4i> lines;       //cost 5ms
	HoughLinesP(equal_dst,lines,1,CV_PI/180,40,25,10);//最小点数为70，,线条不短于30，间隙不大于10  cost 5ms_640*480

	//第三步（1）对霍夫变换检测的直线进行聚类
	vector<vector<Vec4i>> out_line_clustered;
	getMyClusterLine(lines,out_line_clustered);        

	//第三步(2)：对聚类后的直线进行滤波，使得每个类只有一条直线，并求出线段对应垂直直线的交点坐标 
	vector<intersectionPoint_s> my_intersectionPoint_vec(0);   
	filterLine(src,out_line_clustered,my_intersectionPoint_vec);   

	//第三步(3)：对求得的垂直线的交点坐标进行滤波，使得每个交点都是独立的交点   
	vector<intersectionPoint_s> filter_intersectionPoint_vec(0);
	intersectionPointFilter(my_intersectionPoint_vec,filter_intersectionPoint_vec);  

	for (int k=0;k<filter_intersectionPoint_vec.size();k++)
	{
		circle(src,filter_intersectionPoint_vec[k].intersectionPoint,3,Scalar(255,0,0),5); //当前帧检测到的垂直交点
	}

	//第四步:利用最近的4帧的交点，估计全局的运动量
	historyPointVec historyPoint;
	historyPoint.filter_intersectionPoint_vec = filter_intersectionPoint_vec;
	historyPoint.index = frame_count;
	intersectionPointVecVec.push_back(historyPoint);
	if(intersectionPointVecVec.size()>5)
	{
		intersectionPointVecVec.erase(intersectionPointVecVec.begin());
	}

	Point2f temp(0.0,0.0);
	if(frame_count>0&&intersectionPointVecVec.size()>0)
	{
		temp = motionEstimate(src,intersectionPointVecVec,historyFlow,frame_count,tracknum);
	}

	//第五步：反馈报警，并输出估计的全局运动量

	//反馈报警类型一：聚类直线太多
	if(out_line_clustered.size()>7)
		line_count++;
	else
	{
		if(line_count>0)
		{
			line_count--;
		}
		else
			line_count = 0;
	}

	if(line_count>4)
	{
		detect_staus = 0;
		if(line_count>7)
			line_count = 8;
	}

	//反馈报警消息类型二，多帧没有跟踪到交点
	if(tracknum == 0)
	{
		track_count++;
	}
	else
	{
		track_count>0?track_count--:track_count=0;
	}
	if(track_count>50)
	{
		detect_staus = 0;
		if(track_count>65)
		track_count = 65;
	}

	//输出估计的全局运动量
	line_info_tmp->pixel_flow_x = temp.x;
	line_info_tmp->pixel_flow_y = temp.y;

	//保留历史的光流信息
	historyFlow.push_back(temp);
	if(historyFlow.size()>5)
	{
		historyFlow.erase(historyFlow.begin());
	}

	frame_count++;

	return detect_staus;
}