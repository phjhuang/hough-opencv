
#include <opencv2/opencv.hpp> 

using namespace std;
using namespace cv;

typedef struct Line_info_s
{
	float pixel_flow_x ;
	float pixel_flow_y; 
} Line_info_t ;

typedef struct MyIntersectionPoint_s
{
	Point2f intersectionPoint;
	Vec4i Line1;
	Vec4i Line2;
    float intersectionPoint2LineDis;

}intersectionPoint_s;

typedef struct trackedPoint_s
{
	Point2f pre_intersection;
	float pre_dis;
	Point2f cur_intersection;
	float cur_dis;
}trackedPoint;

typedef struct historyPointVec_s
{
	vector<intersectionPoint_s> filter_intersectionPoint_vec;
	int index;
}historyPointVec;

void getMyClusterLine(vector<Vec4i> lines,vector<vector<Vec4i>> &out_line_cluster); //直线的聚类

void myLineMerge(vector<vector<Vec4i>> in_line_clustered,vector<Vec4i> &out_merged_line);  //聚类后的直线的拟合（暂时还没用到）

void mycvGetPoint(Vec4i line1, Vec4i line2,Point &out_inter_point);   //求两直线交点

float intersection_theta(Vec4i line1, Vec4i line2);  //求两直线的夹角（锐角）

void filterLine(Mat img,vector<vector<Vec4i>> out_line_cluster,vector<intersectionPoint_s> &my_intersectionPoint_vec_tmp);  //对直线进行滤波

float twoLineDis(Vec4i line1,Vec4i line2,Point2f intersection);  //求两垂直相交线段到交点的距离和

//当前帧交点的滤波，主要是去除相隔很近的点
void intersectionPointFilter(vector<intersectionPoint_s> in_my_intersectionPoint_vec,vector<intersectionPoint_s>&out_my_intersectionPoint_vec);

//得到最终的全局运动估计flow_x,flow_y
Point2f trackedIntersectionPointFlow(Mat img,vector<intersectionPoint_s> preIntersectionPointVec,vector<intersectionPoint_s> curIntersectionPointVec,int &tracknum);

void addHistoryPoint(vector<vector<intersectionPoint_s>> intersectionPointVecVec,vector<vector<intersectionPoint_s>> &out_temp_vec);
void addHistoryPoint1(vector<vector<intersectionPoint_s>> intersectionPointVecVec,vector<historyPointVec> &out_temp_vec);

Point2f motionEstimate(Mat src,vector<historyPointVec> intersectionPointVecVec,vector<Point2f>historyFlow,int frame_count,int &tracknum);

int lineDetect_run(Mat Image,Line_info_t *line_info_tmp);
