
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

void getMyClusterLine(vector<Vec4i> lines,vector<vector<Vec4i>> &out_line_cluster); //ֱ�ߵľ���

void myLineMerge(vector<vector<Vec4i>> in_line_clustered,vector<Vec4i> &out_merged_line);  //������ֱ�ߵ���ϣ���ʱ��û�õ���

void mycvGetPoint(Vec4i line1, Vec4i line2,Point &out_inter_point);   //����ֱ�߽���

float intersection_theta(Vec4i line1, Vec4i line2);  //����ֱ�ߵļнǣ���ǣ�

void filterLine(Mat img,vector<vector<Vec4i>> out_line_cluster,vector<intersectionPoint_s> &my_intersectionPoint_vec_tmp);  //��ֱ�߽����˲�

float twoLineDis(Vec4i line1,Vec4i line2,Point2f intersection);  //������ֱ�ཻ�߶ε�����ľ����

//��ǰ֡������˲�����Ҫ��ȥ������ܽ��ĵ�
void intersectionPointFilter(vector<intersectionPoint_s> in_my_intersectionPoint_vec,vector<intersectionPoint_s>&out_my_intersectionPoint_vec);

//�õ����յ�ȫ���˶�����flow_x,flow_y
Point2f trackedIntersectionPointFlow(Mat img,vector<intersectionPoint_s> preIntersectionPointVec,vector<intersectionPoint_s> curIntersectionPointVec,int &tracknum);

void addHistoryPoint(vector<vector<intersectionPoint_s>> intersectionPointVecVec,vector<vector<intersectionPoint_s>> &out_temp_vec);
void addHistoryPoint1(vector<vector<intersectionPoint_s>> intersectionPointVecVec,vector<historyPointVec> &out_temp_vec);

Point2f motionEstimate(Mat src,vector<historyPointVec> intersectionPointVecVec,vector<Point2f>historyFlow,int frame_count,int &tracknum);

int lineDetect_run(Mat Image,Line_info_t *line_info_tmp);
