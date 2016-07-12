/****************************************************************************
 *
 *   Copyright (C) 2013 PX4 Development Team. All rights reserved.
 *   Author: Petri Tanskanen <tpetri@inf.ethz.ch>
 *   		 Lorenz Meier <lm@inf.ethz.ch>
 *   		 Samuel Zihlmann <samuezih@ee.ethz.ch>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/



#include "stdafx.h"
#include "flow_linux.h"
#include <stdlib.h>
#include <stdio.h>
//#include <stdbool.h>
#include <math.h>
#include <algorithm> 

using namespace std;
#define sign(x) (( x > 0 ) - ( x < 0 ))

//#include "mavlink_bridge_header.h"
//#include <mavlink.h>
//#include "dcmi.h"
//#include "debug.h"

#define __INLINE inline
#define __ASM asm
//#include "core_cm4_simd.h"
//目前所用图像是 800*540
//#define FRAME_SIZE	200//global_data.param[PARAM_IMAGE_WIDTH]
//#define SEARCH_SIZE	8//global_data.param[PARAM_MAX_FLOW_PIXEL] // maximum offset to search: 4 + 1/2 pixels
#define SEARCH_SIZE	8//global_data.param[PARAM_MAX_FLOW_PIXEL] // maximum offset to search: 4 + 1/2 pixels
#define TILE_SIZE	8               						// x & y tile size
#define NUM_BLOCKS	10// x & y number of tiles to check
#define BOTTOM_FLOW_FEATURE_THRESHOLD  30 //需要确定
//#define IMAGE_WIDTH  FRAME_SIZE //图像宽度需要确定
#define BOTTOM_FLOW_VALUE_THRESHOLD  400//需要确定
#define  BOTTOM_FLOW_HIST_FILTER 1//启用柱状图分析
#define  FOCAL_LENGTH_MM 12
#define  BOTTOM_FLOW_GYRO_COMPENSATION 0//启用补偿
#define  get_time_between_images() 100 //帧时间间隔
#define  GYRO_COMPENSATION_THRESHOLD 10//角速度阈值
#define FRAME_SIZE 120
#define IMAGE_WIDTH 160

/**
 * @brief Compute the average pixel gradient of all horizontal and vertical steps
 *
 * TODO compute_diff is not appropriate for low-light mode images
 *
 * @param image ...
 * @param offX x coordinate of upper left corner of 8x8 pattern in image
 * @param offY y coordinate of upper left corner of 8x8 pattern in image
 */
static unsigned int compute_diff(unsigned char *image, unsigned short offX, unsigned short offY, unsigned short row_size)
{
	/* calculate position in image buffer */
	unsigned int off = (offY + 2) * row_size + (offX + 2); // we calc only the 4x4 pattern
	unsigned int acc;
	unsigned int i,j;

	/*计算行梯度*/             
	acc = 0;
	for(i =0; i < 4; i++)
	{
		for(j = 0; j < 4; j ++)
		{
			acc += abs((int)(image[off + j + i * row_size]) - (int)(image[off + j + (i + 1)* row_size]));
		}
	}

	/*计算列梯度*/

	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j ++)
		{
			acc += abs(image[off + j + i * row_size] - image[off + j +1 + i* row_size]);
		}
	}

	return acc;
}

/**
 * @brief Compute SAD distances of subpixel shift of two 8x8 pixel patterns.
 *
 * @param image1 ...
 * @param image2 ...
 * @param off1X x coordinate of upper left corner of pattern in image1
 * @param off1Y y coordinate of upper left corner of pattern in image1
 * @param off2X x coordinate of upper left corner of pattern in image2
 * @param off2Y y coordinate of upper left corner of pattern in image2
 * @param acc array to store SAD distances for shift in every direction
 */
static  int compute_subpixel(unsigned char *image1, unsigned char *image2, unsigned short off1X, unsigned short off1Y, unsigned short off2X, unsigned short off2Y, unsigned int *acc, unsigned short row_size)
{
	/* calculate position in image buffer */
	unsigned int off1 = off1Y * row_size + off1X; // image1
	unsigned int off2 = off2Y * row_size + off2X; // image2

	unsigned int s0, s1, s2, s3, s4, s5, s6, s7, t1, t3, t5, t7;
	unsigned int gray;
	unsigned short i;
	unsigned short j;

	for (i = 0; i < 8; i++)
	{
		acc[i] = 0;
	}


	/*
	 * calculate for each pixel in the 8x8 field with upper left corner (off1X / off1Y)
	 * every iteration is one line of the 8x8 field.
	 *
	 *  + - + - + - + - + - + - + - + - +
	 *  |   |   |   |   |   |   |   |   |
	 *  + - + - + - + - + - + - + - + - +
	 *
	 *
	 */

	for (i = 0; i < 4; i++)
	{

		for(j = 0; j < 4; j ++)
		{
			s0 = (image2[off2 +  j  + 0 + (i+0) * row_size] + image2[off2 + j + 1 + (i+0) * row_size])/2;
			s1 = (image2[off2 +  j  + 0 + (i+1) * row_size] + image2[off2 + j + 1 + (i+1) * row_size])/2;
			s2 = (image2[off2 +  j  + 0 + (i+0) * row_size] + image2[off2 + j + 0 + (i+1) * row_size])/2;
			s3 = (image2[off2 +  j  + 0 + (i+1) * row_size] + image2[off2 + j - 1 + (i+1) * row_size])/2;
			s4 = (image2[off2 +  j  + 0 + (i+0) * row_size] + image2[off2 + j - 1 + (i+0) * row_size])/2;
			s5 = (image2[off2 +  j  + 0 + (i-1) * row_size] + image2[off2 + j - 1 + (i-1) * row_size])/2;
			s6 = (image2[off2 +  j  + 0 + (i+0) * row_size] + image2[off2 + j + 0 + (i-1) * row_size])/2;
			s7 = (image2[off2 +  j  + 0 + (i-1) * row_size] + image2[off2 + j + 1 + (i-1) * row_size])/2;
			t1 = (s0 + s1)/2;
			t3 = (s3 + s4)/2;
			t5 = (s4 + s5)/2;
			t7 = (s7 + s0)/2;
			gray = image1[off1 + j + 0 + i * row_size];
			acc[0] += abs((int)(gray) - (int)(s0));
			acc[1] += abs((int)(gray) - (int)(t1));
			acc[2] += abs((int)(gray) - (int)(s2));
			acc[3] += abs((int)(gray) - (int)(t3));
			acc[4] += abs((int)(gray) - (int)(s4));
			acc[5] += abs((int)(gray) - (int)(t5));
			acc[6] += abs((int)(gray) - (int)(s6));
			acc[7] += abs((int)(gray) - (int)(t7));
		}
	}

	return 0;
}

/**
 * @brief Compute SAD of two 8x8 pixel windows.
 *
 * @param image1 ...
 * @param image2 ...
 * @param off1X x coordinate of upper left corner of pattern in image1
 * @param off1Y y coordinate of upper left corner of pattern in image1
 * @param off2X x coordinate of upper left corner of pattern in image2
 * @param off2Y y coordinate of upper left corner of pattern in image2
 */
static unsigned int compute_sad_8x8(unsigned char *image1, unsigned char *image2, unsigned short off1X, unsigned short off1Y, unsigned short off2X, unsigned short off2Y, unsigned short row_size)
{
	/* calculate position in image buffer */
	unsigned int off1 = off1Y * row_size + off1X; // image1
	unsigned int off2 = off2Y * row_size + off2X; // image2

	int image1_p,image2_p;

	unsigned int acc;
	unsigned short i ;
	unsigned short j;

	acc = 0;
	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j++)
		{
			acc += abs((int)(image1[off1 + i + j * row_size]) - (int)(image2[off2 + i + j * row_size]));
		}
	}
	return acc;
}

/**
 * @brief Computes pixel flow from image1 to image2
 *
 * Searches the corresponding position in the new image (image2) of max. 64 pixels from the old image (image1)
 * and calculates the average offset of all.
 *
 * @param image1 previous image buffer
 * @param image2 current image buffer (new)
 * @param x_rate gyro x rate
 * @param y_rate gyro y rate
 * @param z_rate gyro z rate
 *
 * @return quality of flow calculation
 */
unsigned char compute_flow2(unsigned char * image1, unsigned char * image2, float x_rate, float y_rate, float z_rate, float *pixel_flow_x, float *pixel_flow_y) {

	/* constants */
	const short winmin = -SEARCH_SIZE;
	const short winmax = SEARCH_SIZE;
	const unsigned short hist_size = 2*(winmax-winmin+1)+1;

	/* variables */
    unsigned short pixLo = 2*SEARCH_SIZE + 1;
    unsigned short pixHi = FRAME_SIZE - (2*SEARCH_SIZE + 1) - TILE_SIZE;
	unsigned short pixLo1 = 3*SEARCH_SIZE + 1;
	unsigned short pixHi1 = IMAGE_WIDTH - (3*SEARCH_SIZE + 1) - TILE_SIZE;
	unsigned short pixStep1 = (pixHi1 - pixLo1) / (NUM_BLOCKS+2)+1;
    unsigned short pixStep = (pixHi - pixLo) / NUM_BLOCKS;
	unsigned short i, j;
	unsigned int acc[8]; // subpixels
	unsigned short histx[hist_size]; // counter for x shift
	unsigned short histy[hist_size]; // counter for y shift
	signed char  dirsx[800]; // shift directions in x
	signed char  dirsy[800]; // shift directions in y
	unsigned char  subdirs[800]; // shift directions of best subpixels
	float meanflowx = 0.0f;
	float meanflowy = 0.0f;
	unsigned int meancount = 0;
	float histflowx = 0.0f;
	float histflowy = 0.0f;
	unsigned int dir  = 0;
	unsigned char qual; 

	/* initialize with 0 */
	for (j = 0; j < hist_size; j++) { histx[j] = 0; histy[j] = 0;}
	memset(dirsx,0,sizeof(dirsx));
	memset(dirsy,0,sizeof(dirsy));
	memset(subdirs,0,sizeof(subdirs));

	/* iterate over all patterns
	 */
	for (j = pixLo; j < pixHi; j += pixStep)
	{
		for (i = pixLo1; i < pixHi1; i += pixStep1)
		{
			unsigned int diff;
			unsigned int dist = 0xFFFFFFFF; // set initial distance to "infinity"
			signed char sumx = 0;
			signed char sumy = 0;
			signed char ii, jj;
			/* test pixel if it is suitable for flow tracking */
			diff = compute_diff(image1, i, j, IMAGE_WIDTH);//求出行列间总梯度4*4  
			
			if (diff < BOTTOM_FLOW_FEATURE_THRESHOLD)     //也就是找纹理丰富的位置
			{
				//小于阈值就跳过
				continue;
			}
			for (jj = winmin; jj <= winmax; jj++)
			{
				for (ii = winmin; ii <= winmax; ii++)
				{
					unsigned int temp_dist;
        		    temp_dist = compute_sad_8x8(image1, image2, i, j, i + ii, j + jj, IMAGE_WIDTH);
					if (temp_dist < dist)
					{
						sumx = ii;
						sumy = jj;
						dist = temp_dist;
					}
				}
			}
			
			if (dist < BOTTOM_FLOW_VALUE_THRESHOLD)  //SAD dist的接受阈值：为400
			{
				unsigned int mindist;
				unsigned char mindir;
				unsigned char k;
				unsigned char hist_index_x;
				unsigned char hist_index_y ;

				meanflowx += (float) sumx;
				meanflowy += (float) sumy;

				compute_subpixel(image1, image2, i, j, i + sumx, j + sumy, acc, IMAGE_WIDTH);
				mindist = dist; // best SAD until now
				mindir = 8; // direction 8 for no direction
				for( k = 0; k < 8; k++)
				{
					if (acc[k] < mindist)
					{
						// SAD becomes better in direction k
						mindist = acc[k];
						mindir = k;
						dir = mindir;
					}
				}

				dirsx[meancount] = sumx;
				dirsy[meancount] = sumy;
				subdirs[meancount] = mindir;

				/* feed histogram filter*/
			    hist_index_x = 2*sumx + (winmax-winmin+1);  // 2*SEARCH_SIZE + 1  若没有偏移，hist柱为17，向左偏1为16，向右为18
				if (subdirs[meancount] == 0 || subdirs[meancount] == 1 || subdirs[meancount] == 7) hist_index_x += 1;
				if (subdirs[meancount] == 3 || subdirs[meancount] == 4 || subdirs[meancount] == 5) hist_index_x += -1;
				hist_index_y = 2*sumy + (winmax-winmin+1);
				if (subdirs[meancount] == 5 || subdirs[meancount] == 6 || subdirs[meancount] == 7) hist_index_y += -1;
				if (subdirs[meancount] == 1 || subdirs[meancount] == 2 || subdirs[meancount] == 3) hist_index_y += 1;

				histx[hist_index_x]++;
				histy[hist_index_y]++;
				meancount++;
			}
		}
	}

	//******************主for循环结束***************************

	/* evaluate flow calculation *///X,Y偏移求平均
	if (meancount > 6)     //符合条件的点大于10个
	{
		short maxpositionx = 0;
		short maxpositiony = 0;
		unsigned short maxvaluex = 0;
		unsigned short maxvaluey = 0;

		meanflowx /= meancount;
		meanflowy /= meancount;

		/* position of maximal histogram peek *///寻找 X,Y分布最大
		for (j = 0; j < hist_size; j++)
		{
			if (histx[j] > maxvaluex)
			{
				maxvaluex = histx[j];
				maxpositionx = j;
			}
			if (histy[j] > maxvaluey)
			{
				maxvaluey = histy[j];
				maxpositiony = j;
			}
		}

		/* check if there is a peak value in histogram *///检查柱状图中的估计值 ( 概率统计)
		if (1) //(histx[maxpositionx] > meancount / 6 && histy[maxpositiony] > meancount / 6)
		{
			//if (global_data.param[PARAM_BOTTOM_FLOW_HIST_FILTER])
			if (BOTTOM_FLOW_HIST_FILTER)             //else直接平均
			{
				/* use histogram filter peek value */
				unsigned short hist_x_min = maxpositionx;
				unsigned short hist_x_max = maxpositionx;
				unsigned short hist_y_min = maxpositiony;
				unsigned short hist_y_max = maxpositiony;

				float hist_x_value = 0.0f;
				float hist_x_weight = 0.0f;

				float hist_y_value = 0.0f;
				float hist_y_weight = 0.0f;

				unsigned short i;

				/* x direction */
				if (maxpositionx > 1 && maxpositionx < hist_size-2)
				{
					hist_x_min = maxpositionx - 2;
					hist_x_max = maxpositionx + 2;
				}
				else if (maxpositionx == 0)
				{
					hist_x_min = maxpositionx;
					hist_x_max = maxpositionx + 2;
				}
				else  if (maxpositionx == hist_size-1)
				{
					hist_x_min = maxpositionx - 2;
					hist_x_max = maxpositionx;
				}
				else if (maxpositionx == 1)
				{
					hist_x_min = maxpositionx - 1;
					hist_x_max = maxpositionx + 2;
				}
				else  if (maxpositionx == hist_size-2)
				{
					hist_x_min = maxpositionx - 2;
					hist_x_max = maxpositionx + 1;
				}

				/* y direction */
				if (maxpositiony > 1 && maxpositiony < hist_size-2)
				{
					hist_y_min = maxpositiony - 2;
					hist_y_max = maxpositiony + 2;
				}
				else if (maxpositiony == 0)
				{
					hist_y_min = maxpositiony;
					hist_y_max = maxpositiony + 2;
				}
				else if (maxpositiony == hist_size-1)
				{
					hist_y_min = maxpositiony - 2;
					hist_y_max = maxpositiony;
				}
				else if (maxpositiony == 1)
				{
					hist_y_min = maxpositiony - 1;
					hist_y_max = maxpositiony + 2;
				}
				else if (maxpositiony == hist_size-2)
				{
					hist_y_min = maxpositiony - 2;
					hist_y_max = maxpositiony + 1;
				}

				//MIN 小2  MAX 大2

				for (i = hist_x_min; i < hist_x_max+1; i++)
				{
					hist_x_value += (float) (i*histx[i]);//面积
					hist_x_weight += (float) histx[i];//高度                   //这里的结果可以用于返回的质量评价！！！！！！！
				}

				for (i = hist_y_min; i<hist_y_max+1; i++)
				{
					hist_y_value += (float) (i*histy[i]);
					hist_y_weight += (float) histy[i];
				}
                
				histflowx = (hist_x_value/hist_x_weight - (winmax-winmin+1)) / 2.0f ;
				histflowy = (hist_y_value/hist_y_weight - (winmax-winmin+1)) / 2.0f;
			}
			else
			{
				/* use average of accepted flow values *///直接平均值
				unsigned int meancount_x = 0;
				unsigned int meancount_y = 0;
				unsigned short i;

				for (i = 0; i < meancount; i++)
				{
					float subdirx = 0.0f;
					float subdiry = 0.0f;

					if (subdirs[i] == 0 || subdirs[i] == 1 || subdirs[i] == 7) subdirx = 0.5f;  //这里就包含了可以使用亚像素，因为如果亚像素匹配低的话，这里的subdirs[i]=8
					if (subdirs[i] == 3 || subdirs[i] == 4 || subdirs[i] == 5) subdirx = -0.5f;
					histflowx += (float)dirsx[i] + subdirx;
					meancount_x++;

					if (subdirs[i] == 5 || subdirs[i] == 6 || subdirs[i] == 7) subdiry = -0.5f;
					if (subdirs[i] == 1 || subdirs[i] == 2 || subdirs[i] == 3) subdiry = 0.5f;
					histflowy += (float)dirsy[i] + subdiry;
					meancount_y++;
				}

				histflowx /= meancount_x;
				histflowy /= meancount_y;
			}

			if (BOTTOM_FLOW_GYRO_COMPENSATION)
			{

			} else
			{
				/* without gyro compensation */
				*pixel_flow_x = histflowx;
				*pixel_flow_y = histflowy;
			}

		}
		else                           //if(1)
		{
			*pixel_flow_x = 0.0f;
			*pixel_flow_y = 0.0f;
			return 0;
		}
	}
	else                  //count<10
	{
		*pixel_flow_x = 0.0f;
		*pixel_flow_y = 0.0f;
		return 0;
	}

	/* calc quality */
	qual = (unsigned char)(meancount * 255 / (NUM_BLOCKS*NUM_BLOCKS+44));

	return qual;
}


int flow_run(unsigned char* _ucImage,Flow_info_t *flow_info_tmp)
{
	static int nFrmNum=0;
	unsigned char return_res = 0;

	int sad_flag = 1;
	static int sad_not_fit_count=0;

	float x_rate=0.0, y_rate=0.0, z_rate=0.0,pixel_flow_x_t=0.0,pixel_flow_y_t=0.0;

	if(nFrmNum>0)
	{
		//主处理函数，第一个参数为上一帧图像数据，第二个为当前帧图像数据，pixel_flow_x，pixel_flow_y为输出结果
		return_res = compute_flow2(pre_img, _ucImage, x_rate, y_rate, z_rate, &pixel_flow_x_t, &pixel_flow_y_t);
	}
	flow_info_tmp->pixel_flow_x = pixel_flow_x_t;
	flow_info_tmp->pixel_flow_y = pixel_flow_y_t;

	memcpy(pre_img,_ucImage,160*120*sizeof(unsigned char));

	nFrmNum++;

	//报警消息
	if((int)return_res<140)
		sad_not_fit_count++;
	else
	{
		if(sad_not_fit_count>0)
			sad_not_fit_count--;
		else
			sad_not_fit_count=0;
	}
	if(sad_not_fit_count>100)
	{
		sad_flag=0;
		if(sad_not_fit_count>120)
		        sad_not_fit_count = 120;
	}	
	return sad_flag;
}

