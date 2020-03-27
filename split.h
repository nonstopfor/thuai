/*************************************************************************
	> File Name: split.h
	> Author: kxz
	> Mail: 15068701650@163.com 
	> Created Time: Fri 27 Mar 2020 09:20:05 PM CST
 ************************************************************************/

#include<iostream>
#include<cmath>
#include<vector>
using namespace std;

double dist(double x1, double y1, double x2, double y2) {
	double deltaX = x1 - x2;
	double deltaY = y1 - y2;
	return sqrt(deltaX*deltaX + deltaY*deltaY);
}

double MINBOUND = 0.32;
double MAXDIST = 425;
int DISASTERROUND
int splitCheck(std::vector<CellInfo> cells, int maxCell, int curCell, int round) {
	double maxR = cells[maxCell].r;
	if ((cells[curCell].r > MINBOUND*maxR && round < DISASTERROUND) ||
		(round >= DISASTERROUND && curCell != maxCell))
		return -1;
	double x = cells[curCell].x, y = cells[curCell].y;
	int target = -1;
	double minDist = MAXDIST;
	for (int i = 0; i < cells.size(); ++i) {
		double distance = dist(x, y, cells[i].x, cells[i].y);
		if (i != curCell && distance < minDist) {
			target = i;
			minDist = distance;
		}
	}
	return target;
}
