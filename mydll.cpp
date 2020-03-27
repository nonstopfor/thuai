#include <stdio.h>
#include "mydll.h"
#include "definition.h"
#include <math.h>
#include <time.h>
#include <vector>
#include <utility>
#include <algorithm>
#define MI 23
#define NI 19991
using std::pair;
using std::make_pair;
using std::max;
typedef pair<double, double> PAIR;

/*
�������������ָ�������
info.myCommandList.addCommand(Division,aim_cell_id,direction);//��������ڶ��������Ƿ��ѷ���
info.myCommandList.addCommand(Move,aim_cell_id,direction);//�ƶ�����ڶ����������ƶ�����
info.myCommandList.addCommand(Spit,aim_cell_id,direction);//��������ڶ������������·���
*/
#include<iostream>
#include<cmath>
#include<vector>
using namespace std;

double dist(double x1, double y1, double x2, double y2) {
	double deltaX = x1 - x2;
	double deltaY = y1 - y2;
	return sqrt(deltaX * deltaX + deltaY * deltaY);
}

double MINBOUND = 0.32;
double MAXDIST = 425;
int DISASTERROUND = 700;
int splitCheck(std::vector<CellInfo> cells, int maxCell, int curCell, int round) {
	double maxR = cells[maxCell].r;
	if ((cells[curCell].r > MINBOUND * maxR && round < DISASTERROUND) ||
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

bool Judis(PAIR P1, PAIR P2, PAIR yuan, double R) {
	double A, B, C, dist1, dist2, angle1, angle2;//Ax+By+C=0;//(y1-y2)x +(x2-x1)y +x1y2-y1x2=0
	if (P1.first == P2.first)
		A = 1, B = 0, C = -P1.first;
	else if (P1.second == P2.second)
		A = 0, B = 1, C = -P1.second;
	else
	{
		A = P1.second - P2.second;
		B = P2.first - P1.first;
		C = P1.first * P2.second - P1.second * P2.first;
	}
	dist1 = A * yuan.first + B * yuan.second + C;
	dist1 *= dist1;
	dist2 = (A * A + B * B) * R * R;
	if (dist1 > dist2) return false;//�㵽ֱ�߾�����ڰ뾶r  ���ཻ
	angle1 = (yuan.first - P1.first) * (P2.first - P1.first) + (yuan.second - P1.second) * (P2.second - P1.second);
	angle2 = (yuan.first - P2.first) * (P1.first - P2.first) + (yuan.second - P2.second) * (P1.second - P2.second);
	if (angle1 > 0 && angle2 > 0) return true;//���Ҷ�Ϊ����������� �ཻ
	return false;//���ཻ
	
}
bool safe(Info& info, double x1, double y1, double r, double x2, double y2) {
	TPlayerID myid = info.myID;
	auto p1 = make_pair(x1, y1);
	auto p2 = make_pair(x2, y2);

	for (int i = 0; i < info.cellInfo.size(); ++i) {
		auto& cell = info.cellInfo[i];
		if (cell.ownerid == myid) continue;
		if (cell.r <= r) continue;
		double x = cell.x, y = cell.y;
		double ar = cell.r;
		auto p3 = make_pair(x, y);
		if (Judis(p1, p2, p3, max(0.0, ar - r))) return false;
	}
	for (int i = 0; i < info.spikyballInfo.size(); ++i) {
		auto& t = info.spikyballInfo[i];
		double x = t.sx, y = t.sy;
		double ar = t.sr;
		auto p3 = make_pair(x, y);
		if (Judis(p1, p2, p3, max(0.0,ar-r))) return false;
	}
	return true;
}

void player_ai(Info& info)
{
	vector<CellInfo> myCell;
	int maxCell = 0;
	int MincellNum = -1;
	int MincellRadius = info.cellInfo[0].r;
	int Mincell_x, Mincell_y = -1;
	TPlayerID myid = info.myID;
	for (int i = 0; i < info.cellInfo.size(); i++)
		if (info.cellInfo[i].ownerid == myid)
			myCell.push_back(info.cellInfo[i]);

	for(int i=0;i<myCell.size();i++)
		if(myCell[i].r>myCell[maxCell].r)
			maxCell = i;

	for(int i=0;i<myCell.size();i++)
	{
		int split = splitCheck(myCell, maxCell, i, info.round);
		int targetX=10000,targetY=10000;
		if(split!=-1){
			targetX = myCell[split].x;
			targetY = myCell[split].y;
		}else{
			for(auto k:info.nutrientInfo){
				if((targetX-myCell[i].x)*(targetX-myCell[i].x)+(targetY-myCell[i].y)*(targetY-myCell[i].y)<
					(k.nux-myCell[i].x)*(k.nux-myCell[i].x)+(k.nuy-myCell[i].y)*(k.nuy-myCell[i].y)){
					if(safe(info,myCell[i].x,myCell[i].y,myCell[i].r,k.nux,k.nuy)){
						targetX = k.nux;
						targetY = k.nuy;
					}
				}
			}
			for(auto k:info.cellInfo){
				if(k.ownerid!=myid && k.r<myCell[i].r && (targetX-myCell[i].x)*(targetX-myCell[i].x)+(targetY-myCell[i].y)*(targetY-myCell[i].y)<
					(k.x-myCell[i].x)*(k.x-myCell[i].x)+(k.y-myCell[i].y)*(k.y-myCell[i].y)){
					if(safe(info,myCell[i].x,myCell[i].y,myCell[i].r,k.x,k.y)){
						targetX = k.x;
						targetY = k.y;
					}
				}
			}
		}
		int direction = 0;
		double pi = 3.14159265;
		int dx, dy = 0;
		if (targetX != 10000)
		{
			dx = targetX - myCell[i].x;
			dy = targetY - myCell[i].y;
			direction = (int)(atan2(dy, dx) / pi * 180 + 360) % 360;
			info.myCommandList.addCommand(Move, myCell[i].id, direction);
		}
		else
		{
			dx = 150 - myCell[i].x;
			dy = 150 - myCell[i].y;
			direction = (int)(atan2(dy, dx) / pi * 180 + 360) % 360;
			info.myCommandList.addCommand(Move, myCell[i].id, direction);
		}
	}

}
