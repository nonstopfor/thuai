#include <stdio.h>
#include "mydll.h"
#include "definition.h"
#include <time.h>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
#define MI 23
#define NI 19991
const double N = 300;
const double lam = 0.9;//小细胞与大细胞半径比值小于这个时被吞噬
typedef pair<double, double> PAIR;
/*
基本的三个添加指令的命令
info.myCommandList.addCommand(Division,aim_cell_id,direction);//分裂命令，第二个参数是分裂方向
info.myCommandList.addCommand(Move,aim_cell_id,direction);//移动命令，第二个参数是移动方向
info.myCommandList.addCommand(Spit,aim_cell_id,direction);//吞吐命令，第二个参数是吞吐方向
*/


double dist(double x1, double y1, double x2, double y2) {
	double deltaX = x1 - x2;
	double deltaY = y1 - y2;
	return sqrt(deltaX * deltaX + deltaY * deltaY);
}

double distCell(CellInfo c1, CellInfo c2, bool removeRadius = false) {
	double distRaw = dist(c1.x, c1.y, c2.x, c2.y);
	if (removeRadius) distRaw -= c1.r + c2.r;
	return distRaw;
}

double MINBOUND = 0.32;
double MAXDIST = 425;
double DISTFACTOR = 1.1;
int DISASTERROUND = 700;
int splitCheck(std::vector<CellInfo> cells, int maxCell, int curCell,
	int round, CellInfo nearestEnemy) {

	double maxR = cells[maxCell].r;
	double enemyR = nearestEnemy.r;
	bool minboundJudge = cells[curCell].r > MINBOUND * maxR;
	bool roundJudge = round < DISASTERROUND || curCell == maxCell;
	bool enemyJudge = maxR > 0.9 * enemyR;
	if (!enemyJudge) {
		double distToEnemy = dist(
			cells[curCell].x, cells[curCell].y,
			nearestEnemy.x, nearestEnemy.y
		);
		enemyJudge = distToEnemy > DISTFACTOR * (maxR + enemyR); //max cell is still far from enemy
	}
	if (minboundJudge && roundJudge && enemyJudge)
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
	if (dist1 > dist2) return false;//点到直线距离大于半径r  不相交
	angle1 = (yuan.first - P1.first) * (P2.first - P1.first) + (yuan.second - P1.second) * (P2.second - P1.second);
	angle2 = (yuan.first - P2.first) * (P1.first - P2.first) + (yuan.second - P2.second) * (P1.second - P2.second);
	if (angle1 > 0 && angle2 > 0) return true;//余弦都为正，则是锐角 相交
	return false;//不相交

}
bool safe(Info& info, double x1, double y1, double r, double x2, double y2) {
	TPlayerID myID = info.myID;
	auto p1 = make_pair(x1, y1);
	auto p2 = make_pair(x2, y2);

	for (int i = 0; i < info.cellInfo.size(); ++i) {
		auto& cell = info.cellInfo[i];
		if (cell.ownerid == myID) continue;
		if (cell.r <= r) continue;
		double x = cell.x, y = cell.y;
		double ar = cell.r;
		auto p3 = make_pair(x, y);
		if (Judis(p1, p2, p3, ar + r)) return false;
	}
	for (int i = 0; i < info.spikyballInfo.size(); ++i) {
		auto& t = info.spikyballInfo[i];
		double x = t.sx, y = t.sy;
		double ar = t.sr;
		auto p3 = make_pair(x, y);
		if (Judis(p1, p2, p3, ar + r)) return false;
	}
	return true;
}
int compute_dir(double tx, double ty, double sx, double sy) {
	double dx = tx - sx;
	double dy = ty - sy;
	double pi = 3.14159265;
	int direction = (int)(atan2(dy, dx) / pi * 180 + 360) % 360;
	return direction;
}

bool catchable(CellInfo me, CellInfo enemy) {
	double distance = dist(me.x, me.y, enemy.x, enemy.y);
	distance = distance - 2.0 / 3.0 * me.r;
	double dist_hat_dir = compute_dir(enemy.x, enemy.y, me.x, me.y);
	double mySpeed = me.v * cos(me.d - dist_hat_dir);
	double enemySpeed = enemy.v * cos(enemy.v - dist_hat_dir);
	double myAcc = 10/me.r, enAcc = 10/enemy.r, 
		   myTop = 20/me.r;
	double t = (myTop - mySpeed)/myAcc,
		   t_limit = (myTop - enemySpeed)/enAcc;
	if (t >= t_limit) {
		double runDist = (mySpeed - enemySpeed) * t_limit +
						 0.5 * (myAcc - enAcc) * t_limit * t_limit;
		return runDist > distance;
	} else {
		double deltaT = t_limit - t;
		double runDist = (mySpeed - enemySpeed) * t +
						 0.5 * (myAcc - enAcc) * t * t +
						 (myTop - enemySpeed - enAcc * t) * deltaT -
						 0.5 * enAcc * deltaT * deltaT;
		return runDist > distance;
	}
}

void player_ai(Info& info)
{
	vector<CellInfo> myCell;
	int maxCell = 0;// index in myCell

	TPlayerID myID = info.myID;
	for (int i = 0; i < info.cellInfo.size(); i++)
		if (info.cellInfo[i].ownerid == myID)
			myCell.push_back(info.cellInfo[i]);

	for (int i = 0; i < myCell.size(); i++)
		if (myCell[i].r > myCell[maxCell].r)
			maxCell = i;
	if (myCell.empty()) return;
	int nearestEnemy = 0;// index in all cells
	while (nearestEnemy < info.cellInfo.size() && info.cellInfo[nearestEnemy].ownerid == myID) nearestEnemy++;
	for (int k = nearestEnemy + 1; k < info.cellInfo.size(); k++) {
		if (info.cellInfo[k].ownerid == myID) continue;
		CellInfo e = info.cellInfo[k];
		CellInfo m = myCell[maxCell];
		if (distCell(e, m, true) < distCell(info.cellInfo[nearestEnemy], m, true))
			nearestEnemy = k;
	}
	if (nearestEnemy >= info.cellInfo.size()) return;


	for (int i = 0; i < myCell.size(); i++)
	{
		int split = splitCheck(myCell, maxCell, i, info.round, info.cellInfo[nearestEnemy]);
		double targetX = 10000, targetY = 10000;
		if (split != -1) {
			targetX = myCell[split].x;
			targetY = myCell[split].y;
		}
		else {
			for (auto& k : info.nutrientInfo) {
				double t = 1 - sqrt(2) / 3;
				if (min(abs(k.nux), abs(N - k.nux)) <= myCell[i].r * t) continue;
				if (min(abs(k.nuy), abs(N - k.nuy)) <= myCell[i].r * t) continue;
				if ((targetX - myCell[i].x) * (targetX - myCell[i].x) + (targetY - myCell[i].y) * (targetY - myCell[i].y) >
					(k.nux - myCell[i].x) * (k.nux - myCell[i].x) + (k.nuy - myCell[i].y) * (k.nuy - myCell[i].y)) {
					if (safe(info, myCell[i].x, myCell[i].y, myCell[i].r, k.nux, k.nuy)) {
						targetX = k.nux;
						targetY = k.nuy;
					}
				}
			}

			for (auto& k : info.cellInfo) {
				if (k.ownerid == myID) continue;
				double t = 1 - sqrt(2) / 3;
				if (min(abs(k.x), abs(N - k.x)) <= myCell[i].r * t) continue;
				if (min(abs(k.y), abs(N - k.y)) <= myCell[i].r * t) continue;
				if (k.r < myCell[i].r * lam && (targetX - myCell[i].x) * (targetX - myCell[i].x) + (targetY - myCell[i].y) * (targetY - myCell[i].y) >
					(k.x - myCell[i].x) * (k.x - myCell[i].x) + (k.y - myCell[i].y) * (k.y - myCell[i].y)) {
					if (safe(info, myCell[i].x, myCell[i].y, myCell[i].r, k.x, k.y)) {
						targetX = k.x;
						targetY = k.y;
					}
				}
			}
		}
		int direction = 0;
		double pi = 3.14159265;
		//cout << "targetX:" << targetX << " targetY:" << targetY << endl;
		//cout << "myCellX:" << myCell[i].x << " myCellY:" << myCell[i].y << endl;
		if (info.round > 800) {
			direction = compute_dir(150, 150, myCell[i].x, myCell[i].y);
			info.myCommandList.addCommand(Move, myCell[i].id, direction);
		}
		else {
			if (targetX < N + 1)
			{
				direction = compute_dir(targetX, targetY, myCell[i].x, myCell[i].y);
				info.myCommandList.addCommand(Move, myCell[i].id, direction);
			}
			else
			{
				direction = compute_dir(150, 150, myCell[i].x, myCell[i].y);
				info.myCommandList.addCommand(Move, myCell[i].id, direction);
			}
		}

	}

}
