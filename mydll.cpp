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
const double MINR = 6;
const double lam = 0.9;//小细胞与大细胞半径比值小于这个时被吞噬
double PI = 3.14159265;
Info* globalInfo;
typedef pair<double, double> PAIR;
/*
基本的三个添加指令的命令
info.myCommandList.addCommand(Division,aim_cell_id,direction);//分裂命令，第二个参数是分裂方向
info.myCommandList.addCommand(Move,aim_cell_id,direction);//移动命令，第二个参数是移动方向
info.myCommandList.addCommand(Spit,aim_cell_id,direction);//吞吐命令，第二个参数是吞吐方向
*/


double maxSpeed(CellInfo& cell) {
	return 20 / cell.r;
}

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

double MINBOUND = 0.1;
double MAXDIST = 425;
double DISTFACTOR = 1.1;
int DISASTERROUND = 700;
int HELP_RANGE = 50;
int splitCheck(std::vector<CellInfo>& cells, int maxCell, std::vector<int>& cellsIndanger, int curCell,
	int round) {

	double maxR = cells[maxCell].r;
	//double enemyR = nearestEnemy.r;
	bool minboundJudge = cells[curCell].r > MINBOUND * maxR;
	bool roundJudge = round < DISASTERROUND || curCell == maxCell;
	bool enemyJudge = cellsIndanger.size() == 0;//maxR > 0.9 * enemyR;
	if (!enemyJudge) {
		enemyJudge = true;
		for (int i = 0; i < cellsIndanger.size(); ++i) {
			if (curCell == cellsIndanger[i]) {
				enemyJudge = true;
				break;
			}
			CellInfo& tarCell = cells[cellsIndanger[i]];
			double circleDist = distCell(cells[curCell], tarCell);
			if (circleDist > HELP_RANGE)
				continue;
			else {
				//found target which needed help
				enemyJudge = false;
			}
		}
		/*double distToEnemy = dist(
			cells[curCell].x, cells[curCell].y,
			nearestEnemy.x, nearestEnemy.y
		);
		enemyJudge = distToEnemy > DISTFACTOR * (maxR + enemyR); //max cell is still far from enemy*/
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
int safe(Info& info, double x1, double y1, double r, double x2, double y2) {
	//-2表示路径上有其他细胞
	//-1表示路径安全
	//0-x 表示路径上有刺球，返回刺球idx
	TPlayerID myID = info.myID;
	auto p1 = make_pair(x1, y1);
	auto p2 = make_pair(x2, y2);

	for (int i = 0; i < info.cellInfo.size(); ++i) {
		auto& cell = info.cellInfo[i];
		if (cell.ownerid == myID) continue;
		if (r / cell.r > lam) continue;
		double x = cell.x, y = cell.y;
		double ar = cell.r;
		auto p3 = make_pair(x, y);
		if (Judis(p1, p2, p3, r + ar)) return -2;
	}
	for (int i = 0; i < info.spikyballInfo.size(); ++i) {
		auto& t = info.spikyballInfo[i];
		double x = t.sx, y = t.sy;
		double ar = t.sr;
		if (ar >= r) continue;
		auto p3 = make_pair(x, y);
		if (Judis(p1, p2, p3, 1.1 * ar + r)) return i;
	}
	return -1;
}
int compute_dir(double tx, double ty, double sx, double sy, double r = -1) {//算绕路就加r，是自己的半径
	double dx = tx - sx;
	double dy = ty - sy;
	double pi = 3.14159265;
	int spike = -1;
	for (int i = 0; i < globalInfo->spikyballInfo.size(); ++i) {
		if (r == -1) break;
		auto& t = globalInfo->spikyballInfo[i];
		double x = t.sx, y = t.sy;
		double ar = t.sr;
		auto p3 = make_pair(x, y);
		auto p1 = make_pair(sx, sy);
		auto p2 = make_pair(tx, ty);
		if (Judis(p1, p2, p3, ar + r)) {// find nearest spike
			if (spike == -1)
				spike = i;
			else {
				if (dist(sx, sy, globalInfo->spikyballInfo[spike].sx, globalInfo->spikyballInfo[spike].sy)
					- globalInfo->spikyballInfo[spike].sr >
					dist(sx, sy, x, y) - ar)
					spike = i;
			}
		}
	}
	int direction = (int)(atan2(dy, dx) / pi * 180 + 360) % 360;
	if (spike != -1) {
		direction += 90;
		direction %= 360;
	}
	return direction;
}

double distAndTime(CellInfo me, CellInfo enemy, bool time = false) {
	double distance = dist(me.x, me.y, enemy.x, enemy.y);
	distance = distance - 2.0 / 3.0 * me.r;
	double dist_hat_dir = compute_dir(enemy.x, enemy.y, me.x, me.y);
	double mySpeed = me.v * cos(me.d - dist_hat_dir);
	double enemySpeed = enemy.v * cos(enemy.v - dist_hat_dir);
	double myAcc = 10 / me.r, enAcc = 10 / enemy.r,
		myTop = 20 / me.r;
	double t = (myTop - mySpeed) / myAcc,
		t_limit = (myTop - enemySpeed) / enAcc;
	if (time) {
		double v_0 = mySpeed - enemySpeed, a = myAcc - enAcc;
		double delta = v_0 * v_0 + 2 * a * distance;
		if (delta < 0) return -1;
		return (-v_0 + sqrt(delta)) / a;
	}
	if (t >= t_limit) {
		double runDist = (mySpeed - enemySpeed) * t_limit +
			0.5 * (myAcc - enAcc) * t_limit * t_limit;
		return runDist - distance;
	}
	else {
		double deltaT = t_limit - t;
		double runDist = (mySpeed - enemySpeed) * t +
			0.5 * (myAcc - enAcc) * t * t +
			(myTop - enemySpeed - enAcc * t) * deltaT -
			0.5 * enAcc * deltaT * deltaT;
		return runDist - distance;
	}
}
double timeConsume(CellInfo me, CellInfo enemy) {
	return distAndTime(me, enemy, true);
}

double gain_nutrient(CellInfo& mycell, NutrientInfo& nut) {
	//吃营养物质的收益
	double d = dist(mycell.x, mycell.y, nut.nux, nut.nuy) - mycell.r * 2 / 3;
	return PI * nut.nur * nut.nur / (d / 20 * mycell.r);
}
double gain_cell(CellInfo& mycell, CellInfo& enemy) {
	//double d = dist(mycell.x, mycell.y, enemy.x, enemy.y) - mycell.r * 2 / 3;
	return (PI * enemy.r * enemy.r + 500) / timeConsume(mycell, enemy);
}

bool catchable(CellInfo me, CellInfo enemy) {
	double reach = distAndTime(me, enemy);
	return reach > 0;
}

vector<int>getdangeridx(Info& info) {
	vector<int>res;
	TPlayerID myID = info.myID;
	for (int i = 0; i < info.cellInfo.size(); ++i) {
		if (info.cellInfo[i].ownerid != myID) continue;
		int nearestEnemy = 0;// index in all cells
		while (nearestEnemy < info.cellInfo.size() && info.cellInfo[nearestEnemy].ownerid == myID) nearestEnemy++;
		for (int k = nearestEnemy + 1; k < info.cellInfo.size(); k++) {
			if (info.cellInfo[k].ownerid == myID) continue;
			CellInfo& e = info.cellInfo[k];
			if (distCell(e, info.cellInfo[i], true) < distCell(info.cellInfo[nearestEnemy], info.cellInfo[i], true))
				nearestEnemy = k;
		}
		if (nearestEnemy >= info.cellInfo.size()) continue;
		double r = info.cellInfo[i].r;
		double enemyr = info.cellInfo[nearestEnemy].r;
		double dist2enemy = dist(info.cellInfo[i].x, info.cellInfo[i].y, info.cellInfo[nearestEnemy].x, info.cellInfo[nearestEnemy].y);

		if (r / enemyr < 0.9 && dist2enemy <= DISTFACTOR * (r + enemyr)) {
			res.push_back(i);
		}
	}
	return res;
}



void player_ai(Info& info)
{
	globalInfo = &info;

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

	vector<bool>vis(info.nutrientInfo.size(), false);
	vector<int>dangercell_idx = getdangeridx(info);

	for (int cur = 0; cur < myCell.size(); cur++)
	{
		CellInfo& curCell = myCell[cur];
		int split = splitCheck(myCell, maxCell, dangercell_idx, cur, info.round);
		//int split = splitCheck(myCell, maxCell, cur, info.round, info.cellInfo[nearestEnemy]);
		double targetX = 10000, targetY = 10000;
		if (split != -1) {
			targetX = myCell[split].x;
			targetY = myCell[split].y;
		}
		else {
			vector<int>nutrient_idx;//将营养物质按平均收益大小排序

			for (int j = 0; j < info.nutrientInfo.size(); ++j) {
				if (vis[j]) continue;
				auto& k = info.nutrientInfo[j];
				if (k.nur >= curCell.r) continue;
				double t = 1 - sqrt(2) / 3;
				if (min(abs(k.nux), abs(N - k.nux)) <= curCell.r * t) continue;
				if (min(abs(k.nuy), abs(N - k.nuy)) <= curCell.r * t) continue;
				int w = safe(info, curCell.x, curCell.y, curCell.r, k.nux, k.nuy);
				if (w != -2) {
					//如果不是路径上有其他细胞
					nutrient_idx.push_back(j);
				}
			}
			sort(nutrient_idx.begin(), nutrient_idx.end(), [&](int a, int b) {
				double g1 = gain_nutrient(curCell, info.nutrientInfo[a]);
				double g2 = gain_nutrient(curCell, info.nutrientInfo[b]);
				return g1 > g2;
				/*
				double d1 = dist(curCell.x, curCell.y, info.nutrientInfo[a].nux, info.nutrientInfo[a].nuy);
				double d2 = dist(curCell.x, curCell.y, info.nutrientInfo[b].nux, info.nutrientInfo[b].nuy);
				return d1 < d2;
				*/
				});
			vector<int>cell_idx;//将其他小细胞按平均收益排序
			for (int j = 0; j < info.cellInfo.size(); ++j) {
				auto& k = info.cellInfo[j];
				if (k.ownerid == myID) continue;
				if (k.r >= curCell.r * lam) continue;
				double t = 1 - sqrt(2) / 3;
				if (min(abs(k.x), abs(N - k.x)) <= curCell.r * t) continue;
				if (min(abs(k.y), abs(N - k.y)) <= curCell.r * t) continue;
				if (!catchable(curCell, info.cellInfo[j])) continue;
				int w = safe(info, curCell.x, curCell.y, curCell.r, k.x, k.y);
				if (w != -2) {
					//如果不是路径上有其他细胞
					cell_idx.push_back(j);
				}
			}
			sort(cell_idx.begin(), cell_idx.end(), [&](int a, int b) {
				double g1 = gain_cell(curCell, info.cellInfo[a]);
				double g2 = gain_cell(curCell, info.cellInfo[b]);
				return g1 > g2;
				/*
				double d1 = dist(curCell.x, curCell.y, info.cellInfo[a].x, info.cellInfo[a].y);
				double d2 = dist(curCell.x, curCell.y, info.cellInfo[b].x, info.cellInfo[b].y);
				return d1 < d2;
				*/
				});
			double gmax = -1;
			if (cell_idx.empty() && nutrient_idx.empty()) {

			}
			else if (cell_idx.empty()) {
				if (nutrient_idx.size() >= 2 && myCell.size() < 6 && curCell.r > sqrt(2) * MINR && info.round < 150) {
					int dir1 = compute_dir(info.nutrientInfo[nutrient_idx[1]].nux, info.nutrientInfo[nutrient_idx[1]].nuy, curCell.x, curCell.y);
					info.myCommandList.addCommand(Division, curCell.id, dir1);
					continue;
				}
				else {
					vis[nutrient_idx[0]] = true;

					targetX = info.nutrientInfo[nutrient_idx[0]].nux;
					targetY = info.nutrientInfo[nutrient_idx[0]].nuy;

				}
			}
			else if (nutrient_idx.empty()) {
				targetX = info.cellInfo[cell_idx[0]].x;
				targetY = info.cellInfo[cell_idx[0]].y;
			}
			else {
				if (gain_cell(curCell, info.cellInfo[cell_idx[0]]) > gain_nutrient(curCell, info.nutrientInfo[nutrient_idx[0]])) {
					targetX = info.cellInfo[cell_idx[0]].x;
					targetY = info.cellInfo[cell_idx[0]].y;
				}
				else {
					vis[nutrient_idx[0]] = true;
					targetX = info.nutrientInfo[nutrient_idx[0]].nux;
					targetY = info.nutrientInfo[nutrient_idx[0]].nuy;

				}

			}

			/*
			double dn0 = 10000, dn1 = 10000, dc0 = 10000, dc1 = 10000;
			if (nutrient_idx.size() >= 2) dn1 = dist(curCell.x, curCell.y, info.nutrientInfo[nutrient_idx[1]].nux, info.nutrientInfo[nutrient_idx[1]].nuy);
			if (nutrient_idx.size() >= 1) dn0 = dist(curCell.x, curCell.y, info.nutrientInfo[nutrient_idx[0]].nux, info.nutrientInfo[nutrient_idx[0]].nuy);
			int i0 = 0;
			while (i0 < cell_idx.size() && !catchable(curCell, info.cellInfo[cell_idx[i0]])) ++i0;
			int i1 = i0 + 1;
			while (i1 < cell_idx.size() && !catchable(curCell, info.cellInfo[cell_idx[i1]])) ++i1;

			if (i1 < cell_idx.size()) dc1 = dist(curCell.x, curCell.y, info.cellInfo[cell_idx[i1]].x, info.cellInfo[cell_idx[i1]].y);
			if (i0 < cell_idx.size()) dc0 = dist(curCell.x, curCell.y, info.cellInfo[cell_idx[i0]].x, info.cellInfo[cell_idx[i0]].y);
			double dd = 20;//距离差的阙值
			if (dn1 < dc0) {
				//最近的两个是营养物质
				//cout << "两个营养物质" << endl;
				if (dn1 - dn0 < dd && myCell.size() < 6 && curCell.r > sqrt(2) * MINR && info.round < 150) {
					//int dir0 = compute_dir(curCell.x, curCell.y, info.nutrientInfo[nutrient_idx[0]].nux, info.nutrientInfo[nutrient_idx[0]].nuy);
					int dir1 = compute_dir(curCell.x, curCell.y, info.nutrientInfo[nutrient_idx[1]].nux, info.nutrientInfo[nutrient_idx[1]].nuy);
					info.myCommandList.addCommand(Division, curCell.id, dir1);
					//cout << "分裂" << endl;
					continue;
				}
				else {
					//cout << "追近的营养物质" << endl;
					if (dc0 < 10000 && info.round>200) {
						targetX = info.cellInfo[cell_idx[i0]].x;
						targetY = info.cellInfo[cell_idx[i0]].y;
					}
					else {
						vis[nutrient_idx[0]] = true;
						targetX = info.nutrientInfo[nutrient_idx[0]].nux;
						targetY = info.nutrientInfo[nutrient_idx[0]].nuy;
					}

				}
			}
			else if (dc1 < dn0) {
				//最近的两个是细胞
				targetX = info.cellInfo[cell_idx[i0]].x;
				targetY = info.cellInfo[cell_idx[i0]].y;
			}
			else {
				//最近的两个，一个是细胞，一个是营养物质
				//由于细胞可追，故优先追细胞
				if (i0 < cell_idx.size()) {
					targetX = info.cellInfo[cell_idx[i0]].x;
					targetY = info.cellInfo[cell_idx[i0]].y;
				}

			}
			*/
		}
		int direction = 0;
		double pi = 3.14159265;
		//cout << "targetX:" << targetX << " targetY:" << targetY << endl;
		//cout << "myCellX:" << curCell.x << " myCellY:" << curCell.y << endl;
		if (info.round > 800 && cur == maxCell) {
			int nearest = -1;//最近敌人id
			for (int k = 0; k < info.cellInfo.size(); k++) {
				if (info.cellInfo[k].ownerid == myID) continue;
				if (info.cellInfo[k].r * 0.9 <= curCell.r) continue;
				if (nearest == -1) nearest = k;
				else if (distCell(curCell, info.cellInfo[k]) < distCell(myCell[nearest], info.cellInfo[k]))
					nearest = k;
			}
			if (nearest != -1 && distCell(curCell, info.cellInfo[nearest], true) < 0.5 * curCell.r) {
				direction = compute_dir(curCell.x, curCell.y,
					info.cellInfo[nearest].x, info.cellInfo[nearest].y);
				info.myCommandList.addCommand(Move, curCell.id, direction);
			}
			else {
				direction = compute_dir(150, 150, curCell.x, curCell.y);
				info.myCommandList.addCommand(Move, curCell.id, direction);
			}

		}
		else {
			if (targetX < N + 1)
			{
				direction = compute_dir(targetX, targetY, curCell.x, curCell.y, curCell.r);
				info.myCommandList.addCommand(Move, curCell.id, direction);
			}
			else
			{
				// check if enemy too near
				int nearest = -1;//最近敌人id
				for (int k = 0; k < info.cellInfo.size(); k++) {
					if (info.cellInfo[k].ownerid == myID) continue;
					if (info.cellInfo[k].r * 0.9 <= curCell.r) continue;
					if (nearest == -1) nearest = k;
					else if (distCell(curCell, info.cellInfo[k]) < distCell(myCell[nearest], info.cellInfo[k]))
						nearest = k;
				}
				if (nearest != -1 && distCell(curCell, info.cellInfo[nearest], true) < 0.5 * curCell.r) {
					direction = compute_dir(curCell.x, curCell.y,
						info.cellInfo[nearest].x, info.cellInfo[nearest].y);

					//开始判断撞边
					double predictX = curCell.x + (maxSpeed(curCell) + curCell.r) * cos(direction / 360 * 2 * pi);
					double predictY = curCell.y + (maxSpeed(curCell) + curCell.r) * sin(direction / 360 * 2 * pi);
					if (predictX <= 0) {
						if (direction < 180) direction = 180 - direction;
						else if (direction > 180) direction = 540 - direction;
						else direction = 90;
					}
					else if (predictX >= N) {
						if (direction < 180) direction = 180 - direction;
						else if (direction > 180) direction = 540 - direction;
						else direction = 90;
					}
					else if (predictY <= 0) {
						if (direction < 270) direction = 360 - direction;
						else if (direction > 270) direction = 360 - direction;
						else direction = 0;
					}
					else if (predictY >= N) {
						if (direction < 90) direction = 360 - direction;
						else if (direction > 90) direction = 360 - direction;
						else direction = 0;
					}

					info.myCommandList.addCommand(Move, curCell.id, direction);
				}
				else {
					for (double angle = 0; angle < 360; angle += 1) {
						double dx = cos(angle / 360 * 2 * pi) * N;
						double dy = sin(angle / 360 * 2 * pi) * N;
						if (safe(info, curCell.x, curCell.y, curCell.r, curCell.x + dx, curCell.y + dy) == -1) {
							direction = compute_dir(curCell.x + dx, curCell.y + dy, curCell.x, curCell.y, curCell.r);
							info.myCommandList.addCommand(Move, curCell.id, direction);
							break;
						}
					}
				}
			}
		}

	}

}
