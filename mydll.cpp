#include <stdio.h>
#include "mydll.h"
#include "definition.h"
#include <time.h>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
using namespace std;

//#define DEBUG

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
info.myCommandList.addCommand(spit,aim_cell_id,direction);//吞吐命令，第二个参数是吞吐方向
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
double UNIONDIST = 60; // if too far away from ally big cells, don't union
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
	int target = -1;
	double minDist = MAXDIST;
	for (int i = 0; i < cells.size(); ++i) {
		double distance = distCell(cells[curCell], cells[i], true);
		if (i != curCell && distance < minDist) {
			target = i;
			minDist = distance;
		}
	}
	if (minDist < UNIONDIST)
		return target;
	else return -1;
}

bool Judge(PAIR P, PAIR yuan, double R) {
	//判断点是否在圆内
	if ((P.first - yuan.first) * (P.first - yuan.first) + (P.second - yuan.second) * (P.second - yuan.second) <= R * R) return true;
	return false;
}
bool Judis(PAIR P1, PAIR P2, PAIR yuan, double R) {
	if (Judge(P1, yuan, R) || Judge(P2, yuan, R)) return true;
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
		if (Judis(p1, p2, p3, r + 20 / r + 2 * ar / 3 + min(20 / cell.r, cell.v + 10 / cell.r))) return -2;
	}
	for (int i = 0; i < info.spikyballInfo.size(); ++i) {
		auto& t = info.spikyballInfo[i];
		double x = t.sx, y = t.sy;
		double ar = t.sr;
		if (ar >= r) continue;
		auto p3 = make_pair(x, y);
		if (Judis(p1, p2, p3, ar + r * 2 / 3)) return i;
	}
	return -1;
}
int compute_dir(double tx, double ty, double sx, double sy, double r = -1) {//算绕路就加r，是自己的半径

	double dx = tx - sx;
	double dy = ty - sy;
	double PI = 3.14159265;
	int spike = -1;
	for (int i = 0; i < globalInfo->spikyballInfo.size(); ++i) {
		if (r < 0) break;
		auto& t = globalInfo->spikyballInfo[i];
		double x = t.sx, y = t.sy;
		double ar = t.sr;
		if (ar > r) continue;
		auto p3 = make_pair(x, y);
		auto p1 = make_pair(sx, sy);
		auto p2 = make_pair(tx, ty);
		if (Judis(p1, p2, p3, ar + r * 1.1)) {// find nearest spike
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
	int direction = (int)(atan2(dy, dx) / PI * 180 + 360) % 360;
	if (spike != -1) {
		int direction2 = (int)(atan2(globalInfo->spikyballInfo[spike].sx - sx, globalInfo->spikyballInfo[spike].sy - sy) / PI * 180 + 360) % 360;
		if (direction2 < direction - 180) direction2 += 360;
		else if (direction2 > direction + 180) direction2 -= 360;
		if (direction2 < direction) direction = direction2 + 90;
		else direction = direction2 - 90;
		direction = (direction + 720) % 360;
	}
	//cout << "compute dir:" << direction << endl;
	return direction;
}

double INF = 1e10;
double LOOSEBOUND = 10; //如果距离减10刚好追上，也尝试去追

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
		double catchUpTime = 0;
		if (delta < 0 || (catchUpTime = -v_0 + sqrt(delta) < 0))
			return INF;
		return catchUpTime;
	}
	if (t >= t_limit) {
		double runDist = (mySpeed - enemySpeed) * t_limit +
			0.5 * (myAcc - enAcc) * t_limit * t_limit;
		return runDist - distance + LOOSEBOUND;
	}
	else {
		double deltaT = t_limit - t;
		double runDist = (mySpeed - enemySpeed) * t +
			0.5 * (myAcc - enAcc) * t * t +
			(myTop - enemySpeed - enAcc * t) * deltaT -
			0.5 * enAcc * deltaT * deltaT;
		return runDist - distance + LOOSEBOUND;
	}
}
double timeConsume(CellInfo me, CellInfo enemy) {
	return distAndTime(me, enemy, true);
}

double gain_nutrient(CellInfo& mycell, NutrientInfo& nut) {
	//吃营养物质的收益
	double d = max(1e-5, dist(mycell.x, mycell.y, nut.nux, nut.nuy) - mycell.r * 2 / 3);
	return PI * nut.nur * nut.nur / (d / 20 * mycell.r);
}

double gain_increase = 10;
double gain_cell(CellInfo& mycell, CellInfo& enemy, int num) {
	//double d = dist(mycell.x, mycell.y, enemy.x, enemy.y) - mycell.r * 2 / 3;
	double t = timeConsume(mycell, enemy);
	//cout << "timeConsume: " << t << endl;
	return (PI * enemy.r * enemy.r + 500) / timeConsume(mycell, enemy) + exp(num * num) * gain_increase;
}

bool catchable(CellInfo me, CellInfo enemy) {
	double reach = distAndTime(me, enemy);
	return reach > 0;
}

vector<int>getdangeridx(Info& info, vector<CellInfo>& myCell) {
	vector<int>res;
	TPlayerID myID = info.myID;
	for (int i = 0; i < myCell.size(); ++i) {
		int nearestEnemy = 0;// index in all cells
		while (nearestEnemy < info.cellInfo.size() && info.cellInfo[nearestEnemy].ownerid == myID) nearestEnemy++;
		for (int k = nearestEnemy + 1; k < info.cellInfo.size(); k++) {
			if (info.cellInfo[k].ownerid == myID) continue;
			CellInfo& e = info.cellInfo[k];
			if (distCell(e, myCell[i], true) < distCell(info.cellInfo[nearestEnemy], myCell[i], true))
				nearestEnemy = k;
		}
		if (nearestEnemy >= info.cellInfo.size()) continue;
		double r = myCell[i].r;
		double enemyr = info.cellInfo[nearestEnemy].r;
		double dist2enemy = dist(myCell[i].x, myCell[i].y, info.cellInfo[nearestEnemy].x, info.cellInfo[nearestEnemy].y);

		if (r / enemyr < 0.9 && dist2enemy <= DISTFACTOR * (r + enemyr)) {
			res.push_back(i);
		}
	}
	return res;
}

double safe_factor_round(int round) {
	if (round < 200) return 1.0;
	else if (round < 500) return 1.2;
	else if (round < 800) return 1.5;
	else return 10000;
}
bool safe_cell(CellInfo me, Info& info) {
	//判断分裂后是否安全
	TPlayerID myID = info.myID;
	me.r /= sqrt(2.0);
	me.v = 0;
	for (auto& cell : info.cellInfo) {
		if (cell.ownerid == myID) continue;
		if (me.r / cell.r > lam) continue;
		double d = distCell(me, cell, true);

		if (d < 0.85 * safe_factor_round(info.round) * me.r || (d < 1.3 * safe_factor_round(info.round) * me.r && catchable(cell, me))) {
			//cout << info.round << " d and me.r cell.r: " << d << ' ' << me.r << ' ' << cell.r << endl;
			return false;
		}

	}
	return true;
}
double get_danger_dist(CellInfo me, CellInfo enemy) {
	return distCell(me, enemy) - 2.5 * (20 / enemy.r) - 2 * enemy.r / 3;
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

#ifdef DEBUG
	vector<stringstream> debugInfo(myCell.size());
	for (int i = 0; i < debugInfo.size(); i++)
		debugInfo[i] << "[DEBUG] Round " << info.round << " Cell_sum: " << myCell.size() << " Cell " << i << " x=" << myCell[i].x << " y=" << myCell[i].y << " r=" << myCell[i].r << endl;
#endif

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
	vector<int>dangercell_idx = getdangeridx(info, myCell);
	vector<int>cell_clamp(info.cellInfo.size(), 0);
	//先计算出对方细胞被几个我方细胞包围
	for (int i = 0; i < info.cellInfo.size(); ++i) {
		auto& cell = info.cellInfo[i];
		if (cell.ownerid == myID) continue;
		for (int j = 0; j < info.cellInfo.size(); ++j) {
			auto& me = info.cellInfo[j];
			if (me.ownerid != myID) continue;
			if (cell.r / me.r >= lam) continue;
			if (distCell(me, cell) <= cell.r) cell_clamp[i]++;
		}
	}

	for (int cur = 0; cur < myCell.size(); cur++)
	{
		CellInfo& curCell = myCell[cur];
		int direction = 0;

		//先检查是否需要逃跑
		{
			int nearest = -1;//最近敌人id
			for (int k = 0; k < info.cellInfo.size(); k++) {
				if (info.cellInfo[k].ownerid == myID) continue;
				if (info.cellInfo[k].r * 0.9 <= curCell.r) continue;
				if (get_danger_dist(curCell, info.cellInfo[k]) > 0) continue;
				if (nearest == -1) nearest = k;
				else if (distCell(curCell, info.cellInfo[k]) < distCell(curCell, info.cellInfo[nearest]))
					nearest = k;
			}
			int nearest2 = -1;//次近敌人id
			for (int k = 0; k < info.cellInfo.size(); k++) {
				if (k == nearest) continue;
				if (info.cellInfo[k].ownerid == myID) continue;
				if (info.cellInfo[k].r * 0.9 <= curCell.r) continue;
				if (get_danger_dist(curCell, info.cellInfo[k]) > 0) continue;
				if (nearest2 == -1) nearest2 = k;
				else if (distCell(curCell, info.cellInfo[k]) < distCell(curCell, info.cellInfo[nearest2]))
					nearest2 = k;
			}
#ifdef DEBUG
			debugInfo[cur] << "\ttargetX >= N + 1, nearest = " << nearest << " nearest2 = " << nearest2 << endl;
#endif
			if (nearest != -1) {
				direction = compute_dir(curCell.x, curCell.y,
					info.cellInfo[nearest].x, info.cellInfo[nearest].y);
				if (nearest2 != -1) {
					int direction2 = compute_dir(curCell.x, curCell.y,
						info.cellInfo[nearest2].x, info.cellInfo[nearest2].y);
#ifdef DEBUG
					debugInfo[cur] << "\tdirection = " << direction << " direction2 = " << direction2 << endl;
#endif

					if (direction2 < direction - 180) direction2 += 360;
					else if (direction2 > direction + 180) direction2 -= 360;
					direction = ((direction + direction2) / 2 + 360) % 360;
				}
#ifdef DEBUG
				debugInfo[cur] << "\t\tRun Away, direction = " << direction << endl;
#endif
				//开始判断撞边
				double predictX = curCell.x + (maxSpeed(curCell) + 1.1 * curCell.r) * cos((double)direction / 360 * 2 * PI);
				double predictY = curCell.y + (maxSpeed(curCell) + 1.1 * curCell.r) * sin((double)direction / 360 * 2 * PI);
				if (predictX <= 0) {
					if (direction < 135) direction = 180 - direction;
					else if (direction > 225) direction = 540 - direction;
					else if (direction < 180) direction = 90;
					else direction = 270;
				}
				else if (predictX >= N) {
					if (direction > 315) direction = 270;
					else if (direction > 180) direction = 540 - direction;
					else if (direction > 45) direction = 180 - direction;
					else direction = 90;

				}
				else if (predictY <= 0) {
					if (direction < 225) direction = 360 - direction;
					else if (direction > 315) direction = 360 - direction;
					else if (direction < 270) direction = 180;
					else direction = 0;
				}
				else if (predictY >= N) {
					if (direction < 45) direction = 360 - direction;
					else if (direction > 135) direction = 360 - direction;
					else if (direction < 90) direction = 0;
					else direction = 180;
				}
#ifdef DEBUG
				debugInfo[cur] << "\t\tAfter Checking Boundary, direction = " << direction << endl;
#endif

				info.myCommandList.addCommand(Move, curCell.id, direction);
				continue;
			}
		}

		int split = splitCheck(myCell, maxCell, dangercell_idx, cur, info.round);
		//int split = splitCheck(myCell, maxCell, cur, info.round, info.cellInfo[nearestEnemy]);
		double targetX = 10000, targetY = 10000;
#ifdef DEBUG
		debugInfo[cur] << "\tSplit Check = " << split << endl;
#endif
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
				//if (distCell(curCell, info.cellInfo[j]) > 1.5 * info.cellInfo[j].r) continue;
				int w = safe(info, curCell.x, curCell.y, curCell.r, k.x, k.y);
				if (w != -2) {
					//如果不是路径上有其他细胞
					cell_idx.push_back(j);
				}
			}
			sort(cell_idx.begin(), cell_idx.end(), [&](int a, int b) {
				double g1 = gain_cell(curCell, info.cellInfo[a], cell_clamp[a]);
				double g2 = gain_cell(curCell, info.cellInfo[b], cell_clamp[b]);
				return g1 > g2;
				});


			if (myCell.size() < 6 && curCell.r > sqrt(2) * MINR && safe_cell(curCell, info) && cell_idx.size() + nutrient_idx.size() > 1) {
				double gain_1 = 0;//不分裂的最大收益
				double gain_2 = 0;//分裂的最大收益
				double tx1 = 0, ty1 = 0;//不分裂时目标位置
				double tx2 = 0, ty2 = 0;//分裂时冲向的目标位置
				double gain_1_cell = cell_idx.empty() ? -1 : gain_cell(curCell, info.cellInfo[cell_idx[0]], cell_clamp[cell_idx[0]]);
				double gain_1_nut = nutrient_idx.empty() ? -1 : gain_nutrient(curCell, info.nutrientInfo[nutrient_idx[0]]);
				if (gain_1_nut > gain_1_cell) {
					gain_1 = gain_1_nut;
					tx1 = info.nutrientInfo[nutrient_idx[0]].nux;
					ty1 = info.nutrientInfo[nutrient_idx[0]].nuy;
				}
				else {
					gain_1 = gain_1_cell;
					tx1 = info.cellInfo[cell_idx[0]].x;
					ty1 = info.cellInfo[cell_idx[0]].y;
				}
				struct p {
					double gain;
					double x;
					double y;
					p(double _gain, double _x, double _y) :gain(_gain), x(_x), y(_y) {}
					bool operator<(const p& t) {
						return gain > t.gain;
					}
				};
				vector<p>tmp;
				for (int k = 0; k < min((int)nutrient_idx.size(), 2); ++k) {
					auto& nut = info.nutrientInfo[nutrient_idx[k]];
					tmp.push_back(p(gain_nutrient(curCell, nut), nut.nux, nut.nuy));
				}
				int cnt = 0;
				for (int k = 0; k < cell_idx.size(); ++k) {
					auto& cell = info.cellInfo[cell_idx[k]];
					if (cell.r / (curCell.r / sqrt(2)) >= lam) continue;

					tmp.push_back(p(gain_cell(curCell, cell, cell_clamp[cell_idx[k]]), cell.x, cell.y));
					if (++cnt > 1) break;
				}
				sort(tmp.begin(), tmp.end());
				//cout << info.round << ": " << tmp.size() << endl;
				if (tmp.size() > 1) {
					gain_2 = tmp[0].gain + tmp[1].gain;
					if (gain_2 > gain_1) {
						int dir1 = compute_dir(tmp[0].x, tmp[0].y, curCell.x, curCell.y);
						info.myCommandList.addCommand(Division, curCell.id, dir1);
						continue;
					}
					else {
						if (gain_1_nut > gain_1_cell) vis[nutrient_idx[0]] = true;
						targetX = tx1;
						targetY = ty1;
					}
				}
				else {
					if (gain_1_nut > gain_1_cell) vis[nutrient_idx[0]] = true;
					targetX = tx1;
					targetY = ty1;
				}
			}
			else if (cell_idx.size() + nutrient_idx.size() > 0) {
				double max_gain_nut = nutrient_idx.empty() ? -1 : gain_nutrient(curCell, info.nutrientInfo[nutrient_idx[0]]);
				double max_gain_cell = cell_idx.empty() ? -1 : gain_cell(curCell, info.cellInfo[cell_idx[0]], cell_clamp[cell_idx[0]]);
				if (max_gain_nut > max_gain_cell) {
					vis[nutrient_idx[0]] = true;
					targetX = info.nutrientInfo[nutrient_idx[0]].nux;
					targetY = info.nutrientInfo[nutrient_idx[0]].nuy;
				}
				else {
					targetX = info.cellInfo[cell_idx[0]].x;
					targetY = info.cellInfo[cell_idx[0]].y;
				}
			}

			/*
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
			*/
#ifdef DEBUG
			debugInfo[cur] << "\t After nutrient found. targetX = " << targetX << " targetY = " << targetY << endl;
#endif

		}

		if (info.round > 800 && cur == maxCell) {
			int nearest = -1;//最近敌人id
			for (int k = 0; k < info.cellInfo.size(); k++) {
				if (info.cellInfo[k].ownerid == myID) continue;
				if (info.cellInfo[k].r * 0.9 <= curCell.r) continue;
				if (nearest == -1) nearest = k;
				else if (distCell(curCell, info.cellInfo[k]) < distCell(curCell, info.cellInfo[nearest]))
					nearest = k;
			}
			if (nearest != -1 && distCell(curCell, info.cellInfo[nearest]) < info.cellInfo[nearest].r) {
				direction = compute_dir(curCell.x, curCell.y,
					info.cellInfo[nearest].x, info.cellInfo[nearest].y);
				info.myCommandList.addCommand(Move, curCell.id, direction);
			}
			else {
				direction = compute_dir(150, 150, curCell.x, curCell.y);
				info.myCommandList.addCommand(Move, curCell.id, direction);
			}
#ifdef DEBUG
			debugInfo[cur] << "\tinfo.round > 800 && cur == maxCell, nearest = " << nearest << "direction = " << direction << endl;
#endif
		}
		else {
			if (targetX < N + 1)
			{
				direction = compute_dir(targetX, targetY, curCell.x, curCell.y, curCell.r);
				info.myCommandList.addCommand(Move, curCell.id, direction);
#ifdef DEBUG
				debugInfo[cur] << "\ttargetX < N + 1, direction = " << direction << endl;
#endif
			}
			else
			{
				// check if enemy too near
				int nearest = -1;//最近敌人id
				for (int k = 0; k < info.cellInfo.size(); k++) {
					if (info.cellInfo[k].ownerid == myID) continue;
					if (info.cellInfo[k].r * 0.9 <= curCell.r) continue;
					if (nearest == -1) nearest = k;
					else if (distCell(curCell, info.cellInfo[k]) < distCell(curCell, info.cellInfo[nearest]))
						nearest = k;
				}
				int nearest2 = -1;//次近敌人id
				for (int k = 0; k < info.cellInfo.size(); k++) {
					if (k == nearest) continue;
					if (info.cellInfo[k].ownerid == myID) continue;
					if (info.cellInfo[k].r * 0.9 <= curCell.r) continue;
					if (nearest2 == -1) nearest2 = k;
					else if (distCell(curCell, info.cellInfo[k]) < distCell(curCell, info.cellInfo[nearest2]))
						nearest2 = k;
				}
#ifdef DEBUG
				debugInfo[cur] << "\ttargetX >= N + 1, nearest = " << nearest << " nearest2 = " << nearest2 << endl;
#endif
				if (nearest != -1 && distCell(curCell, info.cellInfo[nearest], true) < 1.0 * curCell.r) {
					direction = compute_dir(curCell.x, curCell.y,
						info.cellInfo[nearest].x, info.cellInfo[nearest].y);
					if (nearest2 != -1 && distCell(curCell, info.cellInfo[nearest], true) < 1.3 * curCell.r) {
						int direction2 = compute_dir(curCell.x, curCell.y,
							info.cellInfo[nearest2].x, info.cellInfo[nearest2].y);
#ifdef DEBUG
						debugInfo[cur] << "\tdirection = " << direction << " direction2 = " << direction2 << endl;
#endif

						if (direction2 < direction - 180) direction2 += 360;
						else if (direction2 > direction + 180) direction2 -= 360;
						direction = ((direction + direction2) / 2 + 360) % 360;
					}
#ifdef DEBUG
					debugInfo[cur] << "\t\tRun Away, direction = " << direction << endl;
#endif
					//开始判断撞边
					double predictX = curCell.x + (maxSpeed(curCell) + 1.1 * curCell.r) * cos((double)direction / 360 * 2 * PI);
					double predictY = curCell.y + (maxSpeed(curCell) + 1.1 * curCell.r) * sin((double)direction / 360 * 2 * PI);
					if (predictX <= 0) {
						if (direction < 135) direction = 180 - direction;
						else if (direction > 225) direction = 540 - direction;
						else if (direction < 180) direction = 90;
						else direction = 270;
					}
					else if (predictX >= N) {
						if (direction > 315) direction = 270;
						else if (direction > 180) direction = 540 - direction;
						else if (direction > 45) direction = 180 - direction;
						else direction = 90;

					}
					else if (predictY <= 0) {
						if (direction < 225) direction = 360 - direction;
						else if (direction > 315) direction = 360 - direction;
						else if (direction < 270) direction = 180;
						else direction = 0;
					}
					else if (predictY >= N) {
						if (direction < 45) direction = 360 - direction;
						else if (direction > 135) direction = 360 - direction;
						else if (direction < 90) direction = 0;
						else direction = 180;
					}
#ifdef DEBUG
					debugInfo[cur] << "\t\tAfter Checking Boundary, direction = " << direction << endl;
#endif

					info.myCommandList.addCommand(Move, curCell.id, direction);
				}
				else {
#ifdef DEBUG
					debugInfo[cur] << "\t\tFind Safe, direction = " << direction;
#endif
					bool flag = false;
					for (double angle = 0; angle < 360; angle += 1) {
						double dx = cos(angle / 360 * 2 * PI) * N;
						double dy = sin(angle / 360 * 2 * PI) * N;
						if (safe(info, curCell.x, curCell.y, curCell.r, curCell.x + dx, curCell.y + dy) == -1) {
							direction = compute_dir(curCell.x + dx, curCell.y + dy, curCell.x, curCell.y, curCell.r);
							info.myCommandList.addCommand(Move, curCell.id, direction);
							flag = true;
							break;
						}
					}
#ifdef DEBUG
					debugInfo[cur] << " After: direction = " << direction;
					if (!flag) debugInfo[cur] << "    not move" << endl;
					else debugInfo[cur] << endl;
#endif

				}
			}
		}
#ifdef DEBUG
		cout << debugInfo[cur].str();
#endif
	}

}
