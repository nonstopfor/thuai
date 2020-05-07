#include <stdio.h>
#include "mydll.h"
#include "definition.h"
#include <math.h>
#include <time.h>
#include <algorithm>
#include <cstring>
#include <queue>
#include <map>
using namespace std;

#define MI 23
#define NI 19991
#define PI 3.1415926
#define LAM 0.9
#define N 300
#define MAX_SCORE 1e10
TPlayerID myID;
/*
info.myCommandList.addCommand(Division,aim_cell_id,direction);
info.myCommandList.addCommand(Move,aim_cell_id,direction);
info.myCommandList.addCommand(Spit,aim_cell_id,direction);
*/
double dist(double x1, double y1, double x2, double y2) {
	double deltaX = x1 - x2;
	double deltaY = y1 - y2;
	return sqrt(deltaX * deltaX + deltaY * deltaY);
}

double distCell(CellInfo& c1, CellInfo& c2, bool removeRadius = false) {
	double distRaw = dist(c1.x, c1.y, c2.x, c2.y);
	if (removeRadius) distRaw -= c1.r + c2.r;
	return distRaw;
}

double maxSpeed(CellInfo& cell) {
	return 20 / cell.r;
}
double accelerate(CellInfo& cell) {
	return 10 / cell.r;
}
double brakeLen(CellInfo& cell) {
	double maxSpd = maxSpeed(cell),
		acc = accelerate(cell);
	return maxSpd * maxSpd / 2 / acc;
}
double threatenR(CellInfo& myCell, CellInfo& cell) {
	//4帧的原因是2帧减速，1帧防止当前帧他直接吃我，1帧防止是他先动
	//直接用3/4半径用于防止他脸上刷营养，别的细胞上去送了这类的情况导致的半径突增
	double newEnemyR = cell.r / sqrt(2);
	bool canDivEatMe = newEnemyR * LAM > myCell.r;
	double divDist = canDivEatMe ? 1.2 * newEnemyR : 0;
	return 2.0 / 3 * cell.r + divDist + 4.0 * maxSpeed(cell);
}
int point_dir(double x1, double y1, double x2, double y2) {
	double dx = x2 - x1;
	double dy = y2 - y1;
	int dir = (int)(atan2(dy, dx) / PI * 180 + 360) % 360;
	return dir;
}
bool eat_nut(CellInfo& cell, NutrientInfo& nut) {
	if (nut.nur >= cell.r) return false;
	if (dist(cell.x, cell.y, nut.nux, nut.nuy) > 2 * cell.r / 3.0) return false;
	return true;
}
bool eat_cell(CellInfo& cell, CellInfo& target) {
	if (target.r / cell.r >= LAM) return false;
	if (dist(cell.x, cell.y, target.x, target.y) > 2 * cell.r / 3.0) return false;
	return true;
}


struct status;
vector<status>all_status;

struct status {
	//A*搜索过程中的一个状态
	double x;
	double y;
	double r;
	double v;
	int d;
	int step = 0;//搜索步数
	double score = 0;//相当于g
	int fa = -1;//父亲在vector中的下标
	int num = 0;//在vector中的下标
	double k = 0.3;//奖励衰减因子
	double h = 0;//估价值
	status(double _x, double _y, double _r, double _v, int _d, int _step = 0) :x(_x), y(_y), r(_r), v(_v), d(_d), step(_step) {

	}
	double get_ave_score() const {
		return score / step + h;
	}
	bool operator<(const status& t) const {
		return get_ave_score() < t.get_ave_score();
	}
	int get_root() {
		int cur_num = num;
		//cout << cur_num << endl;
		while (all_status[cur_num].fa != 0) {
			cur_num = all_status[cur_num].fa;
		}
		return cur_num;
	}
	status move(int dir) {
		status t = *this;
		const double maxSpeed = 20.0 / t.r;
		const double acc = 10.0 / t.r;
		t.v *= cos(PI * (dir - t.d) / 180.0);
		double deltaX = 0;
		if (t.v + acc > maxSpeed) {
			double accTime = (maxSpeed - t.v) / acc;
			deltaX = t.v * accTime + 0.5 * acc * accTime * accTime;
			deltaX += maxSpeed * (1.0 - accTime);
			t.v = maxSpeed;
		}
		else {
			deltaX = t.v + 0.5 * acc;
			t.v += acc;
		}
		t.x += deltaX * cos(dir * PI / 180.0);
		t.y += deltaX * sin(dir * PI / 180.0);
		t.d = dir;
		if (t.x < t.r) t.x = t.r;
		else if (t.x > N - t.r) t.x = N - t.r;
		if (t.y < t.r) t.y = t.r;
		else if (t.y > N - t.r) t.y = N - t.r;
		t.step++;
		t.fa = num;
		return t;
	}

	void update_score(vector<NutrientInfo>& nut_info, vector<CellInfo>& cell_info) {
		CellInfo cell;
		cell.x = x;
		cell.y = y;
		cell.r = r;
		//先对细胞进行判断，因为有可能发现不安全，直接return，可以省时间
		for (auto& enemy : cell_info) {
			if (enemy.ownerid == myID) continue;
			if (eat_cell(cell, enemy)) {
				score += (PI * enemy.r * enemy.r + 500) / exp(k * (step - 1));
				cell.r = sqrt(cell.r * cell.r + enemy.r * enemy.r);

			}
			else if (cell.r / enemy.r < LAM && dist(cell.x, cell.y, enemy.x, enemy.y) < threatenR(cell, enemy) + brakeLen(cell)) {
				score -= (PI * cell.r * cell.r + 500) / exp(k * (step - 1)) * 10; //+ 100000;
				//break;
				//score = -MAX_SCORE;
				//return;
			}
		}

		for (auto& nut : nut_info) {
			if (eat_nut(cell, nut)) {
				score += (PI * nut.nur * nut.nur) / exp(k * (step - 1));
				cell.r = sqrt(cell.r * cell.r + nut.nur * nut.nur);
			}
		}

		//估计h值
		for (auto& nut : nut_info) {
			if (nut.nur >= r) continue;
			int dir = point_dir(x, y, nut.nux, nut.nuy);
			if (abs(dir - d) < 90) {
				if (dist(x, y, nut.nux, nut.nuy) * sin(abs(dir - d) * PI / 180) < r * 2 / 3) {
					h += PI * nut.nur * nut.nur;
				}
			}
		}

		for (auto& enemy : cell_info) {
			if (enemy.ownerid == myID) continue;
			if (enemy.r / r >= LAM && r / enemy.r >= LAM) continue;
			int dir = point_dir(x, y, enemy.x, enemy.y);
			if (abs(dir - d) < 90) {
				if (enemy.r / r < LAM) {
					if (dist(x, y, enemy.x, enemy.y) * sin(abs(dir - d) * PI / 180) < r * 2 / 3) {
						h += PI * enemy.r * enemy.r + 500;
					}
				}
				if (r / enemy.r < LAM) {
					if (dist(x, y, enemy.x, enemy.y) * sin(abs(dir - d) * PI / 180) < enemy.r * 2 / 3) {
						h -= PI * r * r + 500;
					}
				}
			}
		}
	}
};

vector<int>get_dirs(status s0, Info& info) {
	vector<int>dirs;
	//暂时按10度一间隔
	for (int i = 0; i < 360; i += 5) {
		dirs.push_back(i);
	}
	return dirs;
}
int times = 0;
int get_best_move_dir(status s0, Info& info, double start_time) {
	//只考虑一定范围内的细胞和营养物质
	vector<NutrientInfo>newnutinfo;
	vector<CellInfo>newcellinfo;

	double minDistOfEligibleNut = 1e10,
		distLimit = 10 * s0.r;
	NutrientInfo mostCloseNut;//之后最好改成下标，如果之后对nut有更改的话
	for (auto& nut : info.nutrientInfo) {
		double distance = dist(nut.nux, nut.nuy, s0.x, s0.y);
		if (s0.r > nut.nur && distance < minDistOfEligibleNut) {
			minDistOfEligibleNut = distance;
			mostCloseNut = nut;
		}
		if (distance > distLimit) continue;
		newnutinfo.push_back(nut);
	}
	if (newnutinfo.size() == 0) newnutinfo.push_back(mostCloseNut);

	for (auto& cell : info.cellInfo) {
		if (dist(cell.x, cell.y, s0.x, s0.y > 10 * s0.r)) continue;
		newcellinfo.push_back(cell);
	}
	//cout << "newnutinfo.size():" << newnutinfo.size() << endl;
	//cout << "newcellinfo.size():" << newcellinfo.size() << endl;
	all_status.clear();
	all_status.push_back(s0);
	priority_queue<status>q;
	//queue<status>q;
	vector<int>dirs = get_dirs(s0, info);
	for (auto& dir : dirs) {
		status w = s0.move(dir);
		w.update_score(newnutinfo, newcellinfo);
		w.num = all_status.size();
		q.push(w);
		all_status.push_back(w);
	}
	double max_ave_score = -1;
	int best_status_num = 0;
	int max_step = 0;
	while (!q.empty()) {
		++times;
		double tmp_t = clock();
		status st = q.top(); q.pop();
		//cout << st.step << endl;
		//cout << st.get_root() << endl;
		max_step = max(max_step, st.step);
		double score = st.get_ave_score();
		if (score > max_ave_score) {
			max_ave_score = score;
			best_status_num = st.num;
		}
		double end_time = clock();
		if ((end_time - start_time) / CLOCKS_PER_SEC * 1000 > 140) break;
		if (st.step > 6) continue;//大于6步就不再搜了
		vector<int>tmp_dirs = get_dirs(st, info);
		for (auto& dir : tmp_dirs) {
			status w = st.move(dir);
			w.update_score(newnutinfo, newcellinfo);
			w.num = all_status.size();
			q.push(w);
			all_status.push_back(w);
		}
		//cout << q.size() << endl;
		//cout << "one status time:" << (clock() - tmp_t) / CLOCKS_PER_SEC * 1000 << endl;
		//cout << times << endl;
	}
	int best_dir = 0;
	while (best_status_num != 0) {
		best_dir = all_status[best_status_num].d;
		best_status_num = all_status[best_status_num].fa;
	}
	cout << "best_dir:" << best_dir << endl;
	cout << "max step:" << max_step << endl;
	return best_dir;
}

void player_ai(Info& info)
{
	double start_time = clock();
	cout << "round: " << info.round << " my score and rank: " << info.playerInfo.score << " " << info.playerInfo.rank << endl;

	vector<CellInfo> myCell;

	myID = info.myID;
	for (int i = 0; i < info.cellInfo.size(); i++)
		if (info.cellInfo[i].ownerid == myID)
			myCell.push_back(info.cellInfo[i]);

	if (myCell.empty()) return;
	times = 0;
	int cell_num = myCell.size();

	for (int cur = 0; cur < myCell.size(); cur++)
	{
		CellInfo& curCell = myCell[cur];
		//先判断能否分裂直接吃细胞
		if (cell_num < 6) {
			bool div_eat = false;
			int div_select = -1;
			double max_r = 0.0;
			for (int i = 0; i < info.cellInfo.size(); ++i) {
				auto& cell = info.cellInfo[i];
				if (cell.ownerid == myID) continue;
				if (cell.invincibleround) continue;
				if (cell.r / (curCell.r / (sqrt(2.0))) >= LAM) continue;
				double d = distCell(curCell, cell);
				if (d < (1.2 + (2.0 / 3)) * curCell.r / sqrt(2.0)) {
					if (cell.r > max_r) {
						max_r = cell.r;
						div_select = i;
						div_eat = true;
					}

				}
			}
			if (div_eat) {
				auto& cell = info.cellInfo[div_select];

				double dx = cell.x - curCell.x;
				double dy = cell.y - curCell.y;
				int direction = (int)(atan2(dy, dx) / PI * 180 + 360) % 360;
				info.myCommandList.addCommand(Division, curCell.id, direction);
				cell_num++;
				continue;

			}
		}
		status s0(curCell.x, curCell.y, curCell.r, curCell.v, curCell.d);
		int dir = get_best_move_dir(s0, info, start_time);
		info.myCommandList.addCommand(Move, curCell.id, dir);

		cout << "times: " << times << endl;
		double end_time = clock();
		cout << "time(ms): " << (end_time - start_time) / CLOCKS_PER_SEC * 1000 << endl;
	}

}