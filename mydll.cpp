#include <stdio.h>
#include "mydll.h"
#include "definition.h"
#include <math.h>
#include <time.h>
#include <algorithm>
#include <cstring>
#include <queue>
#include <map>
#include <assert.h>
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
	return 20.0 / cell.r;
}
double accelerate(CellInfo& cell) {
	return 10.0 / cell.r;
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
double gain_nut(CellInfo& me, NutrientInfo& nut) {
	double u = max(0.1, dist(me.x, me.y, nut.nux, nut.nuy) - 2 * me.r / 3);
	double gain = PI * nut.nur * nut.nur * maxSpeed(me) / u;
	return gain;
}
double gain_cell_eat(CellInfo& me, CellInfo& enemy) {
	//吃对方
	double u = max(0.1, dist(me.x, me.y, enemy.x, enemy.y) - 2 * me.r / 3);
	double gain = (PI * enemy.r * enemy.r + 500) * maxSpeed(me) / u;
	return gain;
}
double gain_cell_run(CellInfo& me, CellInfo& enemy) {
	//逃跑
	double u = dist(me.x, me.y, enemy.x, enemy.y);
	double gain = 0;
	if (u < threatenR(me, enemy) + brakeLen(me)) {
		gain = (PI * me.r * me.r + 500) * 10000 * (maxSpeed(me) + maxSpeed(enemy)) / u;
	}
	//double gain = (PI * me.r * me.r + 500) * 2 * 5 * 10 * 10 * (maxSpeed(me) + maxSpeed(enemy)) / u;
	return -gain;
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
double factor[8] = { 0,1,0.9,0.81,0.729,0.6561,0.59049,0.531441 };
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
	double k = 1.0;//奖励衰减因子
	double h = 0;//估价值
	bool end = false;//为true的时候不再拓展
	bool safe = true;//是否安全
	double ave_score = 0;

	status(double _x, double _y, double _r, double _v, int _d, int _step = 0) :x(_x), y(_y), r(_r), v(_v), d(_d), step(_step) {

	}
	double get_ave_score() const {
		return score + h;
		//return (score + h) / step;
		if (score + h > 0) return (score + h) * factor[step];
		else return (score + h) / factor[step];
	}
	double get_dist_score(status s0) {
		double dis = dist(x, y, s0.x, s0.y);
		return (score + h) / dis;
	}
	double get_sum_score() {
		return score + h;
	}
	bool operator<(const status& t) const {
		return ave_score < t.ave_score;
	}

	int get_root(vector<status>& v) {
		int cur_num = num;
		while (v[cur_num].fa != 0) {
			cur_num = v[cur_num].fa;
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
		t.h = 0;//g继承父节点的，h另外算
		return t;
	}

	void checkSafe(vector<CellInfo>& cell_info) {
		CellInfo cell;
		cell.x = x;
		cell.y = y;
		cell.r = r;
		cell.v = v;
		cell.d = d;
		safe = true;
		for (auto& enemy : cell_info) {
			if (enemy.ownerid == myID) continue;
			if (cell.r / enemy.r < LAM && gain_cell_run(cell, enemy) < 0) {
				safe = false;
				return;
			}
		}
	}

	void update_score(vector<NutrientInfo>& nut_info, vector<CellInfo>& cell_info, vector<status>& all_status, Info& info) {
		CellInfo cell;
		cell.x = x;
		cell.y = y;
		cell.r = r;
		cell.v = v;
		cell.d = d;
		//checkSafe(cell_info);
		//先对细胞进行判断，因为有可能发现不安全，直接return，可以省时间
		for (auto& enemy : cell_info) {
			if (enemy.ownerid == myID) continue;
			if (eat_cell(cell, enemy)) {
				if (info.round > 300 || step == 1) {
					score += (PI * enemy.r * enemy.r + 500) / step;
					cell.r = sqrt(cell.r * cell.r + enemy.r * enemy.r);
					end = true;
				}
				
			}else if (eat_cell(enemy, cell)) {
                score = -MAX_SCORE;
                end = true;
                return;
            }

			for (auto& nut : nut_info) {
				if (eat_nut(cell, nut)) {
					score += (PI * nut.nur * nut.nur) / step;
					cell.r = sqrt(cell.r * cell.r + nut.nur * nut.nur);
					end = true;
				}
			}

			//估计h值
			int count = 0;
			double threatenh = 0;
			for (auto& nut : nut_info) {
				if (nut.nur >= r) continue;
				h += gain_nut(cell, nut);
				++count;
			}

			for (auto& enemy : cell_info) {
				if (enemy.ownerid == myID) continue;
				if (enemy.r / r >= LAM && r / enemy.r >= LAM) continue;
				if (enemy.r / r < LAM) {
					if (info.round > 300) {
						h += gain_cell_eat(cell, enemy);
						++count;
					}
				}
				else if (r / enemy.r < LAM) {
					//威胁算在g里面
					//h += gain_cell_run(cell, enemy);
					double gain = gain_cell_run(cell, enemy);
					if (gain < 0) safe = false;
					threatenh += gain;
					//score += gain;
					//++count;
				}
			}
			if (count) h /= count;
			h += threatenh;

			if (fa != -1 && all_status[fa].safe && !safe) {//进入危险
				score -= 10000;
			}
			if (fa != -1 && !all_status[fa].safe && safe) {//离开危险
				score += 7000; //进入过危险
			}
			ave_score = get_ave_score();
		}
	};

};
vector<int>get_dirs(status s0, status st, Info& info) {
	vector<int>dirs;
	//离当前方向较近的更密

	int d = st.step >= 1 ? point_dir(s0.x, s0.y, st.x, st.y) : s0.d;
	for (int i = 0; i < 60; i += 5) {
		dirs.push_back((i + d) % 360);
		dirs.push_back((d - i + 360) % 360);
	}
	if (st.step < 1) {
		for (int i = 60; i < 180; i += 10) {
			dirs.push_back((i + d) % 360);
			dirs.push_back((d - i + 360) % 360);
		}
	}
	else {
		for (int i = 60; i < 90; i += 10) {
			dirs.push_back((i + d) % 360);
			dirs.push_back((d - i + 360) % 360);
		}
	}
	return dirs;
}
int times = 0;
void clear(priority_queue<status>& q) {
	//priority_queue<status>().swap(q);
	while (!q.empty()) q.pop();
}
void clear(vector<status>& v) {
	//vector<status>().swap(v);
	while (!v.empty()) v.clear();
}

int get_best_move_dir(status s0, Info& info, double start_time, double max_time) {
	//只考虑一定范围内的细胞和营养物质
	vector<NutrientInfo>newnutinfo;
	vector<CellInfo>newcellinfo;

	double minDistOfEligibleNut = 1e10,
		distLimit = 10 * s0.r,
		t = 1 - sqrt(2) / 3;
	NutrientInfo mostCloseNut;//之后最好改成下标，如果之后对nut有更改的话
	for (auto& nut : info.nutrientInfo) {
		double distance = dist(nut.nux, nut.nuy, s0.x, s0.y);
		if (s0.r > nut.nur && distance < minDistOfEligibleNut) {
			minDistOfEligibleNut = distance;
			mostCloseNut = nut;
		}
		if (distance > distLimit ||
			fmin(fabs(nut.nux), fabs(BFSIZE - nut.nux)) <= s0.r * t ||//地图卡边
			fmin(fabs(nut.nuy), fabs(BFSIZE - nut.nuy)) <= s0.r * t)
			continue;
		newnutinfo.push_back(nut);
	}
	if (newnutinfo.size() == 0) newnutinfo.push_back(mostCloseNut);

	for (auto& cell : info.cellInfo) {
		if (dist(cell.x, cell.y, s0.x, s0.y) > 10 * s0.r) continue;
		double t = 1 - sqrt(2) / 3;
		if (min(abs(cell.x), abs(N - cell.x)) <= s0.r * t) continue;
		if (min(abs(cell.y), abs(N - cell.y)) <= s0.r * t) continue;
		newcellinfo.push_back(cell);
	}


	int size = 1;
	priority_queue<status>q;
	vector<status>all_status;

	s0.checkSafe(newcellinfo);
	//s0.update_score(newnutinfo, newcellinfo, all_status);
	all_status.push_back(s0);

	vector<int>dirs = get_dirs(s0, s0, info);
	for (auto& dir : dirs) {
		status w = s0.move(dir);
		w.update_score(newnutinfo, newcellinfo, all_status, info);
		w.num = size;

		q.push(w);
		all_status.push_back(w);
		++size;
	}
	double max_score = -1e10;
	int best_status_num = 0;
	int max_step = 0;
	//map<int, int>cnt;
	while (!q.empty()) {
		//++times;
		status st = q.top(); q.pop();
		//cout << "step/score/h/ave_score:" << st.step << "/" << st.score << "/" << st.h << "/" << st.get_ave_score() << endl;

		//cnt[st.get_root(all_status)]++;
		//max_step = max(max_step, st.step);
		//double score = st.get_dist_score(s0);
		//double score = st.ave_score;
		double score = st.get_sum_score();
		if (score > max_score) {
			max_score = score;
			best_status_num = st.num;
		}
		double end_time = clock();
		if ((end_time - start_time) / CLOCKS_PER_SEC * 1000 > max_time) break;
		if (st.end) continue;//不再拓展该节点
		if (st.step > 6) continue;//大于6步就不再搜了
		vector<int>tmp_dirs = get_dirs(s0, st, info);
		for (auto& dir : tmp_dirs) {
			status w = st.move(dir);
			w.update_score(newnutinfo, newcellinfo, all_status, info);
			w.num = size;

			q.push(w);
			all_status.push_back(w);
			++size;
		}

	}
	int best_dir = 0;
	/*
	for (int i = 1; i <= 48; ++i) {
		cout << i << ":" << cnt[i] << endl;
		//cout << "score/h/ave_score:" << all_status[i].score << "/" << all_status[i].h << "/" << all_status[i].get_ave_score() << endl;
	}
	*/
	//cout << "before while:best_status_num/all_status.size():" << best_status_num << "/" << all_status.size() << endl;
	while (best_status_num != 0) {
		//cout << all_status[best_status_num].score << "/" << all_status[best_status_num].h << endl;
		best_dir = all_status[best_status_num].d;
		best_status_num = all_status[best_status_num].fa;

	}
	//cout << "best_dir:" << best_dir << endl;
	//cout << "max step:" << max_step << endl;
	return best_dir;
}

bool div_eat(Info& info, CellInfo& myCell, CellInfo& enemy) {
	bool div_eat = false;
	if (enemy.invincibleround) return false;
	if (enemy.r / (myCell.r / (sqrt(2.0))) >= LAM) return false;
	double d = distCell(myCell, enemy);
	if (d < (1.2 + (2.0 / 3)) * myCell.r / sqrt(2.0))
		div_eat = true;
	return div_eat;
}

bool div_safe(Info& info, CellInfo& me, double tx, double ty,
	double add_r) {
	//add_r是被吃目标的半径
	//判断向(tx,ty)位置分裂是否安全
	double dx = tx - me.x;
	double dy = ty - me.y;
	TPlayerID myID = info.myID;

	//接着判断是否分裂后仍然安全
	CellInfo stay, rush;    //分裂出的两个细胞
	stay.x = me.x; stay.y = me.y; stay.r = me.r / sqrt(2);
	double rushRatio = 1.2 * stay.r / sqrt(dx * dx + dy * dy);//1.2是rush距离比新半径的倍数
	rush.x = dx * rushRatio + stay.x;
	rush.y = dy * rushRatio + stay.y;
	rush.r = sqrt(stay.r * stay.r + add_r * add_r);
	bool bothAreSafe = true;
	for (int i = 0; i < info.cellInfo.size(); ++i) {
		if (info.cellInfo[i].ownerid == myID) continue;
		auto& enemy = info.cellInfo[i];
		if (dist(stay.x, stay.y, enemy.x, enemy.y) < threatenR(stay, enemy) ||
			dist(rush.x, rush.y, enemy.x, enemy.y) < threatenR(rush, enemy)) {
			bothAreSafe = false;
			break;
		}
	}
	return bothAreSafe;
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

	int cell_num = myCell.size();

	for (int cur = 0; cur < myCell.size(); cur++)
	{
		times = 0;
		CellInfo& curCell = myCell[cur];
		//先判断能否分裂直接吃细胞
		//注意最大细胞数
		if (cell_num < MAXCELLNUM && curCell.r > sqrt(2) * minR) {
			bool div = false;
			int direction = -1;
			int tarIdx = -1;
			double maxEatR = 0;
			for (int i = 0; i < info.cellInfo.size(); ++i) {
				auto& cell = info.cellInfo[i];
				if (cell.ownerid == info.myID) continue;
				if (div_eat(info, curCell, cell) && cell.r > maxEatR &&
					div_safe(info, curCell, cell.x, cell.y, cell.r)) {
					div = true;
					tarIdx = i;
					maxEatR = cell.r;
                    cout<<"Eat\n";
				}
			}
			if (div) {
				auto& cell = info.cellInfo[tarIdx];

				double dx = cell.x - curCell.x;
				double dy = cell.y - curCell.y;
				int direction = (int)(atan2(dy, dx) / PI * 180 + 360) % 360;
				info.myCommandList.addCommand(Division, curCell.id, direction);
				cell_num++;
				continue;
			}
		}
		/*
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
		*/

		status s0(curCell.x, curCell.y, curCell.r, curCell.v, curCell.d);
		double tmp_start_time = clock();
		double max_time = 140.0 / cell_num;
		int dir = get_best_move_dir(s0, info, tmp_start_time, max_time);
		info.myCommandList.addCommand(Move, curCell.id, dir);

		//cout << "times: " << times << endl;
		double end_time = clock();
		//cout << "time(ms): " << (end_time - start_time) / CLOCKS_PER_SEC * 1000 << endl;
	}

}
