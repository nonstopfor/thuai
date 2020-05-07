#include <stdio.h>
#include "mydll.h"
#include "definition.h"
#include <math.h>
#include <time.h>
#include <algorithm>
#include <cstring>
#define MI 23
#define NI 19991
#define PI 3.1415926
#define LAM 0.9
#define N 300
/*
ָ
info.myCommandList.addCommand(Division,aim_cell_id,direction);
info.myCommandList.addCommand(Move,aim_cell_id,direction);
info.myCommandList.addCommand(Spit,aim_cell_id,direction);
*/

struct status {
	//A*搜索过程中的一个状态
	double x;
	double y;
	double r;
	double v;
	int d;
	int step = 0;//搜索步数
	double score = 0;
	status(double _x, double _y, double _r, double _v, int _d, int _step) :x(_x), y(_y), r(_r), v(_v), step(_step) {

	}
	status move(int dir) {
		status t = *this;

	}
};

void player_ai(Info& info)
{
	vector<CellInfo> myCell;

	TPlayerID myID = info.myID;
	for (int i = 0; i < info.cellInfo.size(); i++)
		if (info.cellInfo[i].ownerid == myID)
			myCell.push_back(info.cellInfo[i]);

	if (myCell.empty()) return;
	for (int cur = 0; cur < myCell.size(); cur++)
	{
		CellInfo& curCell = myCell[cur];

	}
}

