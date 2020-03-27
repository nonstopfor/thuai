#include <stdio.h>
#include "mydll.h"
#include "definition.h"
#include <math.h>
#include <time.h>
#include <vector>
#define MI 23
#define NI 19991
/*
基本的三个添加指令的命令
info.myCommandList.addCommand(Division,aim_cell_id,direction);//分裂命令，第二个参数是分裂方向
info.myCommandList.addCommand(Move,aim_cell_id,direction);//移动命令，第二个参数是移动方向
info.myCommandList.addCommand(Spit,aim_cell_id,direction);//吞吐命令，第二个参数是吞吐方向
*/

using namespace std;

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
		int split = splitCheck(myCell,maxCell,i,info.round)
		int targetX=10000,targetY=10000;
		if(split!=-1){
			targetX = myCell[split].x;
			targetY = myCell[split].y;
		}else{
			for(auto k:info.nutrientInfo){
				if((targetX-myCell[i].x)*(targetX-myCell[i].x)+(targetY-myCell[i].y)*(targetY-myCell[i].y)<
					(k.x-myCell[i].x)*(k.x-myCell[i].x)+(k.y-myCell[i].y)*(k.y-myCell[i].y)){
					if(safe(info,myCell[i].x,myCell[i].y,myCell[i].r,k.x,k.y)){
						targetX = k.x;
						targetY = k.y;
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
