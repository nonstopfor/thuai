#include <stdio.h>
#include "mydll.h"
#include "definition.h"
#include <math.h>
#include <time.h>
#define MI 23
#define NI 19991
/*
�������������ָ�������
info.myCommandList.addCommand(Division,aim_cell_id,direction);//��������ڶ��������Ƿ��ѷ���
info.myCommandList.addCommand(Move,aim_cell_id,direction);//�ƶ�����ڶ����������ƶ�����
info.myCommandList.addCommand(Spit,aim_cell_id,direction);//��������ڶ������������·���
*/

void player_ai(Info& info)
{
	
	int Mycellnum = -1;
	int Mycell_x, Mycell_y = -1;
	int MincellNum = -1;
	int MincellRadius = info.cellInfo[0].r;
	int Mincell_x, Mincell_y = -1;
	TPlayerID myid = info.myID;
	for (int i = 0; i < info.cellInfo.size(); i++)
	{
		if (info.cellInfo[i].ownerid == myid)
		{
			Mycellnum = info.cellInfo[i].id;
			Mycell_x = info.cellInfo[i].x;
			Mycell_y = info.cellInfo[i].y;
		}
		if (info.cellInfo[i].r <= MincellRadius && Mycellnum != i)
		{
			Mincell_x = info.cellInfo[i].x;
			Mincell_y = info.cellInfo[i].y;
			MincellNum = i;
		}
	}
	if (Mycellnum != -1)
	{
		int direction = 0;
		double pi = 3.14159265;
		int dx, dy = 0;
		if (MincellNum != -1)
		{
			dx = Mincell_x - Mycell_x;
			dy = Mincell_y - Mycell_y;
			direction = (int)(atan2(dy, dx) / pi * 180 + 360) % 360;
			info.myCommandList.addCommand(Move, Mycellnum, direction);
		}
		else
		{
			dx = 150 - Mycell_x;
			dy = 150 - Mycell_y;
			direction = (int)(atan2(dy, dx) / pi * 180 + 360) % 360;
			info.myCommandList.addCommand(Move, Mycellnum, direction);
		}
	}

}
