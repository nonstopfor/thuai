#ifndef DEFINITION_H
#define DEFINITION_H
#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>
using std::map;
using std::vector;
using std::string;

typedef int TPlayerID;
typedef int TCellID;
typedef int TNutrientID;
typedef int TSpikyBallID;
typedef int TID;
typedef int TWindID;
typedef int TFireNetID;

typedef double TCredit;
typedef int TRank;
typedef int TRound;
typedef int TCellNum;

typedef double TPosition;//坐标值
typedef double TRadius;//半径

class Player;
class BattleField;
class Cell;
class Nutrient;
class Game;
class FireNet;
class Wind;

const int RESISTANCE = 4;//当没有给细胞任何操作的时候，细胞会以（RESISTANCE/R）的加速度减速，R为当前细胞半径

const int NUTRIENTMOVEROUND = 3;//喷吐营养物质回合
const int MAXCELLNUM = 6;//玩家的最大细胞数
const int PERTIME = 1;

const int BFSIZE = 300;//地图的大小
const int NURMAX = 5;//营养物质半径上限
const int RCELLMIN = 6;//细胞最小半径
const int minR = 6;//RCELLMIN保持一致
const int RViewR = 4;//不需要

const int SPIKYNUM = 2;
const int SBallNUM = 10;//地图中刺球数量

const int FirenetID = -1;//火网的ID

//四个角落飓风的ID
const int Wind1ID = -21;
const int Wind2ID = -22;
const int Wind3ID = -23;
const int Wind4ID = -24;

//处于飓风处或火网外的细胞每回合面积缩小百分比
const int shrinkpercentage = 4;
//飓风移动速度
const int WindVelocity = 10;
//飓风半径
const int WindRadius = 30;
//飓风移动回合数
const int WINDROUNDMAX = 50;
//每次重生后的细胞无敌回合数
const int NOHARMROUNDMAX = 10;
//下雨持续回合数
const int RAINROUNDMAX = 50;
//下雨过程中整个地图所增加的营养物质数目
const int RAINNUNUM = 30;

//命令类型：移动、分裂、吞吐
enum CommandType
{
	Move = 1,
	Division,
	Spit
};
//result.txt输出类型：球的产生、死亡、改变
enum FoutType
{
	NewCircle = 1, //产生新细胞（类型1 编号 从属 半径 位置x 位置y）
	CircleDie, //已经存在的细胞死亡（类型2 编号 从属0 半径0 位置0 位置0）
	CircleChange //求改变大小或者半径（类型3 编号 新半径 新位置x 新位置y）
};
//玩家信息
struct PlayerInfo
{
	TPlayerID id;//玩家ID

	bool live;//是否完全死亡
	int hasRevive;//是否已经复活过

	TRank rank;          //该选手排名/从1开始
	TCredit score;//积分

	TCellNum cell_num;//细胞的数量
	int kill_num;//击杀数目
	int survival_round;//存活回合数
};

//细胞信息
struct CellInfo {
	TCellID id;//细胞ID
	TPlayerID ownerid;//所属玩家的ID

	bool live;//是否死亡
	//细胞坐标
	double x;
	double y;
	//细胞半径
	double r;
	//细胞移动速度
	double v;
	//细胞移动角度
	int d;

};

//火网信息
struct FirenetInfo {
	//火网坐标
	TPosition ex;
	TPosition ey;
	//火网半径
	TRadius er;
	//火网ID
	TFireNetID eID;
};

//飓风信息
struct WindInfo {
	//飓风坐标
	TPosition wx;
	TPosition wy;
	//飓风半径
	TRadius wr;
	//飓风ID
	TWindID wID;
};

//营养物质信息
struct NutrientInfo {
	//营养物质坐标
	TPosition nux;
	TPosition nuy;
	//营养物质半径
	TRadius nur;
	//营养物质ID
	TNutrientID nuID;
	//营养物质移动速度
	double nuv;
	//营养物质移动角度
	int nud;
};

//刺型营养物质信息
struct SpikyBallInfo {
	//刺型营养物质坐标
	TPosition sx;
	TPosition sy;
	//刺型营养物质半径
	TRadius sr;
	//刺型营养物质ID
	TSpikyBallID sID;
};

//命令信息
struct Command
{
	Command(CommandType _type, vector<int> _parameters) :
		type(_type), parameters(_parameters) {}
	Command() {}
	CommandType type; //命令类型
	vector<int> parameters;  //参数
};

class CommandList
{
public:
	//添加命令
	void addCommand(CommandType _type, vector<int> _parameters)
	{
		Command c;
		c.type = _type;
		c.parameters = _parameters;
		m_commands.push_back(c);
	}
	void addCommand(CommandType _type, int para1)
	{
		Command c;
		c.type = _type;
		c.parameters.push_back(para1);
		m_commands.push_back(c);
	}
	void addCommand(CommandType _type, int para1, int para2)
	{
		Command c;
		c.type = _type;
		c.parameters.push_back(para1);
		c.parameters.push_back(para2);
		m_commands.push_back(c);
	}
	void addCommand(CommandType _type, int para1, int para2, int para3)
	{
		Command c;
		c.type = _type;
		c.parameters.push_back(para1);
		c.parameters.push_back(para2);
		c.parameters.push_back(para3);
		m_commands.push_back(c);
	}

	//移除命令
	void removeCommand(int n)
	{
		m_commands.erase(m_commands.begin() + n);
	}
	int  size() const { return int(m_commands.size()); }
	vector<Command>::iterator begin() { return m_commands.begin(); }
	vector<Command>::iterator end() { return m_commands.end(); }
	vector<Command>::const_iterator  begin()const { return m_commands.cbegin(); }
	vector<Command>::const_iterator end() const { return m_commands.cend(); }
	Command& operator[](int n)
	{
		if (n < 0 || size() <= n)
			throw std::out_of_range("访问命令时越界");
		return m_commands[n];
	}
public:
	vector<Command> m_commands;
};

//信息
struct Info
{
	int playerSize;  //总玩家数
	int playerAlive; //剩余玩家数
	TPlayerID myID;  //选手ID号

	TRound round;     //回合数
	CommandList myCommandList;  //命令列表
	PlayerInfo playerInfo;   //玩家信息
	vector<CellInfo> cellInfo;  //细胞信息
	vector<NutrientInfo> nutrientInfo;  //营养物质信息
	vector<FirenetInfo> firenetInfo;  //火网信息
	vector<WindInfo> windInfo;  //飓风信息
	vector<SpikyBallInfo> spikyballInfo; //刺型营养物质信息
	//=重构
	Info& operator=(Info& tempinfo) {
		this->myID = tempinfo.myID;
		this->playerAlive = tempinfo.playerAlive;
		this->playerInfo = tempinfo.playerInfo;
		this->round = tempinfo.round;
		this->myCommandList = tempinfo.myCommandList;
		this->cellInfo = tempinfo.cellInfo;
		this->firenetInfo = tempinfo.firenetInfo;
		this->windInfo = tempinfo.windInfo;
		this->nutrientInfo = tempinfo.nutrientInfo;
		return *this;
	}
};

#endif
