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

typedef double TPosition;//����ֵ
typedef double TRadius;//�뾶

class Player;
class BattleField;
class Cell;
class Nutrient;
class Game;
class FireNet;
class Wind;

const int RESISTANCE = 4;//��û�и�ϸ���κβ�����ʱ��ϸ�����ԣ�RESISTANCE/R���ļ��ٶȼ��٣�RΪ��ǰϸ���뾶

const int NUTRIENTMOVEROUND = 3;//����Ӫ�����ʻغ�
const int MAXCELLNUM = 6;//��ҵ����ϸ����
const int PERTIME = 1;

const int BFSIZE = 300;//��ͼ�Ĵ�С
const int NURMAX = 5;//Ӫ�����ʰ뾶����
const int RCELLMIN = 6;//ϸ����С�뾶
const int minR = 6;//RCELLMIN����һ��
const int RViewR = 4;//����Ҫ

const int SPIKYNUM = 2;
const int SBallNUM = 10;//��ͼ�д�������

const int FirenetID = -1;//������ID

//�ĸ�����쫷��ID
const int Wind1ID = -21;
const int Wind2ID = -22;
const int Wind3ID = -23;
const int Wind4ID = -24;

//����쫷紦��������ϸ��ÿ�غ������С�ٷֱ�
const int shrinkpercentage = 4;
//쫷��ƶ��ٶ�
const int WindVelocity = 10;
//쫷�뾶
const int WindRadius = 30;
//쫷��ƶ��غ���
const int WINDROUNDMAX = 50;
//ÿ���������ϸ���޵лغ���
const int NOHARMROUNDMAX = 10;
//��������غ���
const int RAINROUNDMAX = 50;
//���������������ͼ�����ӵ�Ӫ��������Ŀ
const int RAINNUNUM = 30;

//�������ͣ��ƶ������ѡ�����
enum CommandType
{
	Move = 1,
	Division,
	Spit
};
//result.txt������ͣ���Ĳ������������ı�
enum FoutType
{
	NewCircle = 1, //������ϸ��������1 ��� ���� �뾶 λ��x λ��y��
	CircleDie, //�Ѿ����ڵ�ϸ������������2 ��� ����0 �뾶0 λ��0 λ��0��
	CircleChange //��ı��С���߰뾶������3 ��� �°뾶 ��λ��x ��λ��y��
};
//�����Ϣ
struct PlayerInfo
{
	TPlayerID id;//���ID

	bool live;//�Ƿ���ȫ����
	int hasRevive;//�Ƿ��Ѿ������

	TRank rank;          //��ѡ������/��1��ʼ
	TCredit score;//����

	TCellNum cell_num;//ϸ��������
	int kill_num;//��ɱ��Ŀ
	int survival_round;//���غ���
};

//ϸ����Ϣ
struct CellInfo {
	TCellID id;//ϸ��ID
	TPlayerID ownerid;//������ҵ�ID

	bool live;//�Ƿ�����
	//ϸ������
	double x;
	double y;
	//ϸ���뾶
	double r;
	//ϸ���ƶ��ٶ�
	double v;
	//ϸ���ƶ��Ƕ�
	int d;

};

//������Ϣ
struct FirenetInfo {
	//��������
	TPosition ex;
	TPosition ey;
	//�����뾶
	TRadius er;
	//����ID
	TFireNetID eID;
};

//쫷���Ϣ
struct WindInfo {
	//쫷�����
	TPosition wx;
	TPosition wy;
	//쫷�뾶
	TRadius wr;
	//쫷�ID
	TWindID wID;
};

//Ӫ��������Ϣ
struct NutrientInfo {
	//Ӫ����������
	TPosition nux;
	TPosition nuy;
	//Ӫ�����ʰ뾶
	TRadius nur;
	//Ӫ������ID
	TNutrientID nuID;
	//Ӫ�������ƶ��ٶ�
	double nuv;
	//Ӫ�������ƶ��Ƕ�
	int nud;
};

//����Ӫ��������Ϣ
struct SpikyBallInfo {
	//����Ӫ����������
	TPosition sx;
	TPosition sy;
	//����Ӫ�����ʰ뾶
	TRadius sr;
	//����Ӫ������ID
	TSpikyBallID sID;
};

//������Ϣ
struct Command
{
	Command(CommandType _type, vector<int> _parameters) :
		type(_type), parameters(_parameters) {}
	Command() {}
	CommandType type; //��������
	vector<int> parameters;  //����
};

class CommandList
{
public:
	//�������
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

	//�Ƴ�����
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
			throw std::out_of_range("��������ʱԽ��");
		return m_commands[n];
	}
public:
	vector<Command> m_commands;
};

//��Ϣ
struct Info
{
	int playerSize;  //�������
	int playerAlive; //ʣ�������
	TPlayerID myID;  //ѡ��ID��

	TRound round;     //�غ���
	CommandList myCommandList;  //�����б�
	PlayerInfo playerInfo;   //�����Ϣ
	vector<CellInfo> cellInfo;  //ϸ����Ϣ
	vector<NutrientInfo> nutrientInfo;  //Ӫ��������Ϣ
	vector<FirenetInfo> firenetInfo;  //������Ϣ
	vector<WindInfo> windInfo;  //쫷���Ϣ
	vector<SpikyBallInfo> spikyballInfo; //����Ӫ��������Ϣ
	//=�ع�
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
