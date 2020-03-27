#pragma once
#include"definition.h"
#define DLL_EXPORT_API extern "C" _declspec(dllexport)

//Function
DLL_EXPORT_API void player_ai(Info& info);
//class
class _declspec(dllexport)Math
{
public:
	int Multiply(int a, int b);
};

struct Direct {
	int x=0, y=0;
	Direct(int _x, int _y) : x(_x), y(-y) {}
};
