#include"函数.h"
#include"stdafx.h"
#include<iostream>
#include<string>
#include"stdlib.h"
int _tmain(int argc, _TCHAR* argv[])
{
	//string strn = "E:\\单点定位程序\\数据\\brdm3340.16p", stro = "E:\\单点定位程序\\数据\\harb3340.16o";
	//	string strn = "E:\\桌面\\结果\\jfng1800.16n", stro = "E:\\桌面\\结果\\jfng1800.16o";
	//string strn = "E:\\单点定位程序\\数据\\brdm2900.16p", stro = "E:\\单点定位程序\\数据\\harb2900.16o";
	string strn = "E:\\单点定位程序\\数据\\brdm2120.17p", stro = "E:\\单点定位程序\\数据\\twtf2120.17o";
	pof po = new of;
	pnf pn = new nf;
	presult pt = new result;
	readofile(stro, po);
	readnfile(strn, pn);
	spp(po, pn, pt);
	putresult(pt);
	delete po;
	delete pn;
	delete pt;
	return 0;
}