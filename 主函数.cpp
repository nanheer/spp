#include"����.h"
#include"stdafx.h"
#include<iostream>
#include<string>
#include"stdlib.h"
int _tmain(int argc, _TCHAR* argv[])
{
	//string strn = "E:\\���㶨λ����\\����\\brdm3340.16p", stro = "E:\\���㶨λ����\\����\\harb3340.16o";
	//	string strn = "E:\\����\\���\\jfng1800.16n", stro = "E:\\����\\���\\jfng1800.16o";
	//string strn = "E:\\���㶨λ����\\����\\brdm2900.16p", stro = "E:\\���㶨λ����\\����\\harb2900.16o";
	string strn = "E:\\���㶨λ����\\����\\brdm2120.17p", stro = "E:\\���㶨λ����\\����\\twtf2120.17o";
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