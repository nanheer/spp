#ifndef FUNCTION_H
#define FUNCTION_H
#include<string>
#include"结构体.h"
#include"常量.h"
//#include"matrix.h"
using namespace std;
//读n文件
void readnfile(string strn, pnf pnavefile);

//读o文件
void readofile(string stro,pof obsfile);

//GPST到UTC
void gpsttoutc(pgpst gt, ptc ut);

//UTC到GPST
void utctogpst(ptc ut,pgpst gt);

//UTC到JULIANDAY
void utctojulianday(ptc ut,pjulian ju);

//JULIANDAY到UTC
void juliandaytoutc(pjulian ju,ptc ut);

//GPST到JULIANDAY
void gpsttojulianday(pgpst gt, pjulian ju);

//JULIANDAY到GPST
void juliandaytogpst(pjulian ju, pgpst gt);
//求儒略日差值
double deltjulianday(pjulian ju1, pjulian ju2);
//儒略日变换
void transjulian(pjulian ju1, double * dt, pjulian ju2);
//

//XYZ到BLH
void xyztoblh(pxyz px,pblh pb);
//BLH到XYZ
void blhtoxyz(pblh pb, pxyz px);
//XYZ到ENU
void xyztoenu(pxyz pxcenter, pxyz px, penu pe);
//ENU到XZY
void enutoxyz(pxyz pxcenter,penu pe, pxyz px);
//ENU到ENUPOLAR
void enutoenupolar(penu pe, penupolar pep);
//对流层延迟改正
void tropo(double mjd,pxyz px1, pxyz px2,double * dx);
//求气象数据
void get_nominal_metedata(double mjd, double lat, double lon,double dhgt, double* pres, double* temp, double* rhumi, double* undu);
//求较小值
int min(int a, int b);
//接收机天线高改正
void antena_height_correction(pxyz px1, pxyz px2, penu pe, double*p1, double*p2, double*pp1, double*pp2);
//寻找最优星历
bool find_best_ephem(string prn, ptc ut, pnf nfile, pse ephemeris);
//计算卫星坐标
void cal_sate_coor(string prn, ptc ut, pnf pnavefile, bool &flag, pxyz coor, double* dt, double* reletivity);

//迭代卫星坐标
void iterate_sate_coor(pxyz px, ptc ut, pse ephemeris, pxyz coor);

//相对论和卫星钟差改正
void sateclock(pse ephemeris, ptc ut, double clockoff, double delts);

//由平近点角计算偏近点角
double calE(double* M, double* e);
//由偏近点角求真近点角
double calf(double* E, double* e);
//旋转矩阵
//void rotamatrix(string str, double seta,Matrix r);
//输出结果
void putresult(presult pt);
//spp函数
void spp(pof obsfile, pnf pnavefile, presult pt);
#endif