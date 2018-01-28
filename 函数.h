#ifndef FUNCTION_H
#define FUNCTION_H
#include<string>
#include"�ṹ��.h"
#include"����.h"
//#include"matrix.h"
using namespace std;
//��n�ļ�
void readnfile(string strn, pnf pnavefile);

//��o�ļ�
void readofile(string stro,pof obsfile);

//GPST��UTC
void gpsttoutc(pgpst gt, ptc ut);

//UTC��GPST
void utctogpst(ptc ut,pgpst gt);

//UTC��JULIANDAY
void utctojulianday(ptc ut,pjulian ju);

//JULIANDAY��UTC
void juliandaytoutc(pjulian ju,ptc ut);

//GPST��JULIANDAY
void gpsttojulianday(pgpst gt, pjulian ju);

//JULIANDAY��GPST
void juliandaytogpst(pjulian ju, pgpst gt);
//�������ղ�ֵ
double deltjulianday(pjulian ju1, pjulian ju2);
//�����ձ任
void transjulian(pjulian ju1, double * dt, pjulian ju2);
//

//XYZ��BLH
void xyztoblh(pxyz px,pblh pb);
//BLH��XYZ
void blhtoxyz(pblh pb, pxyz px);
//XYZ��ENU
void xyztoenu(pxyz pxcenter, pxyz px, penu pe);
//ENU��XZY
void enutoxyz(pxyz pxcenter,penu pe, pxyz px);
//ENU��ENUPOLAR
void enutoenupolar(penu pe, penupolar pep);
//�������ӳٸ���
void tropo(double mjd,pxyz px1, pxyz px2,double * dx);
//����������
void get_nominal_metedata(double mjd, double lat, double lon,double dhgt, double* pres, double* temp, double* rhumi, double* undu);
//���Сֵ
int min(int a, int b);
//���ջ����߸߸���
void antena_height_correction(pxyz px1, pxyz px2, penu pe, double*p1, double*p2, double*pp1, double*pp2);
//Ѱ����������
bool find_best_ephem(string prn, ptc ut, pnf nfile, pse ephemeris);
//������������
void cal_sate_coor(string prn, ptc ut, pnf pnavefile, bool &flag, pxyz coor, double* dt, double* reletivity);

//������������
void iterate_sate_coor(pxyz px, ptc ut, pse ephemeris, pxyz coor);

//����ۺ������Ӳ����
void sateclock(pse ephemeris, ptc ut, double clockoff, double delts);

//��ƽ����Ǽ���ƫ�����
double calE(double* M, double* e);
//��ƫ�������������
double calf(double* E, double* e);
//��ת����
//void rotamatrix(string str, double seta,Matrix r);
//������
void putresult(presult pt);
//spp����
void spp(pof obsfile, pnf pnavefile, presult pt);
#endif