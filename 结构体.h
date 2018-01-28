#ifndef STRUCT_FILE
#define STRUCT_FILE
#include <vector>
using namespace std;
const int MAX_OBS_SATE_NUM =40;
const int MAX_OBS_TYPE = 2;
//ʱ��ṹ��
typedef struct utctime{
	int year,month,day,hour,minute;
	double second;
}UTC;
typedef UTC * ptc;

typedef struct gpstime{
	int weeknum;
	long secondnum;
	double secondfrac;
}GPST;
typedef GPST * pgpst;

typedef struct juliantime{
	long daynum;
	int secondnum;
	double secondfrac;
}JULIANDAY;
typedef JULIANDAY * pjulian;

typedef struct doytime{
	int yearnum;
	int daynum;
	int secondnum;
	double secondfrac;
}DOY;
typedef DOY * pdoy;

//����ṹ��
typedef struct xyz_coordinate{
	double x;
	double y;
	double z;
}XYZ;
typedef XYZ * pxyz;

typedef struct enu_coordinate{
	double northing;
	double easting;
	double upping;
}ENU;
typedef ENU * penu;

typedef struct enupolar_coordinate{
	double range;
	double azimuth;
	double elevation;;
}ENUPOLAR;
typedef ENUPOLAR * penupolar;

typedef struct blh_coordinate{
	double longitude;
	double latitude;
	double height;
}BLH;
typedef BLH * pblh;

//n�ļ��ṹ��
typedef struct satelite_ephemeris{
	//��0��
	UTC utime_n;//����ʱ��
	double a0, a1, a2;//�����Ӳ�
	//��1��
	double IODE;//�������ݵ�����
	double Crs;//������������ڸ���������
	double deltn;//���ǽ��ٶ��㶯����
	double M0;//�ο�ʱ��ƽ�����
	//��2��
	double Cuc;//�����Ǿ����������
	double e;//ƫ����
	double Cus;//�����Ǿ����������
	double sqrt_a;//������Ŀ���
	//��3��
	double toe;//�ο�����ʱ���������
	double Cic;//�����Ǹ���������
	double omega0;//�ο���������ʼʱ�̣�������ྭ��GAST֮��
	double Cis;//�����Ǹ���������
	//��4��
	double i0;//�ο�ʱ�̹�����
	double Crc;//������������ڸ���������
	double omega;//���ص�Ǿ�
	double omegadot;//������ྭ�仯��
	//��5��
	double idot;//�����Ǳ仯��
	//double L2CODES;L2����ָʾ
	double weekn;//GPSʱ������
	//double L2P;L2P��
	//��6��
	double accuracy;//��������
	double health;//����ָ��
	double Tgd;//Ⱥ�ӳ�
	double todc;//��������������
	//��7��
	double transtime;//����ʱ��
	//double health;
	//double Tgd;
	//double todc;
}se;
typedef se * pse;
typedef struct singlenavefile
{
	int epochn;
	string prn;
	vector<se> sate;
}snf;
typedef struct navigationfile
{
	vector<snf> navefile;
}nf;
typedef nf * pnf;

//o�ļ��ṹ��
typedef struct obs_data{
	UTC utime_o;
	int obs_n;
	string prn_list[MAX_OBS_SATE_NUM];
	double obsvalue[MAX_OBS_SATE_NUM][MAX_OBS_TYPE];
}od;
typedef struct obs_head_data{
	XYZ approx_coordinate;
	ENU antena_height;
	string obsdata_type;
}ohd;

typedef struct obsfile{
	ohd obsheaddata;
	vector<od> obsdata;
}of;
typedef of * pof;


//����ṹ��
typedef struct single_result_out{
	XYZ xyzresult;
	ENU enuresult;
	UTC rtime;
	double pdop;
	int n;
	double res[15];
}sresult;
typedef struct result_out{
	vector<sresult> allresult;
}result;
typedef result * presult;


//������޳��ṹ��
struct errs
{
	double L0;
	double B0, B1, B2, B3,res;
};
#endif