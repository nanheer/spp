#ifndef STRUCT_FILE
#define STRUCT_FILE
#include <vector>
using namespace std;
const int MAX_OBS_SATE_NUM =40;
const int MAX_OBS_TYPE = 2;
//时间结构体
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

//坐标结构体
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

//n文件结构体
typedef struct satelite_ephemeris{
	//第0行
	UTC utime_n;//星历时间
	double a0, a1, a2;//卫星钟差
	//第1行
	double IODE;//星历数据的年龄
	double Crs;//轨道径向方向周期改正正弦项
	double deltn;//卫星角速度摄动参数
	double M0;//参考时刻平近点角
	//第2行
	double Cuc;//升交角距改正余弦项
	double e;//偏心率
	double Cus;//升交角距改正正弦项
	double sqrt_a;//长半轴的开方
	//第3行
	double toe;//参考星历时间的周内秒
	double Cic;//轨道倾角改正余弦项
	double omega0;//参考星历周起始时刻，升交点赤经与GAST之差
	double Cis;//轨道倾角改正正弦项
	//第4行
	double i0;//参考时刻轨道倾角
	double Crc;//轨道径向方向周期改正余弦项
	double omega;//近地点角距
	double omegadot;//升交点赤经变化率
	//第5行
	double idot;//轨道倾角变化率
	//double L2CODES;L2代码指示
	double weekn;//GPS时星期数
	//double L2P;L2P码
	//第6行
	double accuracy;//星历精度
	double health;//健康指标
	double Tgd;//群延迟
	double todc;//卫星钟数据年龄
	//第7行
	double transtime;//传送时间
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

//o文件结构体
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


//结果结构体
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


//中误差剔除结构体
struct errs
{
	double L0;
	double B0, B1, B2, B3,res;
};
#endif