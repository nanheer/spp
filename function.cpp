#include"函数.h"
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
//#include"matrix.h"
#include<iomanip>
#include<string>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using namespace std;
/*#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
#else
typedef matrix Matrix;
#endif

#ifndef _NO_EXCEPTION
#  define TRYBEGIN()	try {
#  define CATCHERROR()	} catch (const STD::exception& e) { \
						cerr << "Error: " << e.what() << endl; }
#else
#  define TRYBEGIN()
#  define CATCHERROR()
#endif*/


void readnfile(string strn, pnf pnavefile)
{
	ifstream nfile(strn, ios::in);
	if (!nfile){
		cout << "n文件打开错误" << endl;
		exit(0);
	}
	cout << "开始读n文件" << endl;
	string str1;
	getline(nfile, str1);
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		getline(nfile, str1);
	}
	//int a = 0;
	while (!nfile.eof())
	{
		string  str2;
		snf singlenf;
		se ephem;
		//第零行
		getline(nfile, str1);
		while (str1.substr(0, 1) != "G"&& nfile.eof() != true)
		{
			getline(nfile, str1);
			//cout << str1.substr(0, 1) << endl;
		}
		 //cout <<"出while"<< endl;
		// if (str1.empty() == true){ cout << "读取n文件终止" << endl; break; }//防止最后有一行空数据，读取o文件时先判断是否为">",再读取数据，所以不用加这一句
		if (nfile.eof() == true) break;//eof函数用来判断是否读到文件末尾，就算读到末尾，依旧可以用getline函数，不过不再往下读，如果最后一行不为空行，用empty函数判断，上面的while循环就出不去了
		str2 = str1.substr(0, 3);
		ephem.utime_n.year = atoi(str1.substr(4, 4).c_str());
		ephem.utime_n.month = atoi(str1.substr(9, 2).c_str());
		ephem.utime_n.day = atoi(str1.substr(12, 2).c_str());
		ephem.utime_n.hour = atoi(str1.substr(15, 2).c_str());
		ephem.utime_n.minute = atoi(str1.substr(18, 2).c_str());
		ephem.utime_n.second = atof(str1.substr(21, 2).c_str());
		ephem.a0 = atof(str1.substr(23, 19).c_str());
		ephem.a1 = atof(str1.substr(42, 19).c_str());
		ephem.a2 = atof(str1.substr(61, 19).c_str());
		//第一行
		getline(nfile, str1);
		ephem.IODE = atof(str1.substr(4, 19).c_str());
		ephem.Crs = atof(str1.substr(23, 19).c_str());
		ephem.deltn = atof(str1.substr(42, 19).c_str());
		ephem.M0 = atof(str1.substr(61, 19).c_str());
		//第二行
		getline(nfile, str1);
		ephem.Cuc = atof(str1.substr(4, 19).c_str());
		ephem.e = atof(str1.substr(23, 19).c_str());
		ephem.Cus = atof(str1.substr(42, 19).c_str());
		ephem.sqrt_a = atof(str1.substr(61, 19).c_str());
		//cout << str2<<"的长半轴的平方根：" << ephem.sqrt_a << endl;
		//第三行
		getline(nfile, str1);
		ephem.toe = atof(str1.substr(4, 19).c_str());
		ephem.Cic = atof(str1.substr(23, 19).c_str());
		ephem.omega0 = atof(str1.substr(42, 19).c_str());
		ephem.Cis = atof(str1.substr(61, 19).c_str());
		//第四行
		getline(nfile, str1);
		ephem.i0 = atof(str1.substr(4, 19).c_str());
		ephem.Crc = atof(str1.substr(23, 19).c_str());
		ephem.omega = atof(str1.substr(42, 19).c_str());
		ephem.omegadot = atof(str1.substr(61, 19).c_str());
		//第五行
		getline(nfile, str1);
		ephem.idot = atof(str1.substr(4, 19).c_str());
		//ephem.e = atof(str1.substr(24, 19));
		ephem.weekn = atof(str1.substr(42, 19).c_str());
		//ephem.sqrt_a = atof(str1.substr(62, 19));
		//第六行
		getline(nfile, str1);
		ephem.accuracy = atof(str1.substr(4, 19).c_str());
		ephem.health = atof(str1.substr(23, 19).c_str());
		ephem.Tgd = atof(str1.substr(42, 19).c_str());
		ephem.todc = atof(str1.substr(61, 19).c_str());
		//第七行
		getline(nfile, str1);
		ephem.transtime = atof(str1.substr(4, 19).c_str());
		//ephem.e = atof(str1.substr(24, 19));
		//ephem.Cus = atof(str1.substr(43, 19));
		//ephem.sqrt_a = atof(str1.substr(62, 19));

		singlenf.prn = str2;
		singlenf.epochn = 1;
		singlenf.sate.push_back(ephem);
		//vector<int>::iterator it;
		int flag = 0;
		int len = pnavefile->navefile.size();
		for (int i = 0; i<len; i++)
		{
			if (pnavefile->navefile[i].prn == singlenf.prn)
			{
				pnavefile->navefile[i].epochn = pnavefile->navefile[i].epochn + 1;
				pnavefile->navefile[i].sate.push_back(ephem);
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			//cout << "a的值：" << ++a << endl;
			pnavefile->navefile.push_back(singlenf);
		}
	}
	cout << "n文件观测历元数：" << pnavefile->navefile.size() << endl;
	nfile.close();
	cout << "读取n文件成功 " << endl;
}

void readofile(string stro, pof obsfile)
{
	int pos_p1=-1;//L1C在所给数据中的位置
	int pos_p2=-1;//L2W在所给数据中的位置
	ifstream ofile(stro, ios::in);
	if (!ofile){
		cout << "o文件打开错误" << endl;
		exit(0);
	}
	//memset(&obsfile->obsdata, 0, sizeof(obsfile->obsdata));
	//cout << "开始读o文件" << endl;
	od obsd;
	string str1,str2;
	getline(ofile, str1);
	//cout << str1 << endl;
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		//cout << "开始读文件头"<< endl;
		str2 = str1.substr(60, 19);
		//cout <<str2<< endl;
		if (str1.substr(0, 1) == "G" && str1.substr(60, 19) == "SYS / # / OBS TYPES")
		{
			int a = atoi(str1.substr(4, 2).c_str());
			//cout << "a的值：" <<a<< endl;
			for (int i = 0; i < a; i++)
			{
				if (str1.substr(7+4*i, 3) == "C1C") pos_p1 = i;
				if (str1.substr(7+4*i, 3) == "C2W") pos_p2 = i;
			}
			//cout << "p1位置的值：" << pos_p1 << endl;
			//cout << "p2位置的值：" << pos_p2 << endl;
			if (pos_p1 == -1)cout << "没有C1C数据" << endl;
			if (pos_p2 == -1)cout << "没有C2W数据" << endl;
		}

		if (str1.substr(60, 20) == "RINEX VERSION / TYPE"){
			obsfile->obsheaddata.obsdata_type = str1.substr(40, 1);
			//cout << str1.substr(60, 19) << endl;
		}
		if (str1.substr(60, 19) == "APPROX POSITION XYZ"){

			obsfile->obsheaddata.approx_coordinate.x = atof(str1.substr(1, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.y = atof(str1.substr(15, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.z = atof(str1.substr(29, 14).c_str());
			//cout << obsfile->obsheaddata.approx_coordinate.x << endl;
		}
		if (str1.substr(60, 20) == "ANTENNA: DELTA H/E/N"){
			obsfile->obsheaddata.antena_height.upping = atof(str1.substr(7, 7).c_str());
			obsfile->obsheaddata.antena_height.easting = atof(str1.substr(21, 7).c_str());
			obsfile->obsheaddata.antena_height.northing = atof(str1.substr(35, 7).c_str());
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		getline(ofile, str1);
	}
	cout << "读取文件头成功" << endl;
	while (!ofile.eof())
	{
		/*memset(&obsd.utime_o, 0, sizeof(obsd.utime_o));
		memset(&obsd.obs_n, 0, sizeof(int));
		memset(obsd.prn_list, 0, sizeof(string) * MAX_OBS_SATE_NUM);
		memset(obsd.obsvalue, 0, sizeof(double) *(MAX_OBS_TYPE*MAX_OBS_SATE_NUM));*/
		getline(ofile, str1);
		//if (str1.empty() == true) break;
		if (str1.substr(0, 1) == ">")
		{
			int sat = 0;
			obsd.utime_o.year = atoi(str1.substr(2, 4).c_str());
			obsd.utime_o.month = atoi(str1.substr(7, 2).c_str());
			obsd.utime_o.day = atoi(str1.substr(10, 2).c_str());
			obsd.utime_o.hour = atoi(str1.substr(13, 2).c_str());
			obsd.utime_o.minute = atoi(str1.substr(16, 2).c_str());
			obsd.utime_o.second = atof(str1.substr(19, 10).c_str());
			obsd.obs_n = atoi(str1.substr(32, 3).c_str());
			//cout << obsd.obs_n << endl;
			//cout << "观测时间：" << obsd.utime_o.hour << ":" << obsd.utime_o.minute << ":" << obsd.utime_o.second << obsd.obsvalue[sat][0] << endl;
			for (int i = 0; i < obsd.obs_n; i++)
			{
				//cout <<"i的值" <<i<< endl;
				getline(ofile, str1);
				if (str1.substr(0, 1) == "G")
				{
					
					//cout << str1 << endl;
					obsd.prn_list[sat] = str1.substr(0, 3);
					//cout << "obsd.obsvalue[" << sat <<"]的值：" << obsd.prn_list[sat] << endl;
					if (str1.length() < 16 *( pos_p1+1))//如果读取的一个字符串str1后面全是空格，则str1的长度只到最后一个字符串，用atof()函数转换之后的会出错
					{
						obsd.obsvalue[sat][0] = 0; 
					}
					else
					{
						obsd.obsvalue[sat][0] = atof(str1.substr(4 + 16 * pos_p1, 12).c_str());
					}
					//obsd.obsvalue[sat][1] = atof(str1.substr(4 + 16 * pos_p2, 12).c_str());
					if (str1.length() <112)
					{
						obsd.obsvalue[sat][1] = 0;
					}
					else
					{
					    obsd.obsvalue[sat][1] = atof(str1.substr(100, 12).c_str());
					}
					//cout << "obsd.obsvalue["<<sat<<"][0]的值" << obsd.obsvalue[sat][0] << endl;
					//cout << "obsd.obsvalue[" << sat << "][1]的值" << obsd.obsvalue[sat][1] << endl;
					//为了读取任意位置C1C和C2W，改一下

					/*for (int j = 0; j < MAX_OBS_TYPE; j++)
					{
						//cout << "j的值" <<j << endl;
						obsd.obsvalue[sat][j] = atof(str1.substr(4 + 16 * j, 12).c_str());
						//cout << "obsd.obsvalue["<<sat<<"]["<<j<<"]的值："<<obsd.obsvalue[sat][j] << endl;
					 }*/
					sat++;
				}
			}
			obsd.obs_n = sat;
			//cout << "obsd的观测时间：" << obsd.utime_o.hour << ":" << obsd.utime_o.minute << ":" << obsd.utime_o.second << endl;
			if (obsd.obs_n >= 4)
			{
				obsfile->obsdata.push_back(obsd);
				//cout << "新增加一个观测历元" << obsfile->obsdata.size()<< endl;
			}
		}

	}
	cout <<"o文件观测历元数：" <<obsfile->obsdata.size() << endl;
	ofile.close();
	cout << "读取o文件成功 "<< endl;
}

void gpsttoutc(pgpst gt, ptc ut)
{
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	gpsttojulianday(gt,ju);
	juliandaytoutc(ju, ut);
	free(ju);
}

void utctogpst(ptc ut, pgpst gt)
{
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	utctojulianday(ut, ju);
	juliandaytogpst(ju, gt);
	free(ju);
}

void utctojulianday(ptc ut, pjulian ju)
{
	int		m;
	int		y;
	double	dhour;

	dhour = ut->hour + ut->minute / (double)_HOUR_IN_MINUTE
		+ ut->second / (double)_HOUR_IN_SECOND;

	if (ut->month <= 2) {
		y = ut->year - 1;
		m = ut->month + 12;
	}
	else {
		y = ut->year;
		m = ut->month;
	}

	ju->daynum = (long)(365.25*y) + (long)(30.6001*(m + 1))
		+ ut ->day + (long)(dhour / 24 + 1720981.5);
	ju->secondnum = ((ut->hour + 12) % _DAY_IN_HOUR)*_HOUR_IN_SECOND
		+ ut->minute*_MINUTE_IN_SECOND + (long)ut->second;
	ju->secondfrac = ut->second-(long)ut->second;
}

void juliandaytoutc(pjulian ju, ptc ut)
{
	int a, b, c, d, e;
	double JD;
	JD = ju->daynum + (ju->secondnum + ju->secondfrac) / _DAY_IN_SECOND;

	a = static_cast<int>(JD + 0.5);
	b = a + 1537;
	c = static_cast<int>((b - 122.1) / 365.25);
	d = static_cast<int>(365.25*c);
	e = static_cast<int>((b - d) / 30.6001);
	
	double day = b - d-(long)(30.6001*e) + JD + 0.5 - a;
	ut->day = int(day);
	ut->month = e - 1 - 12 * (int)(e / 14);
	ut->year = c - 4715 - (int)((7 + ut->month) / 10);

	ut->hour = int((day - ut->day)*24.0);
	ut->minute = (int)(((day - ut->day)*24.0 - ut->hour)*60.0);
	ut->second = ju->secondnum + ju->secondfrac - (int((ju->secondnum + ju->secondfrac)/60))*60.0;

}

void gpsttojulianday(pgpst gt, pjulian ju)
{
	double JD;
	JD = gt->weeknum * 7 + (gt->secondnum + gt->secondfrac) / _DAY_IN_SECOND + 2444244.5;
	ju->daynum= long(JD);
	
	ju->secondnum = long(gt->secondnum + (gt->weeknum * 7 + 2444244.5 - ju->daynum)*_DAY_IN_SECOND);
	ju->secondfrac = gt->secondfrac;
}


void juliandaytogpst(pjulian ju, pgpst gt)
{
	double JD;
	JD = ju->daynum +( ju->secondnum + ju->secondfrac) / _DAY_IN_SECOND;
	gt->weeknum = int((JD - 2444244.5) / 7);

	gt->secondnum = long((JD - 2444244.5 - gt->weeknum * 7)*_DAY_IN_SECOND);
	gt->secondfrac = ju->secondfrac;

}

//求儒略日差值
double deltjulianday(pjulian ju1, pjulian ju2)
{
	double delt,d1,d2;
	d1 = ju1->daynum + (ju1->secondnum + ju1->secondfrac) / _DAY_IN_SECOND;
	d2 = ju2->daynum + (ju2->secondnum + ju2->secondfrac) / _DAY_IN_SECOND;
	//delt = (ju1->daynum - ju2->daynum)*_DAY_IN_SECOND + (ju1->secondnum - ju2->secondnum) + (ju1->secondfrac - ju2->secondfrac);
	delt = (d1 - d2)*_DAY_IN_SECOND;

	/*if (delt>302400)
		delt -= 604800;
	else if (delt<-302400)
		delt += 604800;
	else
		delt = delt;*/
	return delt;
}

//儒略日变换
void transjulian(pjulian ju1, double * dt, pjulian ju2)
{
	double JDold, JDnew;
	JDold = ju1->daynum + (ju1->secondnum + ju1->secondfrac) / _DAY_IN_SECOND;
	JDnew = JDold-*dt / _DAY_IN_SECOND;

	ju2->daynum = long(JDnew);
	ju2->secondnum = long((JDnew - long(JDnew))*_DAY_IN_SECOND);
	ju2->secondfrac = (JDnew - long(JDnew))*_DAY_IN_SECOND
		- long((JDnew - long(JDnew))*_DAY_IN_SECOND);
}






//坐标函数

void xyztoblh(pxyz px, pblh pb)
{
	double pi = 4.0*atan(1.0);
	double E2 = 2.0*flattening - flattening * flattening;
	double E4 = E2*E2;
	double ALFA = (px->x*px->x + px->y*px->y + (1.0 - E2)*px->z*px->z) / (a*a);
	double BATA = (px->x*px->x + px->y*px->y - (1.0 - E2)*px->z*px->z) / (a*a);
	double Q = 1.0 + 13.50*E4*(ALFA*ALFA - BATA*BATA) / pow(ALFA - E4, 3);
	double A1 = -Q + sqrt(Q*Q - 1.0);
	double AL = (1.0 / 3.0)*log(-A1);
	AL = -exp(AL);
	double A2 = AL + 1.0 / AL;
	double A3 = AL - 1.0 / AL;
	double T23 = (ALFA + E4 / 2.0) / 3.0 - (ALFA - E4)*A2 / 12.0;
	double T32 = sqrt(T23 *T23 + ((ALFA - E4)*A3)*((ALFA - E4)*A3) / 48.0);
	double T1 = -E2*BATA / (4.0*T32);
	double DK = T1 + sqrt(T23 + T32) - (1.0 - E2 / 2.0);
	double EK = (1.0 + DK) / (1.0 + DK - E2);
	/*cout << "T1的值：" << T1 << endl;
	cout << "T32的值：" << T32 << endl;
	cout << "T23的值：" << T23 << endl;
	cout << "E2的值：" << E2 << endl;
	cout << "ALFA的值：" << ALFA << endl;
	cout << "AL的值：" << AL << endl;
	cout << "DK的值：" << DK << endl;
	cout << "EK的值：" << EK << endl;*/
	pb->height = (DK / (1.0 + DK))*sqrt(px->x*px->x + px->y*px->y + (EK*px->z)*(EK*px->z));
	double P = sqrt(px->x*px->x + px->y*px->y);
	pb->latitude = atan(EK*px->z / P);
	double COSFL = px->x / P;
	double SINFL = px->y / P;
	pb->longitude = asin(SINFL);
	if (SINFL > 0.0&&COSFL<0.0) pb->longitude = pi - pb->longitude;
	if (SINFL < 0.0&&COSFL>0.0) pb->longitude = 2.0*pi + pb->longitude;
	if (SINFL <0.0&&COSFL<0.0) pb->longitude = pi - pb->longitude;

	/*double e2;//第一偏心率的平方
	e2 = 2 * flattening - flattening*flattening;

	pb->longitude = atan(px->y / px->x);
	double W, N, N1 = 0, B, B1;
	B1 = atan(px->z / sqrt(px->x*px->x + px->y*px->y));
	while (1)
	{
		W = sqrt(1 - e2*sin(B1)*sin(B1));
		N1 = a / W;
		B = atan((px->z + N1*e2*sin(B1)) / sqrt(px->x*px->x + px->y*px->y));

		if (fabs(B - B1)<delta)
			break;
		else
			B1 = B;
	}

	pb->latitude = B;
	N = a / sqrt(1 - e2*sin(pb->latitude)*sin(pb->latitude));
	pb->height = sqrt(px->x*px->x + px->y*px->y) / cos(B) - N;*/
}

void blhtoxyz(pblh pb, pxyz px)
{
	double e2;//第一偏心率的平方
	double N;//卯酉圈半径
	e2 = 2 * flattening - flattening*flattening;
	N = a / sqrt(1 - e2*sin(pb->latitude)*sin(pb->latitude));

	px->x = (N + pb->height)*cos(pb->latitude)*cos(pb->longitude);
	px->y = (N + pb->height)*cos(pb->latitude)*sin(pb->longitude);
	px->z = (N*(1 - e2) + pb->height)*sin(pb->latitude);
}

void xyztoenu(pxyz pxcenter, pxyz px, penu pe)
{
	double dx, dy, dz;
	dx = px->x - pxcenter->x;
	dy = px->y - pxcenter->y;
	dz = px->z - pxcenter->z;

	pblh pd;
	pd = (pblh)malloc(sizeof(BLH));

	xyztoblh(pxcenter, pd);

	pe->northing = -sin(pd->latitude)*cos(pd->longitude)*dx
		- sin(pd->latitude)*sin(pd->longitude)*dy
		+ cos(pd->latitude)*dz;
	pe->easting = -sin(pd->longitude)*dx
		+ cos(pd->longitude)*dy;
	pe->upping = cos(pd->latitude)*cos(pd->longitude)*dx
		+ cos(pd->latitude)*sin(pd->longitude)*dy
		+ sin(pd->latitude)*dz;
	free(pd);
}

 void enutoxyz(pxyz pxcenter, penu pe, pxyz px)
{
	pblh pd;
	pd = (pblh)malloc(sizeof(BLH));
	xyztoblh(pxcenter, pd);
	MatrixXd H(3, 3), DB(3, 1), DX(3, 1);
	DB(0, 0) = pe->northing;
	DB(1, 0) = pe->easting;
	DB(2, 0) = pe->upping;
	H(0, 0) = -sin(pd->latitude)*cos(pd->longitude);
	H(0, 1) = -sin(pd->latitude)*sin(pd->longitude);
	H(0, 2) = cos(pd->latitude);
	H(1, 0) = -sin(pd->longitude);
	H(1, 1) = cos(pd->longitude);
	H(1, 2) = 0;
	H(2, 0) = cos(pd->latitude)*cos(pd->longitude);
	H(2, 1) = cos(pd->latitude)*sin(pd->longitude);
	H(2, 2) = sin(pd->latitude);
	DX = (H.inverse())*DB;
	double dx, dy, dz;
	dx = DX(0,0 );
	dy = DX(1, 0);
	dz = DX(2, 0);
	px->x = pxcenter->x + dx;
	px->y = pxcenter->y + dy;
	px->z = pxcenter->z + dz;
}
 void enutoenupolar(penu pe, penupolar pep)
 {
	 pep->range = sqrt(pe->northing*pe->northing + pe->easting*pe->easting + pe->upping*pe->upping);
	 pep->azimuth = atan(pe->northing / pe->easting);
	 pep->elevation = atan(pe->upping / sqrt(pe->northing*pe->northing + pe->easting*pe->easting));
 }
 void tropo(double mjd, pxyz px1, pxyz px2, double * dx)//萨斯塔莫宁模型
 {
	 BLH crdSite;
	 xyztoblh(px1,&crdSite);
	 double lat = crdSite.latitude;
	 double lon = crdSite.longitude;
	 double hgt = crdSite.height;
	 ENU  pct;
	 xyztoenu( px1,  px2, &pct);
	 ENUPOLAR site;
	 enutoenupolar(&pct, &site);
	 double elev = site.elevation;
	 double pres, temp, rhumity, undu;
	 get_nominal_metedata(mjd, lat, lon, hgt, &pres, &temp, &rhumity, &undu);
	 double es, E, a, delev;
	 es = undu * exp(-37.2465 + 0.213166 * temp - 0.000256908 * temp *temp);
	 a = 1.16 - 0.00015 * hgt + 0.00000000716 * hgt*hgt;
	 delev = 16.0 / 206264.806247 / temp * (pres + 4810.0 / temp * es) / tan(elev);
	 E = elev + delev;
	 *dx = 0.002277 / sin(E) * (pres + (1255.0 / temp + 0.05) * es - a / (tan(E) *tan(E)));
	 
	/* double T = T0 - 0.0065*(H - H0) + 273.16;
		 double P = P0*pow(1 - 0.0000226*(H - H0), 5.225);
		 double RH = RH0*exp(-0.0006396*(H - H0));

		 double e = RH*exp(-37.2465 + 0.213166*T - 0.000256908*T*T);
		 double hw = 11000;
		 double hd = 40136 + 148.72*(T - 273.16);

		 double Kw = (155.2e-7 * 4810 * e*(hw - H)) / (T*T);
		 double Kd = (155.2e-7*P*(hd-H)) / T;
		 cout << "Kd的值：" << Kd << endl;
		 cout << "Kw的值：" << Kw << endl;
		 cout << "H的值：" << H << endl;
		 cout << "(1 - 0.0000226*(H - H0)的值：" << 1 - 0.0000226*(H - H0) << endl;
		 cout << "P的值：" << P << endl;
		 cout << "pow((1 - 0.0000226*(H - H0)), 5.225)的值：" << pow((1 - 0.0000226*(H - H0)), 5.225) << endl;
		 *dx = Kd /sin(sqrt(E*E + 6.25)) + Kw /sin(sqrt(E*E + 2.25));
		*/

 }


 bool find_best_ephem(string prn, ptc ut, pnf nfile, pse ephemeris)
 {
	 //bool flag =false;
	 //	bool satn = false;
	 JULIANDAY jld1, jld2;//jlld1卫星信号发射时间，jlld2卫星信星历时间
	 utctojulianday(ut, &jld1);
	 /*cout << "卫星信号发射时间：" << ut->year<< ":" << ut->month << ":" << ut->day;
	 cout << ":" << ut->hour << ":" << ut->minute << ":" << ut->second << endl; */
	 int k = nfile->navefile.size();
	 int n = 0;
	 for (int i = 0; i < k; i++)
	 {
		 if (prn == nfile->navefile[i].prn)
		 {
			 // satn = true;
			 for (int j = 0; j < nfile->navefile[i].epochn; j++)
			 {
				 utctojulianday(&nfile->navefile[i].sate[j].utime_n, &jld2);
				 //cout << "星历时间：" << nfile->navefile[i].sate[j].utime_n.year << ":" << nfile->navefile[i].sate[j].utime_n.month << ":" << nfile->navefile[i].sate[j].utime_n.day << ":";
				 //cout << nfile->navefile[i].sate[j].utime_n.hour << ":" << nfile->navefile[i].sate[j].utime_n.minute << ":" << nfile->navefile[i].sate[j].utime_n.second<<endl;
				 double delta;
				 delta = deltjulianday(&jld1, &jld2);
				// cout << "时间差：" << delta << endl;
				 /* if (delta<=7200)
				  {
				  // flag = true;
				  *ephemeris = nfile->navefile[i].sate[j];
				  return true;
				  break;
				  }
				  /* else if (delta<2 * 60 * 60)
				  {
				  // flag = true;
				  *ephemeris = nfile->navefile[i].sate[j];
				  return true;

				  }*/
				 /* if (delta >= 0 && delta<7200)
					 // if (fabs(delta )<= 3600)
					 {
					 *ephemeris = nfile->navefile[i].sate[j];
					 //flag = 0;
					 // cout << "返回找到" << endl;
					 return true;
					 }*/
				 if (delta >= 0 && delta <= 7230)
				 {
					 *ephemeris = nfile->navefile[i].sate[j];
					 return true;
				 }
				 else if (delta > -7200 && delta < 0)
				 {
					 *ephemeris = nfile->navefile[i].sate[j];
					 return true;
				 }
			 }
			 /*	 if (satn == true)
				  {
				  *ephemeris = nfile->navefile[i].sate[n];
				  return true;
				  break;
				  }*/
		 }
	 }
	// cout
	 return false;
 }


 //由平近点角计算偏近点角
 double calE(double* M, double* e)
 {
	 double E0 = *M;
	 double E;
	 E = *M + (*e)*sin(E0);
	 while( fabs(E0 - E) >1e-12)
	 {
		 E0 = E;
		 E = *M + (*e)*sin(E0);
	 }
	 return E;
 }

 //由偏近点角求真近点角
 double calf(double* E, double* e)
 {
	 double a, b,f;
	 a = (cos(*E) - *e) / (1 - (*e)*cos(*E));
	 b = (sqrt(1 - *e*(*e))*sin(*E)) / (1 - (*e)*cos(*E));
	 if (a == 0 && b > 0) f = PI / 2;
	 if (a == 0 && b < 0) f = PI*3 / 2;
	 if (a < 0 && b == 0) f = PI;
	 if (a > 0 && b == 0) f = 2*PI;
	// if (a ==0 && b == 0) f = 0;

	 if (a > 0 && b > 0) f = atan(b/a);
	 if (a > 0 && b < 0) f = atan(b/a)+2*PI;
	 if (a < 0 && b > 0) f = atan(b/a)+PI;
	 if (a < 0 && b < 0) f = atan(b/a)+PI;
	 return f;

	 /*double x,y,z;
	 z = (cos(*E) - *e) / (1 - (*e)*cos(*E));
	 y = (sqrt(1 - *e*(*e))*sin(*E)) / (1 - (*e)*cos(*E));
	 if (z == 0) x = PI / 2;
	 else{
		 if (y == 0) x = PI;
		 else{
			 x = atan(fabs(y / z));
			 if ((y>0) && (z<0)) x = PI - x;
			 else if ((y<0) && (z<0)) x = PI + x;
			 else if ((y<0) && (z>0)) x = 2 * PI - x;
		 }
	 }
	 return x;*/
 }
 
 
 //计算卫星坐标
 void cal_sate_coor(string prn, ptc ut, pnf pnavefile, bool &flag, pxyz coor, double* dt, double* reletivity)
 {
	 pse ephemeris=new se; 
	 if (find_best_ephem(prn, ut, pnavefile, ephemeris))
	 {
	 double tk, n0, n, E, M, f, u0, u, du, r, dr, i, di, l;
	 JULIANDAY t, toe;
	 utctojulianday(ut, &t);
	 utctojulianday(&ephemeris->utime_n, &toe);
	// cout << "findepoch" <<prn<< endl;
	//cout << "观测时间儒略日：" << t.daynum << ":" << t.secondnum << ":" << t.secondfrac<<  endl;
    //cout << "星历时间儒略日：" << toe.daynum << ":" << toe.secondnum << ":" << toe.secondfrac <<endl;
	 //cout << "找到观测时间utc：" << ut->hour << ":" << ut->minute << ":" << ut->second << endl;
	// cout << "找到星历时间utc：" << ephemeris->utime_n.hour << ":" << ephemeris->utime_n.minute << ":" << ephemeris->utime_n.second << endl;
	 GPST tg;
	 utctogpst(ut, &tg);

	 tk = deltjulianday(&t, &toe);
	// cout << "tk的值" << tk << endl;
	 *dt = ephemeris->a0 + ephemeris->a1*tk + ephemeris->a2*tk*tk;
	 n0 = sqrt(GM) / (ephemeris->sqrt_a*ephemeris->sqrt_a*ephemeris->sqrt_a);
	 n = n0 + ephemeris->deltn;
	 M = ephemeris->M0 + n*tk;
	 E = calE(&M, &ephemeris->e);
	 *reletivity = -2 * sqrt(GM)*ephemeris->sqrt_a*ephemeris->e*sin(E) / c;
	 f = calf(&E, &ephemeris->e);
	 u0 = f + ephemeris->omega;
	 du = ephemeris->Cuc*cos(2 * u0) + ephemeris->Cus*sin(2 * u0);
	 dr = ephemeris->Crc*cos(2 * u0) + ephemeris->Crs*sin(2 * u0);
	 di = ephemeris->Cic*cos(2 * u0) + ephemeris->Cis*sin(2 * u0);
	 u = u0 + du;
	 r = (ephemeris->sqrt_a*ephemeris->sqrt_a)*(1 - ephemeris->e*cos(E)) + dr;
	 // cout << "r的值" << r << endl;
	 i = ephemeris->i0 + di + ephemeris->idot*tk;
	 l = ephemeris->omega0 + ephemeris->omegadot*tk - we*(tg.secondnum + tg.secondfrac);
	 //MatrixXd coor1(3, 1), coor2(3, 1), coor3(3, 1), rz(3, 3), rx(3, 3);
	 double x1, y1, z1;
	 x1 = r*cos(u), y1 = r*sin(u), z1 = 0;

	 coor->x = cos(l)*x1 - cos(i)*sin(l)*y1;
	 coor->y = sin(l)*x1 + cos(i)*cos(l)*y1;
	 coor->z = sin(i)*y1;


	/* if (ut->hour == 12 && ut->minute == 15 && ut->second == 0.0)
	 {
		 cout << "tk的值" << tk << endl;
		 cout << "n0的值" << n0 << endl;
		 cout << "n的值" << n << endl;
		 cout << "E的值" << E << endl;
		 cout << "M的值" << M << endl;
		 cout << "f的值" << f << endl;
		 cout << "u0的值" << u0 << endl;
		 cout << "u的值" << u << endl;
		 cout << "du的值" << du << endl;
		 cout << "r的值" << r << endl;
		 cout << "dr的值" << dr << endl;
		 cout << "i的值" << i << endl;
		 cout << "di的值" << di << endl;
		 cout << "l的值" << l << endl;
		 cout << "toe的值" << toe.daynum << ":" << toe.secondnum << ":" << toe.secondfrac << endl;
		 cout << "t的值" << t.daynum << ":" << t.secondnum << ":" << t.secondfrac << endl;
		 cout << "x1的值" << x1 << endl;
		 cout << "y1的值" << y1 << endl;
		 cout << "cos(l)的值" << cos(l) << endl;
		 cout << "sin(l)的值" << sin(l) << endl;
		 cout << "cos(i)的值" << cos(i) << endl;
		 cout << "sin(i)的值" << sin(i) << endl;
		 cout << "coor->x的值" << coor->x << endl;
	 }*/
	 flag= true;
	 }
	 else
	 {
		//cout << "没有找到星历" << prn<<endl;
		//cout << "没有找到星历的观测时间utc：" << ut->hour << ":" << ut->minute << ":" << ut->second << endl;
		 flag= false;
		// cout << "flag的值" << flag<<endl;
		 coor->x = 0.0;
		 coor->y = 0.0;
		 coor->z = 0.0;
		 *dt = 0;
		 *reletivity = 0;
	 }
 }
 
 void putresult(presult pt)
 {
	 ofstream outfile("E:\\单点定位程序\\结果\\结果.txt", ios::out);
	 if (!outfile)
	 {
		 cerr << "文件创建失败" << endl;
		 abort();
	 }
	 cout << "开始写文件" << endl;
	 outfile << " 年 ," << "月," << "日," << "时," << "分," << "    秒,   " << "        X坐标 ," << "        Y坐标 ," << "        Z坐标 ,";
	 outfile << "        N坐标 ," << "        E坐标 ," << "        U坐标 ,";
	 outfile << "        pdop " << "         卫星个数";
	for (int k = 1; k < 10; k++)
	 { 
		 outfile << "  残差" << k<<"   ";
	 }
	for (int k = 10; k < 13; k++)
	{
		outfile << " 残差" << k<<"   ";
		if (k == 12) outfile << endl;
	}
	 cout << "结果个数："<<pt->allresult.size() << endl;
	 for (unsigned int i = 0; i < pt->allresult.size(); i++)
	 {
		 outfile.setf(ios_base::fixed, ios_base::floatfield);
		 outfile.precision(4);
		 outfile << setw(4) << pt->allresult[i].rtime.year << setw(3) << pt->allresult[i].rtime.month << setw(3) << pt->allresult[i].rtime.day;
		 outfile << setw(3) << pt->allresult[i].rtime.hour << setw(3) << pt->allresult[i].rtime.minute << setw(9) << pt->allresult[i].rtime.second;
		 outfile << std::right << setw(15) << pt->allresult[i].xyzresult.x << std::right << setw(15) << pt->allresult[i].xyzresult.y << std::right << setw(15) << pt->allresult[i].xyzresult.z;
		 outfile << std::right << setw(15) << pt->allresult[i].enuresult.northing << std::right << setw(15) << pt->allresult[i].enuresult.easting << std::right << setw(15) << pt->allresult[i].enuresult.upping;
		 outfile << std::right << setw(15) << pt->allresult[i].pdop << std::right << setw(15) << pt->allresult[i].n;
		 for (int k = 0; k < 12; k++)
		 {
			 outfile << std::right << setw(10) << pt->allresult[i].res[k];
			 if (k == 11) outfile << endl;
		 }
	 }
	 outfile.close();
	 cout << "文件写入完毕" << endl;
 }

 void spp(pof obsfile, pnf pnavefile, presult pt)
 {
	 cout << "开始计算" << endl;
	 int k = 0;
	 sresult sr;
	
	 pxyz pxcenter =new XYZ;
	 
	 pxcenter->x = xx;
	 pxcenter->y = yy;
	 pxcenter->z = zz;
	 penu pe = new ENU;
	 
	 pxyz pxy = new XYZ;
	 pjulian ju1 = new JULIANDAY;
	 pjulian ju2 = new JULIANDAY;
	 ptc ut = new UTC;
	 pxyz px = new XYZ;
	 pxyz px1 = new XYZ;

	// vector<errs> ers;
	// errs er;
	 double dt1 = 0;//卫星钟差
	 double dt2 = 0;//接收机钟差

	 double p0;//双频消除伪距
	 double dx;//对流层误差
	 double R;//观测站到卫星距离
	 JULIANDAY tjul;//计算儒略日
	 double mjd;//平儒略日
	 double pp1;//修正天线高之后p1伪距
	 double pp2;//修正天线高之后p2伪距
	 //XYZ antena_coor;//天线中心笛卡尔坐标
	// XYZ vec1;//地面到天线的矢量
	 //XYZ vec2;//天线到卫星的矢量
	 //double pr;//天线到卫星的矢量的长度
	 for (unsigned int i = 0; i < obsfile->obsdata.size(); i++)
	 {
		 double err = 0;
		// double residual=0.0;
		//cout << "i的值：" << i << endl;
		int count = 0;
		/*pxy->x = 0.0;
		pxy->y = 0.0;
		pxy->z = 0.0;*/

		pxy->x = xx;
		pxy->y = yy;
		pxy->z = zz;
		dt2 = 0.0;
		sr.pdop = 0.0;
		for (int r = 0; r < 15; r++)
		{
			sr.res[r] = 0.0;
		}
		double rg = 0.0;

		MatrixXd  Q(4, 4), V(4, 1);
		int len = 0;
		int len1 = 0;
		//cout << "观测时间输入：" << obsfile->obsdata[i].utime_o.day << ":" << obsfile->obsdata[i].utime_o.hour << ":" << obsfile->obsdata[i].utime_o.minute << ":" << obsfile->obsdata[i].utime_o.second << endl;
			do{
				//cout << "进入do循环" << endl;
				vector<errs> ers;
				errs er;
				double sum = 0;
				
				for (int j = 0; j < obsfile->obsdata[i].obs_n; j++)
				{
					//cout << "进入for循环" << endl;
					bool iflag = true;
					if (obsfile->obsdata[i].obsvalue[j][0] == 0 || obsfile->obsdata[i].obsvalue[j][1] == 0)continue;
					p0 = obsfile->obsdata[i].obsvalue[j][0] * ((p1*p1) / (p1*p1 - p2*p2)) - obsfile->obsdata[i].obsvalue[j][1] * ((p2*p2) / (p1*p1 - p2*p2));
					double t1 = 0, t2 = p0 / c, t3 = 0, rele;
					int co = 0;
					while (fabs(t1 - t2) > 1e-9)
					{
						//cout << "进入while循环" << endl;
						co++;
						t1 = t2;
						t3 = t1 + dt2;
						utctojulianday(&obsfile->obsdata[i].utime_o, ju1);
						transjulian(ju1, &t3, ju2);
						juliandaytoutc(ju2, ut);
						//cout << "t1：" << t1 << endl;
						//cout << "观测时间输入：" << obsfile->obsdata[i].utime_o.day << ":"<< obsfile->obsdata[i].utime_o.hour << ":" << obsfile->obsdata[i].utime_o.minute << ":" << obsfile->obsdata[i].utime_o.second << endl;
						//cout << "开始计算卫星坐标" << endl;
						//cout << "zh卫星信号发射时间：" << ut->day << ":" << ut->hour << ":" << ut->minute << ":" << ut->second << endl;
						//cout << "观测时间儒略日：" << ju1->daynum<<endl;
						//cout << "卫星信号发射时间儒略日：" << ju2->daynum << endl;
						cal_sate_coor(obsfile->obsdata[i].prn_list[j], ut, pnavefile, iflag, px, &dt1, &rele);
						if (!iflag)
						{
							//cout << "计算卫星坐标失败" << endl;
							break;

						}
						px1->x = cos(t2*we)*px->x + sin(t2*we)*px->y;
						px1->y = -sin(t2*we)*px->x + cos(t2*we)*px->y;
						px1->z = px->z;

						rg = sqrt((pxy->x - px1->x)*(pxy->x - px1->x) + (pxy->y - px1->y)*(pxy->y - px1->y) + (pxy->z - px1->z)*(pxy->z - px1->z));
						//cout << "rg：" << rg << endl;
						t2 = rg / c;
						if (co > 15)cout << "出不了while循环" << endl;
					}
					if (!iflag)
					{
						//cout << "计算下一颗卫星" << endl;
						continue;
					}

					//cout << "count的值" <<count<< endl;
					/*if (count == 0)
					{
						cout.setf(ios::fixed);
						cout << "观测时间：" << obsfile->obsdata[i].utime_o.hour << ":" << obsfile->obsdata[i].utime_o.minute << ":" << obsfile->obsdata[i].utime_o.second;
						cout << "卫星号：" << obsfile->obsdata[i].prn_list[j]<<endl;
						//cout << "卫星坐标：x:" << px1->x << "y:" << px1->y << "z:" << px1->z << endl;
					}*/
					//修改天线高引起的误差
					antena_height_correction(pxy, px1, &obsfile->obsheaddata.antena_height, &obsfile->obsdata[i].obsvalue[j][0], &obsfile->obsdata[i].obsvalue[j][1], &pp1, &pp2);
					p0 = pp1* ((p1*p1) / (p1*p1 - p2*p2)) - pp2*((p2*p2) / (p1*p1 - p2*p2));
					//地球自转改正
					px1->x = cos(t2*we)*px->x + sin(t2*we)*px->y;
					px1->y = -sin(t2*we)*px->x + cos(t2*we)*px->y;
					px1->z = px->z;

					R = sqrt((pxy->x - px1->x)*(pxy->x - px1->x) + (pxy->y - px1->y)*(pxy->y - px1->y) + (pxy->z - px1->z)*(pxy->z - px1->z));
					utctojulianday(&obsfile->obsdata[i].utime_o, &tjul);
					mjd = tjul.daynum + (tjul.secondnum + tjul.secondfrac)/_DAY_IN_SECOND;
					mjd -= deltjul;
					tropo(mjd,&obsfile->obsheaddata.approx_coordinate, px1, &dx);
					/*if (count == 0)
					{
						BLH crdSite;
						xyztoblh(pxy, &crdSite);
						double lat = crdSite.latitude;
						double lon = crdSite.longitude;
						double hgt = crdSite.height;
						cout.setf(ios::fixed);
						cout << "mjd：" << mjd << endl;
						cout << "纬度：" << crdSite.latitude << endl;
						cout << "经度：" << crdSite.longitude << endl;
						cout << "高度：" << crdSite.height << endl;
						cout << "P0：" << p0 <<endl;
						cout << "R：" << R << endl; 
						cout << "相对论：" << rele << endl;
						cout << "对流层：" << dx << endl;
						cout << "卫星钟差：" << dt1*c << endl;
						cout << "接收机钟差：" << dt2*c << endl;
					}*/
					er.L0 = p0 - R + c*dt1 + rele - c*dt2 - dx;
					//if (er.L0 > 3 * err) continue;
					er.B0 = (pxy->x - px1->x) / R;
					er.B1 = (pxy->y - px1->y) / R;
					er.B2 = (pxy->z - px1->z) / R;
					er.B3 = 1.0;
					er.res = 0;
					ers.push_back(er);
				}
	
					/* for (vector<errs>::iterator iter = ers.begin(); iter != ers.end();)
					 {
						 if (iter->L0 >3*err)
						 {
							 iter = ers.erase(iter);
						 }
						 else
						 {
							 iter++;
						 }
					 }*/
				// }
				 len = ers.size();
				if(len > 4)
				{
					MatrixXd  B(len, 4), L(len, 1),RES(len,1);
					for (int q = 0; q < len; q++)
					{
						L(q, 0) = ers[q].L0;
						B(q, 0) = ers[q].B0;
						B(q, 1) = ers[q].B1;
						B(q, 2) = ers[q].B2;
						B(q, 3) = 1;
					}
					/*if (count == 0)
					{
						cout << L << endl;
					}*/
					Q = (B.transpose()*B).inverse();
					V = ((B.transpose()*B).inverse())*(B.transpose()*L);
				
			
					RES = B*V - L;//残差
					if (count == 0)
					{
						for (int k = 0; k < len; k++)
						{
							sr.res[k] = L(k, 0);
						}
					}
					for (int k = 0; k < len; k++)
					{
						sum += RES(k, 0)*RES(k, 0);
					}
					err = sqrt(sum / len - 4);//残差平均值
				}
				//cout << "剔除之后len的长度：" << len << endl;
				if (len<=4)
				{
					break;
				}
				pxy->x += V(0, 0);
				pxy->y += V(1, 0);
				pxy->z += V(2, 0);
				dt2 += V(3, 0) / c;

				count++;
			} while ((fabs(V(0, 0)) > 1e-5 || fabs(V(1, 0)) > 1e-5 || fabs(V(2, 0)) > 1e-5)&&count<=10);
			//cout << "count的长度：" << count<<endl;
			if (count > 10)
			{
				continue;
			}
			//cout << "count的长度：" << count << endl;
		 if (len > 4)
		 {

			 for (int q = 0; q < 4; q++)
			 {
				 sr.pdop += Q(q, q);
			 }
			 sr.n = len;
			 sr.rtime = obsfile->obsdata[i].utime_o;
			 //cout << "观测时间：" << sr.rtime.hour << ":" << sr.rtime.minute<<":" << sr.rtime.second<< endl;
			 xyztoenu(pxcenter, pxy, pe);
			 sr.enuresult = *pe;
			 pxy->x += -xx;
			 pxy->y += -yy;
			 pxy->z += -zz;
			 sr.xyzresult = *pxy;
			 pt->allresult.push_back(sr);
			// cout << "k=" << k++ << endl;
		 }
	 }
	 delete ju1;
	 delete ju2;
	 delete ut;
	 delete px;
	 delete px1;

	 delete pxcenter ;
	 delete pe;
	 delete pxy;

 }

 //求气象数据
 void get_nominal_metedata(double dmjd, double dlat, double dlon,double dhgt, double* pres, double* temp, double* rhumity, double* undu)
 {
	 double dfac[20], aP[55], bP[55], P[10][10];
	// double ap_mean[55], bp_mean[55], ap_amp[55], bp_amp[55];
	 //double at_mean[55], bt_mean[55], at_amp[55], bt_amp[55];
	 //double a_geoid[55], b_geoid[55];
	 double doy, t, sum, hort;
	 double	 apm, apa, pres0;
	 double	 atm, ata, temp0;
	 int n, m, i, j, k, ir;
	 doy = dmjd-44239.0 + 1 - 28;
	 double a_geoid[55] =
	 { -5.6195e-001,-6.0794e-002, -2.0125e-001, -6.4180e-002, -3.6997e-002,
		 +1.0098e+001, +1.6436e+001, +1.4065e+001, +1.9881e+000, +6.4414e-001,
		 -4.7482e+000, -3.2290e+000, +5.0652e-001, +3.8279e-001, -2.6646e-002, 
		 +1.7224e+000, -2.7970e-001, +6.8177e-001, -9.6658e-002, -1.5113e-002, 
		 +2.9206e-003, -3.4621e+000, -3.8198e-001, +3.2306e-002, +6.9915e-003, 
		 -2.3068e-003, -1.3548e-003, +4.7324e-006, +2.3527e+000, +1.2985e+000, 
		 +2.1232e-001, +2.2571e-002, -3.7855e-003, +2.9449e-005, -1.6265e-004, 
		 +1.1711e-007, +1.6732e+000, +1.9858e-001, +2.3975e-002, -9.0013e-004, 
		 -2.2475e-003, -3.3095e-005, -1.2040e-005, +2.2010e-006, -1.0083e-006, 
		 +8.6297e-001, +5.8231e-001, +2.0545e-002, -7.8110e-003, -1.4085e-004, 
		 -8.8459e-006, +5.7256e-006, -1.5068e-006, +4.0095e-007, -2.4185e-008
	 };
	 double b_geoid[55] =
	 {
		 +0.0000e+000, +0.0000e+000, -6.5993e-002, +0.0000e+000, +6.5364e-002,
		 -5.8320e+000, +0.0000e+000, +1.6961e+000, -1.3557e+000, +1.2694e+000,
		 +0.0000e+000, -2.9310e+000, +9.4805e-001, -7.6243e-002, +4.1076e-002,
		 +0.0000e+000, -5.1808e-001, -3.4583e-001, -4.3632e-002, +2.2101e-003,
		 -1.0663e-002, +0.0000e+000, +1.0927e-001, -2.9463e-001, +1.4371e-003,
		 -1.1452e-002, -2.8156e-003, -3.5330e-004, +0.0000e+000, +4.4049e-001,
		 +5.5653e-002, -2.0396e-002, -1.7312e-003, +3.5805e-005, +7.2682e-005,
		 +2.2535e-006, +0.0000e+000, +1.9502e-002, +2.7919e-002, -8.1812e-003,
		 +4.4540e-004, +8.8663e-005, +5.5596e-005, +2.4826e-006, +1.0279e-006,
		 +0.0000e+000, +6.0529e-002, -3.5824e-002, -5.1367e-003, +3.0119e-005,
		 -2.9911e-005, +1.9844e-005, -1.2349e-006, -7.6756e-009, +5.0100e-008
	 };
	double ap_mean[55] =
	 {
		 +1.0108e+003, +8.4886e+000, +1.4799e+000, -1.3897e+001, +3.7516e-003, 
		 -1.4936e-001, +1.2232e+001, -7.6615e-001, -6.7699e-002, +8.1002e-003, 
		 -1.5874e+001, +3.6614e-001, -6.7807e-002, -3.6309e-003, +5.9966e-004, 
		 +4.8163e+000, -3.7363e-001, -7.2071e-002, +1.9998e-003, -6.2385e-004, 
		 -3.7916e-004, +4.7609e+000, -3.9534e-001, +8.6667e-003, +1.1569e-002, 
		 +1.1441e-003, -1.4193e-004, -8.5723e-005, +6.5008e-001, -5.0889e-001, 
		 -1.5754e-002, -2.8305e-003, +5.7458e-004, +3.2577e-005, -9.6052e-006, 
		 -2.7974e-006, +1.3530e+000, -2.7271e-001, -3.0276e-004, +3.6286e-003, 
		 -2.0398e-004, +1.5846e-005, -7.7787e-006, +1.1210e-006, +9.9020e-008, 
		 +5.5046e-001, -2.7312e-001, +3.2532e-003, -2.4277e-003, +1.1596e-004, 
		 +2.6421e-007, -1.3263e-006, +2.7322e-007, +1.4058e-007, +4.9414e-009 
	 };
	double bp_mean[55] =
	 {
		 +0.0000e+000, +0.0000e+000, -1.2878e+000, +0.0000e+000, +7.0444e-001, 
		 +3.3222e-001, +0.0000e+000, -2.9636e-001, +7.2248e-003, +7.9655e-003, 
		 +0.0000e+000, +1.0854e+000, +1.1145e-002, -3.6513e-002, +3.1527e-003, 
		 +0.0000e+000, -4.8434e-001, +5.2023e-002, -1.3091e-002, +1.8515e-003, 
		 +1.5422e-004, +0.0000e+000, +6.8298e-001, +2.5261e-003, -9.9703e-004, 
		 -1.0829e-003, +1.7688e-004, -3.1418e-005, +0.0000e+000, -3.7018e-001, 
		 +4.3234e-002, +7.2559e-003, +3.1516e-004, +2.0024e-005, -8.0581e-006, 
		 -2.3653e-006, +0.0000e+000, +1.0298e-001, -1.5086e-002, +5.6186e-003, 
		 +3.2613e-005, +4.0567e-005, -1.3925e-006, -3.6219e-007, -2.0176e-008, 
		 +0.0000e+000, -1.8364e-001, +1.8508e-002, +7.5016e-004, -9.6139e-005, 
		 -3.1995e-006, +1.3868e-007, -1.9486e-007, +3.0165e-010, -6.4376e-010
	 };

	 double ap_amp[55] =
	 {
		 -1.0444e-001, +1.6618e-001, -6.3974e-002, +1.0922e+000, +5.7472e-001, 
		 -3.0277e-001, -3.5087e+000, +7.1264e-003, -1.4030e-001, +3.7050e-002, 
		 +4.0208e-001, -3.0431e-001, -1.3292e-001, +4.6746e-003, -1.5902e-004, 
		 +2.8624e+000, -3.9315e-001, -6.4371e-002, +1.6444e-002, -2.3403e-003, 
		 +4.2127e-005, +1.9945e+000, -6.0907e-001, -3.5386e-002, -1.0910e-003, 
		 -1.2799e-004, +4.0970e-005, +2.2131e-005, -5.3292e-001, -2.9765e-001, 
		 -3.2877e-002, +1.7691e-003, +5.9692e-005, +3.1725e-005, +2.0741e-005, 
		 -3.7622e-007, +2.6372e+000, -3.1165e-001, +1.6439e-002, +2.1633e-004, 
		 +1.7485e-004, +2.1587e-005, +6.1064e-006, -1.3755e-008, -7.8748e-008, 
		 -5.9152e-001, -1.7676e-001, +8.1807e-003, +1.0445e-003, +2.3432e-004, 
		 +9.3421e-006, +2.8104e-006, -1.5788e-007, -3.0648e-008, +2.6421e-010
	 };

	 double bp_amp[55] =
	 {
		 +0.0000e+000, +0.0000e+000, +9.3340e-001, +0.0000e+000, +8.2346e-001, 
		 +2.2082e-001, +0.0000e+000, +9.6177e-001, -1.5650e-002, +1.2708e-003, 
		 +0.0000e+000, -3.9913e-001, +2.8020e-002, +2.8334e-002, +8.5980e-004, 
		 +0.0000e+000, +3.0545e-001, -2.1691e-002, +6.4067e-004, -3.6528e-005, 
		 -1.1166e-004, +0.0000e+000, -7.6974e-002, -1.8986e-002, +5.6896e-003, 
		 -2.4159e-004, -2.3033e-004, -9.6783e-006, +0.0000e+000, -1.0218e-001, 
		 -1.3916e-002, -4.1025e-003, -5.1340e-005, -7.0114e-005, -3.3152e-007, 
		 +1.6901e-006, +0.0000e+000, -1.2422e-002, +2.5072e-003, +1.1205e-003, 
		 -1.3034e-004, -2.3971e-005, -2.6622e-006, +5.7852e-007, +4.5847e-008, 
		 +0.0000e+000, +4.4777e-002, -3.0421e-003, +2.6062e-005, -7.2421e-005, 
		 +1.9119e-006, +3.9236e-007, +2.2390e-007, +2.9765e-009, -4.6452e-009
	 };
	 double at_mean[55] =
	 {
		 +1.6257e+001, +2.1224e+000, +9.2569e-001, -2.5974e+001, +1.4510e+000, 
		 +9.2468e-002, -5.3192e-001, +2.1094e-001, -6.9210e-002, -3.4060e-002, 
		 -4.6569e+000, +2.6385e-001, -3.6093e-002, +1.0198e-002, -1.8783e-003, 
		 +7.4983e-001, +1.1741e-001, +3.9940e-002, +5.1348e-003, +5.9111e-003, 
		 +8.6133e-006, +6.3057e-001, +1.5203e-001, +3.9702e-002, +4.6334e-003, 
		 +2.4406e-004, +1.5189e-004, +1.9581e-007, +5.4414e-001, +3.5722e-001, 
		 +5.2763e-002, +4.1147e-003, -2.7239e-004, -5.9957e-005, +1.6394e-006, 
		 -7.3045e-007, -2.9394e+000, +5.5579e-002, +1.8852e-002, +3.4272e-003, 
		 -2.3193e-005, -2.9349e-005, +3.6397e-007, +2.0490e-006, -6.4719e-008, 
		 -5.2225e-001, +2.0799e-001, +1.3477e-003, +3.1613e-004, -2.2285e-004, 
		 -1.8137e-005, -1.5177e-007, +6.1343e-007, +7.8566e-008, +1.0749e-009
	 };
	 double bt_mean[55] =
	 {
		 +0.0000e+000, +0.0000e+000, +1.0210e+000, +0.0000e+000, +6.0194e-001, 
		 +1.2292e-001, +0.0000e+000, -4.2184e-001, +1.8230e-001, +4.2329e-002, 
		 +0.0000e+000, +9.3312e-002, +9.5346e-002, -1.9724e-003, +5.8776e-003, 
		 +0.0000e+000, -2.0940e-001, +3.4199e-002, -5.7672e-003, -2.1590e-003, 
		 +5.6815e-004, +0.0000e+000, +2.2858e-001, +1.2283e-002, -9.3679e-003, 
		 -1.4233e-003, -1.5962e-004, +4.0160e-005, +0.0000e+000, +3.6353e-002, 
		 -9.4263e-004, -3.6762e-003, +5.8608e-005, -2.6391e-005, +3.2095e-006, 
		 -1.1605e-006, +0.0000e+000, +1.6306e-001, +1.3293e-002, -1.1395e-003, 
		 +5.1097e-005, +3.3977e-005, +7.6449e-006, -1.7602e-007, -7.6558e-008, 
		 +0.0000e+000, -4.5415e-002, -1.8027e-002, +3.6561e-004, -1.1274e-004, 
		 +1.3047e-005, +2.0001e-006, -1.5152e-007, -2.7807e-008, +7.7491e-009
	 };
	 double at_amp[55] =
	 {
		 -1.8654e+000, -9.0041e+000, -1.2974e-001, -3.6053e+000, +2.0284e-002, 
		 +2.1872e-001, -1.3015e+000, +4.0355e-001, +2.2216e-001, -4.0605e-003, 
		 +1.9623e+000, +4.2887e-001, +2.1437e-001, -1.0061e-002, -1.1368e-003, 
		 -6.9235e-002, +5.6758e-001, +1.1917e-001, -7.0765e-003, +3.0017e-004, 
		 +3.0601e-004, +1.6559e+000, +2.0722e-001, +6.0013e-002, +1.7023e-004, 
		 -9.2424e-004, +1.1269e-005, -6.9911e-006, -2.0886e+000, -6.7879e-002, 
		 -8.5922e-004, -1.6087e-003, -4.5549e-005, +3.3178e-005, -6.1715e-006, 
		 -1.4446e-006, -3.7210e-001, +1.5775e-001, -1.7827e-003, -4.4396e-004, 
		 +2.2844e-004, -1.1215e-005, -2.1120e-006, -9.6421e-007, -1.4170e-008, 
		 +7.8720e-001, -4.4238e-002, -1.5120e-003, -9.4119e-004, +4.0645e-006, 
		 -4.9253e-006, -1.8656e-006, -4.0736e-007, -4.9594e-008, +1.6134e-009
	 };
	 double bt_amp[55] =
	 {
		 +0.0000e+000, +0.0000e+000, -8.9895e-001, +0.0000e+000, -1.0790e+000, 
		 -1.2699e-001, +0.0000e+000, -5.9033e-001, +3.4865e-002, -3.2614e-002, 
		 +0.0000e+000, -2.4310e-002, +1.5607e-002, -2.9833e-002, -5.9048e-003, 
		 +0.0000e+000, +2.8383e-001, +4.0509e-002, -1.8834e-002, -1.2654e-003, 
		 -1.3794e-004, +0.0000e+000, +1.3306e-001, +3.4960e-002, -3.6799e-003, 
		 -3.5626e-004, +1.4814e-004, +3.7932e-006, +0.0000e+000, +2.0801e-001, 
		 +6.5640e-003, -3.4893e-003, -2.7395e-004, +7.4296e-005, -7.9927e-006, 
		 -1.0277e-006, +0.0000e+000, +3.6515e-002, -7.4319e-003, -6.2873e-004, 
		 -8.2461e-005, +3.1095e-005, -5.3860e-007, -1.2055e-007, -1.1517e-007, 
		 +0.0000e+000, +3.1404e-002, +1.5580e-002, -1.1428e-003, +3.3529e-005, 
		 +1.0387e-005, -1.9378e-006, -2.7327e-007, +7.5833e-009, -9.2323e-009
	 };
	 t = sin(dlat);
	 n = 9;
	 m = 9;

	 dfac[0] = 1;;
	 for (i = 1; i <=(2 * n + 1); i++)
	 {
		 dfac[i] = dfac[i-1] * i;
	 }

	 for (i = 0; i <= n; i++)
	 {
		 for (j = 0; j <= min(i, m); j++)
		 {
			 ir = int((i - j) / 2);
			 sum = 0;
			 for (k = 0; k <= ir; k++)
			 {
				 sum = sum + pow((-1),k)*dfac[2 * i - 2 * k] / dfac[k] / 
					 dfac[i - k] / dfac[i - j - 2 * k ]*pow(t,(i - j - 2 * k));
					
			 }
			  P[i][j] = 1.0 /pow( 2 ,i)*sqrt(pow(1 - pow(t,2),j))*sum;
		 }
	 }
	 i = 0;
	 for (n = 0; n <= 9; n++)
	 {
		 for (m = 0; m <= n; m++)
		 {
			 
			 aP[i] = P[n][m] * cos(m*dlon);
			 bP[i] = P[n][m] * sin(m*dlon);
			 i = i + 1;
		 }
	 }
	 *undu = 0.0;
	 for (i = 1; i <= 55; i++)
	 {
		 *undu = *undu + (a_geoid[i-1] * aP[i-1] + b_geoid[i-1] * bP[i-1]);
	 }
	 hort = dhgt - *undu;
	 apm = 0.0;
	 apa = 0.0;
	 for(i = 1; i <= 55; i++)
	 {
		 apm = apm + (ap_mean[i-1]*aP[i-1] + bp_mean[i-1]*bP[i-1]);
		 apa = apa + (ap_amp[i-1] *aP[i-1] + bp_amp[i-1] *bP[i-1]);
	 }
	 pres0 = apm + apa*cos(doy / 365.25*2.0*PI);
	 *pres = pres0*pow(1.0 - 0.0000226*hort,5.225);
	 atm = 0.0;
	 ata = 0.0;
	 for (i = 1; i <= 55; i++)
	 {
		 atm = atm + (at_mean[i-1]*aP[i-1] + bt_mean[i-1]*bP[i-1]);
		 ata = ata + (at_amp[i-1] *aP[i-1] + bt_amp[i-1] *bP[i-1]);
	 }
	 temp0 = atm + ata*cos(doy / 365.25 * 2 * PI);
	 *temp = temp0 - 0.0065*hort;

	 *rhumity = 0.60;
 }

 int min(int a, int b)
 {
	 if (a >= b){ return b; }
	 else{ return a; }
 }
 void antena_height_correction(pxyz px1, pxyz px2, penu pe, double*p1, double*p2, double*pp1, double*pp2)
 {
	 XYZ antena_coor;//天线中心笛卡尔坐标
	 XYZ vec1;//地面到天线的矢量
	 XYZ vec2;//天线到卫星的矢量
	 XYZ vec3;//天线到卫星的p1改正之后矢量，即长度为p1，方向仍是天线中心指向卫星的方向
	 XYZ vec4;//天线到卫星的p2改正之后矢量，即长度为p2，方向仍是天线中心指向卫星的方向
	 XYZ vec5;//地面到卫星的矢量，即地面到天线的矢量加上天线到卫星的p1改正之后矢量
	 XYZ vec6;//地面到卫星的矢量，即地面到天线的矢量加上天线到卫星的p2改正之后矢量
	 double pr;//天线到卫星的矢量的长度
	 enutoxyz(px1, pe, &antena_coor);
	 vec1.x = antena_coor.x - px1->x, vec1.y = antena_coor.y - px1->y, vec1.z = antena_coor.z - px1->z;
	 vec2.x = px2->x - antena_coor.x, vec2.y = px2->y - antena_coor.y, vec2.z = px2->z - antena_coor.z;
	 pr = sqrt(vec2.x*vec2.x + vec2.y*vec2.y + vec2.z*vec2.z);
	 
	 vec3.x = vec2.x*(*p1) / pr;
	 vec3.y = vec2.y*(*p1) / pr;
	 vec3.z = vec2.z*(*p1) / pr;
	
	 vec4.x = vec2.x*(*p2) / pr;
	 vec4.y = vec2.y*(*p2) / pr;
	 vec4.z = vec2.z*(*p2) / pr;

	 vec5.x = vec1.x + vec3.x;
	 vec5.y = vec1.y + vec3.y;
	 vec5.z = vec1.z + vec3.z;

	 vec6.x = vec1.x + vec4.x;
	 vec6.y = vec1.y + vec4.y;
	 vec6.z = vec1.z + vec4.z;

	 *pp1 = sqrt(vec5.x*vec5.x + vec5.y*vec5.y + vec5.z*vec5.z);
	 *pp2 = sqrt(vec6.x*vec6.x + vec6.y*vec6.y + vec6.z*vec6.z);
 }
