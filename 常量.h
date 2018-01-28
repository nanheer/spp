#ifndef CONSTANT_H
#define CONSTANT_H
const double PI = 3.1415926535898;//圆周率
const double c = 2.99792458e+08;//光速
const double GM = 3.986005e+14;//重力常数
const double we = 7.2921151467e-05;//地球自转角速度
const double F = -4.442807633e-10;
const double f1 = 1575420000.0;
const double f2 = 1227600000.0;
const long _HOUR_IN_MINUTE = 60L;
const long _MINUTE_IN_SECOND = 60L;
const long	_DAY_IN_HOUR = 24L;
const long	_WEEK_IN_DAY = 7L;
const long	_HOUR_IN_SECOND = _HOUR_IN_MINUTE * _MINUTE_IN_SECOND;
const long	_DAY_IN_SECOND = _DAY_IN_HOUR * _HOUR_IN_SECOND;
const long	_WEEK_IN_SECOND = _WEEK_IN_DAY * _DAY_IN_SECOND;

//WGS-84椭球体参数
const double a = 6378136.49;//长半轴
const double flattening = 1 / 298.25642;//扁率
const double delta = 0.0000001;

const double H0 = 0;//参考面上的高度,温度，气压，湿度
const double T0 = 20;
const double P0 = 1013.25;
const double RH0 = 0.50;

const double p1 = 1575420000.0;
const double p2 = 1227600000.0;

const double deltjul = 2400000.5;

/*const double xx = 5084657.6239;
const double yy = 2670325.3492;
const double zz = -2768480.9482;

const double xx = -2279829.6306;
const double yy = 5004706.2836;
const double zz = 3219776.9440;

const double xx = -2279828.96926868;
const double yy = 5004706.49112338;
const double zz = 3219777.42147052;*/
//TWTF
const double xx = -2994428.49808605;
const double yy = 4951309.09172896;
const double zz = 2674496.73881957;


#endif