package com.shbd.ddap.coordinate;


import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import com.shbd.ddap.coordinate.ConstantHolder.ArcFormat;
import com.shbd.ddap.coordinate.ConstantHolder.CoordinateType;
public class CoordTransform {

	/**
	 *  同一参考椭球下，空间直角坐标--> 大地坐标
	 * @param source
	 * @return 大地坐标（弧度）
	 */
	public static Coordinate ECEF2BLH(Coordinate source) {
		Coordinate target = new Coordinate();
		target.setRe(source.getRe());	// 设置转换后的椭球
		target.setCoordType(CoordinateType.GEODETIC);	// 设置转换后的坐标为大地坐标(B/L/H)

		double X = source.getX();
		double Y = source.getY();
		double Z = source.getZ();

		double L = Math.atan2(Y, X); // ? atan(),atan2()
		target.setX(L);

		double a = source.getRe().getA();
		double b = source.getRe().getB();
		double e = source.getRe().getE1();

		// 迭代初始值
		double N0 = a;
		double H0 = Math.sqrt(X * X + Y * Y + Z * Z) - Math.sqrt(a * b);
		double B0 = Math.atan(Z / (Math.sqrt(X * X + Y * Y) * (1 - e * e * N0 / (N0 + H0))));

		// 迭代开始
		while (true) {
			double N1 = a / Math.sqrt(1 - e * e * Math.pow(Math.sin(B0), 2));
			double H1 = Math.sqrt(X * X + Y * Y) / Math.cos(B0) - N1;
			double B1 = Math.atan(Z / (Math.sqrt(X * X + Y * Y) * (1 - e * e * N1 / (N1 + H1))));

			// 满足下列条件时，退出迭代。当B(i)-B(i-1)<=0.00001秒， H(i)-H(i-1)<=0.001米时
			if (((B1 - B0) <= as2Rad(0.00001)) && ((H1 - H0) <= 0.001)) {
				double B = B1;
				double H = H1;
				target.setY(B);
				target.setZ(H);

				return target;
			} else {
				N0 = N1;
				H0 = H1;
				B0 = B1;
			}
		}
	}
	
	/**
	 *  同一参考椭球下，空间直角坐标--> 大地坐标
	 * @param source
	 * @return 大地坐标（度）
	 */ 
	public static Coordinate ECEF2BLH(Coordinate source, ArcFormat arcFormat) {
		Coordinate target = ECEF2BLH(source);
		if (ArcFormat.DEGREE == arcFormat) {
			target.setX(rad2Deg(target.getX()));
			target.setY(rad2Deg(target.getY()));
			target.setZ(target.getZ());
		}

		return target;
	}
	
	/**
	 * 同一参考椭球下，大地坐标(弧度) --> 空间直角坐标
	 * @param source (rad)
	 * @return
	 */
	public static Coordinate BLH2ECEF(Coordinate source) {	
		double a = source.getRe().getA();
		double square_e1 = source.getRe().getSquare_el();
		
		double B = source.getX();
		double L = source.getY();
		double H = source.getZ();
		
		// 卯酉圈曲率半径
		double N = a/Math.sqrt(1-square_e1 * Math.pow(Math.sin(B), 2));
		double X = (N+H) * Math.cos(B) * Math.cos(L);
		double Y = (N+H) * Math.cos(B) * Math.sin(L);
		double Z = (N*(1-square_e1)+H) * Math.sin(B);
		
		Coordinate target = null;	// ?
		target = new Coordinate();
		target.setRe(source.getRe());
		target.setCoordType(CoordinateType.SPATIAL_RECTANGULAR);
		target.setX(X);
		target.setY(Y);
		target.setZ(Z);
		
		return target;		
	}
	
	/**
	 * 大地坐标--> 高斯平面坐标
	 * @param B:	纬度
	 * @param L:	经度
	 * 		  L0:	中央子午线
	 * 		  addY:	Y轴加常数
	 * @return: x,y
	 */
	public static Coordinate BL2xy(Ellipsoid re,double B,double L,double L0,double addY) {
		
		double a = re.getA();	// 椭球长半轴
		double e1 = re.getE1();
		double e2 = re.getE2();
		
		
		double deltaL = L - L0;	
		double N = a/Math.sqrt(1-Math.pow(e1*Math.sin(B), 2));
		double t = Math.tan(B);
		double eta = e2 * Math.cos(B);
		
		double C1 = 1.0 + (3.0/4) * Math.pow(e1, 2) + (45.0/64) * Math.pow(e1, 4) 
				+ (175.0/256) * Math.pow(e1, 6) + (11025.0/16384) * Math.pow(e1, 8);
		
		double C2 = (3.0/4) * Math.pow(e1, 2) + (15.0/16) * Math.pow(e1, 4) 
				+ (525.0/512) * Math.pow(e1, 6) + (2205.0/2048) * Math.pow(e1, 8);		
		
		double C3 = (15.0/64) * Math.pow(e1, 4) + (105.0/256) * Math.pow(e1, 6) 
				+ (2205.0/4096) * Math.pow(e1, 8);
		
		double C4 = (35.0/512) * Math.pow(e1, 6) + (315.0/2048) * Math.pow(e1, 8);
		
		double C5 = (315.0/131072) * Math.pow(e1, 8);
		
		double X = a * (1-e1*e1)*(C1*B - (1.0/2)*C2*Math.sin(2*B) + (1.0/4)*C3*Math.sin(4*B) 
				- (1.0/6)*C4*Math.sin(6*B) + C5*Math.sin(8*B));
		
		double x = X + (1.0/2)*N* Math.sin(B) * Math.cos(B) * Math.pow(deltaL, 2) 
				*(1 + (1.0/12) * Math.pow(deltaL*Math.cos(B),2) * (5 - t*t + 9*eta*eta + 4*Math.pow(eta, 4)) 
						+ (1.0/360)*Math.pow(deltaL*Math.cos(B), 4)*(61-58*t*t + Math.pow(t, 4)));
		
		double y = addY + N*Math.cos(B)*deltaL*(1 + (1.0/6)*Math.pow(deltaL*Math.cos(B),2)*(1-t*t +eta*eta) 
				+ (1.0/120)*Math.pow(deltaL*Math.cos(B), 4)*(5 - 18*t*t +Math.pow(t, 4) - 14*eta*eta - 58*Math.pow(eta*t, 2)));
		
		return new Coordinate(re,CoordinateType.GAUSS_PLAN,x,y,0);
	}
	
	/**
	 * 高斯平面坐标 --> 大地坐标（角度单位都是弧度）
	 * @param re
	 * @param x
	 * @param y
	 * @param L0
	 * @param addY
	 * @return
	 */
	public static Coordinate xy2BL(Ellipsoid re,double x,double y,double L0,double yAdd) {
		y = y -yAdd;
		
		double a = re.getA();	// 椭球长半轴
		double e1 = re.getE1();
		double e2 = re.getE2();
		
		double C1 = 1.0 + (3.0/4) * Math.pow(e1, 2) + (45.0/64) * Math.pow(e1, 4) 
				+ (175.0/256) * Math.pow(e1, 6) + (11025.0/16384) * Math.pow(e1, 8)
				+ (43659.0/65536)*Math.pow(e1, 10);
		
		double C2 = (3.0/4) * Math.pow(e1, 2) + (15.0/16) * Math.pow(e1, 4) 
				+ (525.0/512) * Math.pow(e1, 6) + (2205.0/2048) * Math.pow(e1, 8)
				+ (72765.0/65536)*Math.pow(e1, 10);		
		
		double C3 = (15.0/64) * Math.pow(e1, 4) + (105.0/256) * Math.pow(e1, 6) 
				+ (2205.0/4096) * Math.pow(e1, 8) + (10395.0/16384)*Math.pow(e1, 10);
		
		double C4 = (35.0/512) * Math.pow(e1, 6) + (315.0/2048) * Math.pow(e1, 8)
				+ (31185.0/131072)*Math.pow(e1, 10);
		
		double C5 = (315.0/16384) * Math.pow(e1, 8) + (3645.0/65536)*Math.pow(e1, 10);
		
		// 迭代计算Bf,底点纬度（横坐标y在高斯投影带中央子午线上的垂足点的纬度）
		double Bf0 = x / (a*C1*(1-e1*e1));
		double Bf1 = 0.0;
		while(true) {
			Bf1 = x / (a*C1*(1-e1*e1)) + (1.0/C1)*( (1.0/2)*C2*Math.sin(2*Bf0) 
												  - (1.0/4)*C3*Math.sin(4*Bf0)
				                                  + (1.0/6)*C4*Math.sin(6*Bf0)
				                                  - (1.0/8)*C5*Math.sin(8*Bf0));
			if(Math.abs(Bf1-Bf0) < CoordTransform.as2Rad(0.0001)) break;  // less than 0.0001 second.
			else{
				Bf0 = Bf1;
			}
		}
		
		// 底点纬度Bf处的子午圈曲率半径
		double Mf = a*(1-e1*e1)/Math.sqrt(Math.pow(1-Math.pow(e1*Math.sin(Bf1), 2), 3));
		// 底点纬度Bf处的卯酉圈曲率半径
		double Nf = a/Math.sqrt(1-Math.pow(e1*Math.sin(Bf1), 2));
		double t = Math.tan(Bf1);
//		double eta = Math.pow(e2 * Math.cos(Bf1),2);
		double eta = e2 * Math.cos(Bf1);
	
		double B = Bf1 - y*y*t/(2*Mf*Nf) 
				*(1 - y*y*(5 + eta*eta + 3*t*t - 9*eta*eta*t*t)/(12*Nf*Nf) 
				    + Math.pow(y, 4)*(61+90*t*t+45*t*t)/(360*Math.pow(Nf, 4)) );
		
		double L = L0 + y/(Nf*Math.cos(Bf1))
				*(1 - y*y*(1 + eta*eta + 2*t*t)/(6*Nf*Nf) 
				    + Math.pow(y, 4)*(5 + 6*eta*eta + 28*t*t + 8*eta*eta*t*t + 24*Math.pow(t, 4))/(120*Math.pow(Nf, 4)));
		
		return new Coordinate(re,CoordinateType.GEODETIC,B,L,0);
	}
	
	/**
	 * 空间直角坐标--> 站心坐标（NEU）
	 * @param originPoint：站心点的空间直角坐标,(XYZ)
	 * @param targetPoint：目标点的空间直角坐标,(XYZ)
	 * @return: NEU
	 */
	public static Coordinate ECEF2NEU(Coordinate originPoint,Coordinate targetPoint){
		Coordinate originBLH = ECEF2BLH(originPoint);	// 站心点的经纬度坐标（弧度）
		double B = originBLH.getX();
		double L = originBLH.getY();		
		
		// 原点到目标点的基线向量
		Array2DRowRealMatrix baselineO2T = new Array2DRowRealMatrix(new double[] {
				targetPoint.getX()-originPoint.getX(),
				targetPoint.getY()-originPoint.getY(),
				targetPoint.getZ()-originPoint.getZ()});
		
		// 旋转矩阵
		Array2DRowRealMatrix r = new Array2DRowRealMatrix(new double[][] {
			{-Math.sin(B)*Math.cos(L),-Math.sin(B)*Math.sin(L),Math.cos(B)},
			{-Math.sin(L),Math.cos(L),0},
			{Math.cos(B)*Math.cos(L),Math.cos(B)*Math.sin(L),Math.sin(B)}
		});
		
		// 目标点的站心坐标向量
		Array2DRowRealMatrix v = r.multiply(baselineO2T);			
		
		Coordinate NEU = null;	// ？
		NEU = new Coordinate(v.getEntry(0, 0),v.getEntry(1, 0),v.getEntry(2, 0));
		
		return NEU;
	}
	
	/**
	 * 根据空间两个点的空间直角坐标，计算目标点的方位角(弧度)
	 * @param NEU
	 * @return
	 */
	public static double calculateAzimuth(Coordinate originPoint,Coordinate targetPoint) {
		Coordinate NEU = ECEF2NEU(originPoint, targetPoint);
		return calculateAzimuth(NEU);
	}
	
	/**
	 * 根据空间两个点的空间直角坐标，计算目标点的高度角(弧度)
	 * @param NEU
	 * @return
	 */
	public static double calculateElevation(Coordinate originPoint,Coordinate targetPoint) {
		Coordinate NEU = ECEF2NEU(originPoint, targetPoint);
		return calculateElevation(NEU);
	}
	
	/*
	 *  以站心坐标为原点，计算卫星的方位角(弧度)
	 */
	public static double calculateAzimuth(Coordinate NEU) {
		return Math.atan2(NEU.getY(), NEU.getX());
	}
	
	/*
	 *  以站心坐标为原点，计算卫星的高度角(弧度)
	 */
	public static double calculateElevation(Coordinate NEU) {
	//	return Math.atan( Math.sqrt(NEU.getX()*NEU.getX() + NEU.getY()*NEU.getY()));
		
		double a = calculateAzimuth(NEU);
		return Math.atan(NEU.getZ()/(NEU.getX()*Math.cos(a) + NEU.getY()*Math.sin(a)));
	}
	

	/**
	 * 弧度 --> 到度
	 * 
	 * @param rad
	 * @return
	 */
	public static double rad2Deg(double rad) {
		return (180 / Math.PI) * rad;
	}

	/**
	 * 度 --> 弧度
	 * 
	 * @param rad
	 * @return
	 */
	public static double deg2Rad(double deg) {
		return (Math.PI / 180) * deg;
	}

	/**
	 * 秒 转换为 弧度,arc second to rad
	 * 
	 * @return
	 */
	public static double as2Rad(double as) {
		return deg2Rad(as / 3600);
	}

}
