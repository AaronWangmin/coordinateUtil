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
	 * 空间直角坐标--> 站心坐标（NEU）
	 * @param originPoint：站心点的空间直角坐标,(XYZ)
	 * @param targetPoint：目标点的空间直角坐标,(XYZ)
	 * @return: NEU
	 */
	public static Coordinate calculateNEU(Coordinate originPoint,Coordinate targetPoint){
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
		Coordinate NEU = calculateNEU(originPoint, targetPoint);
		return calculateAzimuth(NEU);
	}
	
	/**
	 * 根据空间两个点的空间直角坐标，计算目标点的高度角(弧度)
	 * @param NEU
	 * @return
	 */
	public static double calculateElevation(Coordinate originPoint,Coordinate targetPoint) {
		Coordinate NEU = calculateNEU(originPoint, targetPoint);
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
		return Math.atan( Math.sqrt(NEU.getX()*NEU.getX() + NEU.getY()*NEU.getY())/NEU.getZ());
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
