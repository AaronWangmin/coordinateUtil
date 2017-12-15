package com.shbd.ddap.coordinate;

import com.shbd.ddap.coordinate.ConstantHolder.ArcFormat;
import com.shbd.ddap.coordinate.ConstantHolder.CoordinateType;

public class CoordTransform {

	/**
	 *  同一参考椭球下，空间直角坐标--> 大地坐标
	 * @param source
	 * @return 大地坐标（弧度）
	 */
	public static Coordinate ECEF2BLH(Coordinate source) {
		// 如果源坐标类型不是地心地固坐标（XYZ）时，函数中止退出。
		if (source.getCoordType() != CoordinateType.SPATIAL_RECTANGULAR) {
			return null;
		}

		Coordinate target = new Coordinate();
		target.setRe(source.getRe()); // 设置转换前后的椭球相同
		target.setCoordType(CoordinateType.GEODETIC); // 设置转换后的坐标为大地坐标(B/L/H)

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
