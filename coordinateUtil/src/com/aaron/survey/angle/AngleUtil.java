package com.aaron.survey.angle;

public class AngleUtil {
	/**
	 * 弧度 --> 度
	 * 
	 * @param rad
	 * @return
	 */
	public static double rad2Deg(double rad) {
		return (180 / Math.PI) * rad;
	}
	
	/**
	 * 弧度 --> 秒 , rad to arc second
	 * 
	 * @return
	 */
	public static double rad2As(double rad) {
		return rad2Deg(rad) * 3600.0;
	}
	
	/**
	 * 弧度 --> 度分秒 , rad to dms
	 * 
	 * @return
	 */
	public static double[] rad2Dms(double rad) {
		double d = (int)rad2Deg(rad);					// 度
		double m = (int)((rad2Deg(rad) - d)*60.0);		// 分
		double s = ((rad2Deg(rad) - d)*60.0 - m)*60;	//秒
		
		return new double[]{d,m,s};
	}

	/**
	 * 度 --> 弧度
	 * 
	 * @param rad
	 * @return
	 */
	public static double deg2Rad(double deg) {
		return (Math.PI / 180.0) * deg;
	}

	/**
	 * 秒--> 弧度,arc second to rad
	 * 
	 * @return
	 */
	public static double as2Rad(double as) {
		return deg2Rad(as / 3600.0);
	}
	
	/**
	 * 度分秒--> 弧度,dms to rad
	 * 
	 * @return
	 */
	public static double dms2Rad(double d,double m,double s) {
		
		return deg2Rad(d + m/60.0 + s/3600.0);
	}
	
}
