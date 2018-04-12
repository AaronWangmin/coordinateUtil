package com.aaron.survey.angle;

public class AngleUtil {
//	private double angleValue;	// 默认为 弧度
//	
//	public AngleUtil() {
//	}
//	
//	// 弧度 --构造
//	public AngleUtil(double angleValue) {
//		this.angleValue = angleValue;
//	}
//	
//	// 度分秒 --构造
//	public AngleUtil(double d,double m,double s) {
//		this.angleValue = dms2Rad(d, m, s);
//	}
//	
//	// 其它单位 --构造
//	public AngleUtil(double angleValue,ConstantHolder.ArcFormat arcFormat) {
//		switch(arcFormat) {
//		case RAD:
//			this.angleValue = angleValue;
//			break;
//		case DEGREE:
//			this.angleValue = deg2Rad(angleValue);
//			break;
//		case DMS:
//			throw new SurveyException("DMS is not use,please choice other units");
////			break;
//		case SECOND:
//			this.angleValue = as2Rad(angleValue);
//			break;
//		}
//	}
//
//	public void setAngleValue(double angleValue) {
//		this.angleValue = angleValue;
//	}
//	
//	public double getRad() {
//		return angleValue;
//	}
//	
//	public double getDegree() {
//		return rad2Deg(angleValue);
//	}
//	
//	public double getSecond() {
//		return rad2As(angleValue);
//	}
//	
//	public double[] getDMS() {
//		return rad2Dms(angleValue);
//	}
//	
//	public String getDmsText() {
//		double[] dms = rad2Dms(angleValue);
//		return dms[0] + " : " + dms[1] + " : " + dms[2] ;
//	}
//	
//	public void showAngle() {
//		System.out.println(this.getDmsText());
//		System.out.println(this.getDegree());
//		System.out.println(this.getSecond());
////		System.out.println(this.getDMS());
//		System.out.println(this.getRad());
//	}
	
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
