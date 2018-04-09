package com.shbd.ddap.coordinate;

public final class ConstantHolder {
	
	private ConstantHolder() {}
	
	// 参考椭球类型：
	public enum ReferenceEllipsoidType{
		WGS84,CGCS2000,BJ54,XA80
	}
	
	// 坐标系统类型：
//	public enum CoordinateSystemType{
//		WGS84,CGCS2000,BJ54,XA80
//	}
	
	// 坐标类型：空间直角坐标、大地坐标、高斯平面坐标
	public enum CoordinateType{
		SPATIAL_RECTANGULAR,GEODETIC,GAUSS_PLAN
	}
	
	// 角度格式：弧度、度、度/分/秒
	public enum ArcFormat{
		RAD,DEGREE,DMS
	}
	
	// 参考椭球
	public static class ReferenceEllipsoid {
		// WGS84
		public static Ellipsoid WGS84 = new Ellipsoid(ReferenceEllipsoidType.WGS84,
				6378137.0,
				1/298.257223563,
				3.986005E14,
				7.292115E-5);
		
		// CGCS2000
		public static Ellipsoid CGCS2000 = new Ellipsoid(ReferenceEllipsoidType.CGCS2000,
				6378137.0,
				1/298.257222101,
				3.986004418E14,
				7.292115E-5);
		
		// BJ54
		public static Ellipsoid BJ54 = new Ellipsoid(ReferenceEllipsoidType.BJ54,
				6378245.0,
				1/298.3,
				0,		//未知
				0);		//未知
	}	

}
