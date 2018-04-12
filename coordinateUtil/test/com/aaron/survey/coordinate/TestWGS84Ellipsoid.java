package com.aaron.survey.coordinate;


import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.junit.Test;

import com.aaron.survey.ConstantHolder.ArcFormat;
import com.aaron.survey.ConstantHolder.CoordinateType;
import com.aaron.survey.coordinate.CoordinateUtil;
import com.aaron.survey.coordinate.Coordinate;

public class TestWGS84Ellipsoid {
	
	@Test
	public void testMath() {
		System.out.println(Math.sqrt(4));
		System.out.println(Math.pow(3, 2));		
		
		double X = -2830150.6998;
		double Y = 4651464.6082;
		double theta = Math.atan2(Y, X);
		
		System.out.println(theta);
	}
	
	@Test
	public void testMatrix() {
		Array2DRowRealMatrix r = new Array2DRowRealMatrix(new double[][] {
			{1,2,3},
			{4,5,6},
			{7,8,9},
			{11,12,13}
		});
		
		System.out.println(r);
		System.out.println(r.getColumnDimension());
		System.out.println(r.getRow(0)[0]);
	}
	
	@Test
	public void testRad2Deg() {
		double X = -2888608.4197;
		double Y = 4655861.6904;
		double theta = Math.atan2(Y, X);
		
		System.out.println(CoordinateUtil.rad2Deg(theta));
	}

	@Test
	public void testECEF2BLH() {		
		
		Coordinate c1 = new Coordinate(-2888608.4197,4655861.6904,3254028.1766);
		System.out.println(c1);
		System.out.println(CoordinateUtil.ECEF2BLH(c1));		
		
		System.out.println(CoordinateUtil.ECEF2BLH(c1,ArcFormat.DEGREE));
		
	}
	
	@Test
	public void testBLH2ECEF() {
		
		double B = CoordinateUtil.deg2Rad(31+16.0/60 + 9.411132/3600);
		double L = CoordinateUtil.deg2Rad(121+28.0/60 + 57.543924/3600);
		double H = 31.537;
		
//		double B = CoordTransform.deg2Rad(30.874999463688056);
//		double L = CoordTransform.deg2Rad(121.816458229654);
//		double H = 42.057353113777936;
		
		Coordinate c1 = new Coordinate(B,L,H);
		System.out.println(c1);
		System.out.println(CoordinateUtil.BLH2ECEF(c1));		
	}
	
	@Test
	public void testCalculateNEU() {		
		Coordinate originPoint = CoordinateUtil.BLH2ECEF(new Coordinate(
				CoordinateUtil.deg2Rad(15),
				CoordinateUtil.deg2Rad(105),
				0));
		Coordinate targetPoint = CoordinateUtil.BLH2ECEF(new Coordinate(
				CoordinateUtil.deg2Rad(15),
				CoordinateUtil.deg2Rad(135),
				0));
		Coordinate NEU = CoordinateUtil.ECEF2NEU(originPoint,targetPoint);
		System.out.println(NEU);
		
		double azimuth = CoordinateUtil.calculateAzimuth(NEU);
		double elevation = CoordinateUtil.calculateElevation(NEU);
		
		System.out.println(CoordinateUtil.rad2Deg(azimuth));
		
		System.out.println(CoordinateUtil.rad2Deg(elevation));
		
	}
	
	

}
