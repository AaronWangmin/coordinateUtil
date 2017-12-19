package com.shbd.ddap.coordinate;


import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.junit.Test;

import com.shbd.ddap.coordinate.ConstantHolder.ArcFormat;
import com.shbd.ddap.coordinate.ConstantHolder.CoordinateType;

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
		
		System.out.println(CoordTransform.rad2Deg(theta));
	}

	@Test
	public void testECEF2BLH() {		
		
		Coordinate c1 = new Coordinate(-2888608.4197,4655861.6904,3254028.1766);
		System.out.println(c1);
		System.out.println(CoordTransform.ECEF2BLH(c1));		
		
		System.out.println(CoordTransform.ECEF2BLH(c1,ArcFormat.DEGREE));
		
	}
	
	@Test
	public void testBLH2ECEF() {
		double B = CoordTransform.deg2Rad(30.874999463688056);
		double L = CoordTransform.deg2Rad(121.816458229654);
		double H = 42.057353113777936;
		
		Coordinate c1 = new Coordinate(B,L,H);
		System.out.println(c1);
		System.out.println(CoordTransform.BLH2ECEF(c1));		
	}
	
	@Test
	public void testCalculateNEU() {
		Coordinate originPoint = new Coordinate(1,2,3);
		Coordinate targetPoint = new Coordinate(3,2,1);
		Coordinate NEU = CoordTransform.calculateNEU(originPoint,targetPoint);
		System.out.println(NEU);
		
		double azimuth = CoordTransform.calculateAzimuth(NEU);
		double elevation = CoordTransform.calculateElevation(NEU);
		
		System.out.println(CoordTransform.rad2Deg(azimuth));
		
		System.out.println(CoordTransform.rad2Deg(elevation));
		
	}
	
	

}
