package com.shbd.ddap.coordinate;


import org.junit.Test;

import com.shbd.ddap.coordinate.ConstantHolder.ArcFormat;

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
	public void testRad2Deg() {
		double X = -2888608.4197;
		double Y = 4655861.6904;
		double theta = Math.atan2(Y, X);
		
		System.out.println(CoordTransform.rad2Deg(theta));
	}

	@Test
	public void testCoordTransform() {		
		
		Coordinate c1 = new Coordinate(-2888608.4197,4655861.6904,3254028.1766);
		System.out.println(CoordTransform.ECEF2BLH(c1));		
		
		System.out.println(CoordTransform.ECEF2BLH(c1,ArcFormat.DEGREE));
		
	}

}
