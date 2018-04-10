package com.shbd.ddap.coordinate;

import org.junit.Test;

public class TestCoordTransform {

	@Test
	public void testBL2xy() {
		
		Ellipsoid re = ConstantHolder.ReferenceEllipsoid.BJ54;
		double B = CoordTransform.deg2Rad(30.8692658056);
		double L = CoordTransform.deg2Rad(121.8836248944);
		double L0 = CoordTransform.deg2Rad(120);
		double yAdd = 500000;
		
		Coordinate result = CoordTransform.BL2xy(re, B, L, L0,yAdd);
		
		assert (result.getX()-3418060.509 <= 0.05) && (result.getY()-680157.668 <= 0.05);
		
		System.out.println(result);	
		
		
	}
	
	@Test
	public void testxy2BL() {
		
		Ellipsoid re = ConstantHolder.ReferenceEllipsoid.BJ54;
		double x = 3418060.509;
		double y = 680157.668;
		double L0 = CoordTransform.deg2Rad(120);
		double yAdd = 500000;
		
		Coordinate result = CoordTransform.xy2BL(re, x, y, L0, yAdd);
		
		double degB = CoordTransform.rad2Deg(result.getX());
		double degL = CoordTransform.rad2Deg(result.getY());
		
		assert (degB-30.8692658056 <= CoordTransform.as2Rad(0.0001)) && (degL-121.8836248944 <= CoordTransform.as2Rad(0.0001));
		
		System.out.println(degB + " , " + degL);
		
	}

}
