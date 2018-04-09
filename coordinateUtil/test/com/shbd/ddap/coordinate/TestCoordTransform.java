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
	
	

}
