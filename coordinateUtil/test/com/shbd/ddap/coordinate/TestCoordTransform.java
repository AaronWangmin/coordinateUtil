package com.shbd.ddap.coordinate;

import java.util.ArrayList;

import org.apache.commons.math3.linear.RealVector;
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
		System.out.println(result);
		
	}
	
	@Test
	public void testCalculateServenParam() {
		
		
		ArrayList<Coordinate> sources = new ArrayList<Coordinate>();
		Coordinate source1 = new Coordinate(-2858951.7086, 4660391.5147, 3273543.0373);
		Coordinate source2 = new Coordinate(-2888608.4197, 4655861.6904, 3254028.1766);
		Coordinate source3 = new Coordinate(-2849958.9961, 4688645.1877, 3241078.5745);
		sources.add(source1);
		sources.add(source2);
		sources.add(source3);
		
		ArrayList<Coordinate> targets = new ArrayList<Coordinate>();
		Coordinate target1 = new Coordinate(-17199.2926, 5744.6724, 0);
		Coordinate target2 = new Coordinate(-39909.6566, 33405.4226, 0);
		Coordinate target3 = new Coordinate(-55021.4263, -16673.3169, 0);
		targets.add(target1);
		targets.add(target2);
		targets.add(target3);
		
		RealVector result = CoordTransform.calculateServenParam(sources, targets);
		System.out.println(result.toString());
		
		
	}

}
