package com.aaron.survey.angle;


import org.junit.Test;

import com.aaron.survey.ConstantHolder.ArcFormat;

public class TestAngleUtil {
	Angle a1 = new Angle(Math.PI/2.0);
	
	Angle a2 = new Angle(30,50,12.009);
	
	Angle a3 = new Angle(122, ArcFormat.DEGREE);
	
	Angle a4 = new Angle(61, ArcFormat.SECOND);
	
	Angle a5 = new Angle(Math.PI/2.0, ArcFormat.RAD);

	@Test
	public void rad2DmsTest() {
		
		a1.showAngle();
		
		a2.showAngle();
		
		a3.showAngle();
		
		a4.showAngle();
		
		a5.showAngle();
	}

}
