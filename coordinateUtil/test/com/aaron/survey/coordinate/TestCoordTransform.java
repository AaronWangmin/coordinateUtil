package com.aaron.survey.coordinate;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.function.Abs;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Assert;
import org.junit.Test;

import com.aaron.survey.ConstantHolder;
import com.aaron.survey.adjust.ParamAdjust;
import com.aaron.survey.angle.AngleUtil;
import com.aaron.survey.coordinate.CoordinateUtil;
import com.aaron.survey.coordinate.Coordinate;
import com.aaron.survey.coordinate.Ellipsoid;

public class TestCoordTransform {

	@Test
	public void testBL2xy() {
		
		Ellipsoid re = ConstantHolder.ReferenceEllipsoid.BJ54;
		double B = AngleUtil.deg2Rad(30.8692658056);
		double L = AngleUtil.deg2Rad(121.8836248944);
		double L0 = AngleUtil.deg2Rad(120);
		double yAdd = 500000;
		
		Coordinate result = CoordinateUtil.BL2xy(re, B, L, L0,yAdd);
		
		assert (result.getX()-3418060.509 <= 0.05) && (result.getY()-680157.668 <= 0.05);
		
		System.out.println(result);	
		
		
	}
	
	@Test
	public void testxy2BL() {
		
		Ellipsoid re = ConstantHolder.ReferenceEllipsoid.BJ54;
		double x = 3418060.509;
		double y = 680157.668;
		double L0 = AngleUtil.deg2Rad(120);
		double yAdd = 500000;
		
		Coordinate result = CoordinateUtil.xy2BL(re, x, y, L0, yAdd);
		
		double degB = AngleUtil.rad2Deg(result.getX());
		double degL = AngleUtil.rad2Deg(result.getY());
		
		assert (degB-30.8692658056 <= AngleUtil.as2Rad(0.0001)) && (degL-121.8836248944 <= AngleUtil.as2Rad(0.0001));
		
		System.out.println(degB + " , " + degL);
		System.out.println(result);
		
	}
	
	@Test
	public void testCalculateServenParam() {
		
		ArrayList<Coordinate> sources = new ArrayList<Coordinate>();
		Coordinate source1 = new Coordinate(-2085738.7757, 5503702.8697, 2892977.6829);
		Coordinate source2 = new Coordinate(-2071267.5135, 5520926.7235, 2883341.8135);
		Coordinate source4 = new Coordinate(-2093693.1744, 5511218.2651, 2869861.8947);
		Coordinate source5 = new Coordinate(-2113681.5062, 5491864.0382, 2896934.4852);
		sources.add(source1);
		sources.add(source2);
		sources.add(source4);
		sources.add(source5);
		
		ArrayList<Coordinate> targets = new ArrayList<Coordinate>();
		Coordinate target1 = new Coordinate(-2085635.1879, 5503757.4154, 2892982.0896);
		Coordinate target2 = new Coordinate(-2071164.1636, 5520981.4653, 2883346.1670);
		Coordinate target4 = new Coordinate(-2093589.3723, 5511272.3144, 2869866.0221);
		Coordinate target5 = new Coordinate(-2113577.7476, 5491917.9895, 2896938.5457);
		targets.add(target1);
		targets.add(target2);
		targets.add(target4);
		targets.add(target5);
		
		ParamAdjust pa = CoordinateUtil.calculateServenParam(sources, targets);
		
		// 获取程序计算所得7参数
		RealVector X = pa.calculateX();
		
		// 将三个旋转参数的单位设置为：秒		
		X.setSubVector(4, 
				new ArrayRealVector(new double[] {AngleUtil.rad2As(X.getEntry(4)),
												  AngleUtil.rad2As(X.getEntry(5)),
												  AngleUtil.rad2As(X.getEntry(6)),}));
		// 7参数的对比值（已知值）
		ArrayRealVector ref_X= new ArrayRealVector(
				new double[] {280.083679038385,57.495009802255,116.647652203164,
							  0.6304723858*Math.pow(10, -6),
							  2.938211271759,3.529785370116,-4.710281764406});	// 秒
		// 若最大值小于 0.001，则成功
		RealVector diff = (X.add(ref_X.mapMultiplyToSelf(-1))).map(new Abs());
		assert (diff.getMaxValue() < Math.pow(10, -3));
		
		System.out.println("diff----------------");
		System.out.println(diff.toString());
		pa.showAdjustResult();
		
		System.out.println(pa.calculateX().toString());
	}
	
	@Test
	public void testECEF2ECEFByServenParam() {
		Coordinate source = new Coordinate(-2079412.5535, 5512450.8800, 2879771.2119);
		
		ArrayRealVector servenParam = new ArrayRealVector(
				new double[] {280.0836790372,57.4950097981,116.6476522069,
							  0.0000006305,
							  0.0000142449,0.0000171129,-0.0000228361});
		
		Coordinate actual = CoordinateUtil.ECEF2ECEFByServenParam(source,
				servenParam,ConstantHolder.ReferenceEllipsoid.BJ54);
		
		System.out.println(actual.toString());
		
		Coordinate expected = new Coordinate(-2079308.984,5512505.3689,2879775.4919);
		System.out.println("X坐标误差: " + (actual.getX() - expected.getX()));
		System.out.println("Y坐标误差: " + (actual.getY() - expected.getY()));
		System.out.println("Z坐标误差: " + (actual.getZ() - expected.getZ()));
		
		Assert.assertTrue("X坐标误差大于0.1米", Math.abs(actual.getX() - expected.getX())<0.1);
		Assert.assertTrue("Y坐标误差大于0.1米", Math.abs(actual.getY() - expected.getY())<0.1);
		Assert.assertTrue("Z坐标误差大于0.1米", Math.abs(actual.getZ() - expected.getZ())<0.1);
	}
	
	@Test
	public void testCalculateFourParam() {
		ArrayList<Coordinate> sources = new ArrayList<Coordinate>();
		Coordinate source1 = new Coordinate(3348724.128, 526679.306, 0);
		Coordinate source2 = new Coordinate(3324702.854, 522575.273, 0);
		sources.add(source1);
		sources.add(source2);
		
		ArrayList<Coordinate> targets = new ArrayList<Coordinate>();
		Coordinate target1 = new Coordinate(81281.89, 91322.285, 0);
		Coordinate target2 = new Coordinate(57266.11, 87186.224, 0);
		targets.add(target1);
		targets.add(target2);
		
		ParamAdjust fourParam = CoordinateUtil.calculateFourParam(sources, targets);
		
		// 获取程序计算所得4参数
				RealVector X = fourParam.calculateX();
				
		// 将一个旋转参数的单位设置为：秒		
		X.setEntry(3, AngleUtil.rad2As(X.getEntry(3)));
		
		X.setEntry(2, 1+X.getEntry(2));
		
		// 4参数的对比值（已知值）
		ArrayRealVector ref_X= new ArrayRealVector(
				new double[] {-3266736.944339,-439821.976053,
							  0.999999998447,
						      AngleUtil.rad2As(0.0013334707)});	// 秒
		
		// 若最大值小于 0.001，则成功
		RealVector diff = (X.add(ref_X.mapMultiplyToSelf(-1))).map(new Abs());
		assert (diff.getMaxValue() < Math.pow(10, -3));
		
		System.out.println("/n---------------------4参数----------------------");
		System.out.println("diff----------------");
		System.out.println(diff.toString());
		
		fourParam.showAdjustResult();
	}
	
}
