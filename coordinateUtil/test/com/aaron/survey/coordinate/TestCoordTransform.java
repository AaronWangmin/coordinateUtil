package com.aaron.survey.coordinate;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.function.Abs;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Test;

import com.aaron.survey.ConstantHolder;
import com.aaron.survey.adjust.ParamAdjust;
import com.aaron.survey.coordinate.CoordTransform;
import com.aaron.survey.coordinate.Coordinate;
import com.aaron.survey.coordinate.Ellipsoid;
import com.aaron.survey.matrix.MatrixUtil;

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
		
		ParamAdjust pa = CoordTransform.calculateServenParam(sources, targets);
		
		// 获取程序计算所得7参数
		RealVector X = pa.calculateX();
		
		// 将三个旋转参数的单位设置为：秒		
		X.setSubVector(4, 
				new ArrayRealVector(new double[] {CoordTransform.rad2As(X.getEntry(4)),
												  CoordTransform.rad2As(X.getEntry(5)),
												  CoordTransform.rad2As(X.getEntry(6)),}));
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
	}
}
