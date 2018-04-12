package com.aaron.survey.adjust;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.aaron.survey.matrix.MatrixUtil;

/**
 *  参数平差 
 *  误差方程： V = B X - L 
 * 最小二乘解：X = Inverse(B'P B) B' P L
 * 单位权方差：sigma =( V'P V /(n-t))^0.5
 *  协因数阵： Qx = Inverse(B'PB) 
 * 
 * 
 * @author Administrator
 *
 */
public class ParamAdjust {
	public Array2DRowRealMatrix B;
	public Array2DRowRealMatrix P;
	public ArrayRealVector L;
	
	public ParamAdjust() {
	}

	public ParamAdjust(Array2DRowRealMatrix B, Array2DRowRealMatrix P, ArrayRealVector L) {
		this.B = B;
		this.P = P;
		this.L = L;
	}
	
	public RealVector calculateX() {
		RealMatrix N = this.calculateN();
		RealVector W = this.calculateW();
		return MatrixUtil.inverseMatrix(N).operate(W);
	}
	
	public RealMatrix calculateQx() {
		RealMatrix N = this.calculateN();
		return MatrixUtil.inverseMatrix(N);
	}
	
	public double calculateSigma() {
		int n = B.getRowDimension();
		int t = B.getColumnDimension();
		
		RealMatrix V = MatrixUtil.vector2Matrix(this.calculateV());
		RealMatrix VTPV = V.transpose().multiply(P).multiply(V);
		return Math.sqrt(VTPV.getEntry(0, 0)/(n-t));
	}
	
	public RealVector calculateV() {
		RealVector X = this.calculateX();
		return B.operate(X).add(L.mapMultiplyToSelf(-1));
	}
	
	public void showAdjustResult() {
		System.out.println("****************已知条件*****************");
		System.out.println("n = " + this.getB().getRowDimension());
		System.out.println("t = " + this.getB().getColumnDimension());
		
		System.out.println("------------------B:---------------------");
		System.out.println(this.getB().toString());
		System.out.println("------------------L:---------------------");
		System.out.println(this.getL().toString());
		System.out.println("------------------P:---------------------");
		System.out.println(this.getP().toString());
		
		System.out.println("****************计算结果*****************");
		System.out.println("------------------X:---------------------");
		System.out.println(this.calculateX().toString());
		System.out.println("------------------sigma:---------------------");
		System.out.println(this.calculateSigma());
		System.out.println("------------------Qx:---------------------");
		System.out.println(this.calculateQx().toString());
		
	}
	
	
	private RealMatrix calculateN() {
		return B.transpose().multiply(P).multiply(B);
	}
	
	private RealVector calculateW() {
		return B.transpose().multiply(P).operate(L);
	}

	public Array2DRowRealMatrix getB() {
		return B;
	}	

	public void setB(Array2DRowRealMatrix b) {
		B = b;
	}

	public Array2DRowRealMatrix getP() {
		return P;
	}

	public void setP(Array2DRowRealMatrix p) {
		P = p;
	}

	public ArrayRealVector getL() {
		return L;
	}

	public void setL(ArrayRealVector l) {
		L = l;
	}

}
