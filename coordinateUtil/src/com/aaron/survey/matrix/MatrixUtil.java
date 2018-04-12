package com.shbd.ddap.math;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class MatrixUtil {
	
	// 构造给定大小的单位矩阵
	public static Array2DRowRealMatrix eye(int size) {
		Array2DRowRealMatrix result = new Array2DRowRealMatrix(size,size);
		for(int i=0;i<size;i++) {			
			result.setEntry(i, i, 1.0);
		}
		return result;		
	}
	
	// 计算给定矩阵的逆矩阵
	public static RealMatrix inverseMatrix(RealMatrix A) {
		RealMatrix result = new LUDecomposition(A).getSolver().getInverse();		
		return result;		
	}
	
	// 将一个向量转换为矩阵
	public static RealMatrix vector2Matrix(RealVector v) {
		return new Array2DRowRealMatrix(v.toArray());
	}
}
