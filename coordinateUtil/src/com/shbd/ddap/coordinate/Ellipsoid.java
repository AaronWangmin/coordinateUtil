package com.shbd.ddap.coordinate;

import com.shbd.ddap.coordinate.ConstantHolder.ReferenceEllipsoidType;

public class Ellipsoid {
	private ReferenceEllipsoidType type;
	private double a; // 椭球长半轴,major radius of ellipsoid
	private double f; // 椭球扁率,flattening of ellipsoid
	private double GM; // 地球引力常数(地球引力与地球质量的乘积),m^3 * s^(-2)
	private double W; // 地球自转角速度,rad/s

	// 以下参数可根据 a,f两个参数进行计算
	// the below parameters can be caculated on 'a' and 'f'
	private double b; // 椭球短半轴,minor radius of ellipsoid
	private double e1; // 椭球第一偏心率
	private double e2; // 椭球第二偏心率	
	private double square_el; // 椭球第一偏心率的平方
	private double square_e2; // 椭球第二偏心率的平方

	public Ellipsoid() {
	}

	public Ellipsoid(ReferenceEllipsoidType type, double a, double f, double GM, double W) {
		super();
		this.type = type;
		this.a = a;
		this.f = f;
		this.GM = GM;
		this.W = W;
		this.initParameters(a, f);
	}

	public ReferenceEllipsoidType getType() {
		return type;
	}

	public void setType(ReferenceEllipsoidType type) {
		this.type = type;
	}

	private void initParameters(double a, double f) {
		this.b = a * (1 - f);
		this.square_el = (a * a - b * b) / (a * a);
		this.square_e2 = (a * a - b * b) / (b * b);
		this.e1 = Math.sqrt(this.square_el);
		this.e2 = Math.sqrt(this.square_e2);
	}

	public double getA() {
		return a;
	}

	public void setA(double a) {
		this.a = a;
	}

	public double getF() {
		return f;
	}

	public void setF(double f) {
		this.f = f;
	}

	public double getGM() {
		return GM;
	}

	public void setGM(double gM) {
		GM = gM;
	}

	public double getW() {
		return W;
	}

	public void setW(double w) {
		W = w;
	}

	public double getB() {
		return b;
	}

	public void setB(double b) {
		this.b = b;
	}

	public double getE1() {
		return e1;
	}

	public void setE1(double e1) {
		this.e1 = e1;
	}

	public double getE2() {
		return e2;
	}

	public void setE2(double e2) {
		this.e2 = e2;
	}

	public double getSquare_el() {
		return square_el;
	}

	public void setSquare_el(double square_el) {
		this.square_el = square_el;
	}

	public double getSquare_e2() {
		return square_e2;
	}

	public void setSquare_e2(double square_e2) {
		this.square_e2 = square_e2;
	}

	@Override
	public String toString() {
		return "Ellipsoid [type=" + type + ", a=" + a + ", f=" + f + ", GM=" + GM + ", W=" + W + ", b=" + b + ", e1="
				+ e1 + ", e2=" + e2 + ", square_el=" + square_el + ", square_e2=" + square_e2 + "]";
	}

	



}
