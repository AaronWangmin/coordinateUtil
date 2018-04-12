package com.shbd.ddap.coordinate;

import com.shbd.ddap.coordinate.ConstantHolder.CoordinateType;

public class Coordinate {
	private Ellipsoid re;
	private CoordinateType coordType;
	private double x;
	private double y;
	private double z;

	public Coordinate() {}

	// 默认椭球：WGS84椭球,默认坐标形式：空间直角坐标
	public Coordinate(double x, double y, double z) {
		super();
		this.re = ConstantHolder.ReferenceEllipsoid.WGS84;
//		this.coordType = CoordinateType.SPATIAL_RECTANGULAR;
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public Coordinate(Ellipsoid re, CoordinateType coordType, double x, double y, double z) {
		super();
		this.re = re;
		this.coordType = coordType;
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public Ellipsoid getRe() {
		return re;
	}

	public void setRe(Ellipsoid re) {
		this.re = re;
	}

	public CoordinateType getCoordType() {
		return coordType;
	}

	public void setCoordType(CoordinateType coordType) {
		this.coordType = coordType;
	}

	public double getX() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}

	public double getY() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}

	public double getZ() {
		return z;
	}

	public void setZ(double z) {
		this.z = z;
	}

	@Override
	public String toString() {
		return "Coordinate [referenceEllipsoid=" + re.getType() + ", coordType=" + coordType + ", x=" + x + ", y=" + y + ", z=" + z + "]";
	}

	
}
