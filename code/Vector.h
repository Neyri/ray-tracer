#pragma once

class Vector {
public:
	Vector(double x = 0, double y = 0, double z = 0) :x(x), y(y), z(z) {};
	double norm2() { return x * x + y * y + z * z; }
	void normalize() {
		double n = sqrt(norm2());
		x /= n;
		y /= n;
		z /= n;
	}
	Vector& operator+=(const Vector& b) {
		x += b.x;
		y += b.y;
		z += b.z;
		return *this;
	}
	double x, y, z;
};

Vector operator+(const Vector& a, const Vector &b) {
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}
Vector operator-(const Vector& a, const Vector &b) {
	return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}
Vector operator-(const Vector& a) {
	return Vector(-a.x, -a.y, -a.z);
}
Vector operator*(const double& a, const Vector &b) {
	return Vector(a * b.x, a * b.y, a * b.z);
}
Vector operator*(const Vector& a, const double &b) {
	return b * a;
}
Vector operator*(const Vector& a, const Vector &b) {
	return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
}
Vector operator/(const Vector& a, const double &b) {
	return Vector(a.x / b, a.y / b, a.z / b);
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
	);
}
double dot(const Vector& a, const Vector& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
std::ostream &operator<<(std::ostream &os, Vector v) {
	std::stringstream ss;
	ss << "Vector(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os << ss.str();
}