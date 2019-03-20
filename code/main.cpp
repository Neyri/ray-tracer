#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>

using namespace std;
using namespace std::chrono;


std::default_random_engine engine;
std::uniform_real_distribution<double> uniform(0, 1);

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "Vector.h"
//#include "Geometry.h"

#define M_PI 3.141592653589793238

double sqr(double x) {
	return x * x;
}

class Ray {
public:
	Ray(const Vector& C, const Vector& u) : C(C), u(u) {};
	Vector u, C;
};

class Object {
public:
	Object()  {};
	virtual bool intersect(const Ray& r, Vector& P, Vector& N, Vector &color) = 0;
	
	Vector albedo;
	bool mirror;
	bool transparency;
	bool light;
	double I;
};

class Sphere: public Object {
	public:
	Sphere(
		const Vector& O = Vector(0, 0, 0),
		double R = 0,
		const Vector& albedo = Vector(0, 0, 0),
		const bool& mirror = false,
		const bool& transparency = false,
		const bool& light = false
	) : O(O), R(R){
		this->albedo =albedo;
		this->mirror = mirror;
		this->transparency = transparency;
		this->light = light;
	};
	virtual bool intersect(const Ray& r, Vector& P, Vector& N, Vector& color) {
		// résolution de t²|u|² + 2t*(u scalaire C-O) + |C-O|² = R²
		double a = 1;
		double b = 2 * dot(r.u, r.C - O);
		double c = (r.C - O).norm2() - R * R;
		double delta = b * b - 4 * a*c;
		if (delta >= 0) {
			double t1 = (-b - sqrt(delta)) / (2 * a);
			double t2 = (-b + sqrt(delta)) / (2 * a);
			if (t1 > 0) {
				P = t1 * r.u + r.C;
				N = P - O;
				N.normalize();
			}
			else {
				if (t2 > 0) {
					P = t2 * r.u + r.C;
					N = P - O;
					N.normalize();
				}
				else {
					return false;
				}
			}
		}
		color = this->albedo;
		return delta >= 0;
	};
	Vector O;
	double R;
};

class Light {
public:
	Light(const Sphere& S, const double& I) : S(S), I(I) {};
	Sphere S;
	double I;
};

class Triangle : public Object {
public:
	Triangle(const Vector& A,const Vector& B, const Vector& C,
		const Vector& albedo = Vector(0, 0, 0),
		const bool& mirror = false,
		const bool& transparency = false,
		const bool& light = false) : 
		A(A), B(B), C(C) {
		this->albedo = albedo;
		this->mirror = mirror;
		this->transparency = transparency;
		this->light = light;
	};

	virtual bool intersect(const Ray& r, Vector& P, Vector& N, Vector& color) {
		color = this->albedo;
		N = cross(B - A, C - A);
		N.normalize();
		if (dot(r.u, N) > 0) {
			N = -N;
		}
		double t = dot(A - r.C, N) / dot(r.u, N);
		if (t < 0) return false;
		P = r.C + t * r.u;
		double alpha, beta, gamma;
		getBarycentric(P, alpha, beta, gamma);
		if (alpha <= 0 || alpha >= 1 ||
			beta  <= 0 || beta  >= 1 ||
			gamma <= 0 || gamma >= 1) {
			return false;
		}
		return true;
	};

	void getBarycentric(const Vector& P, double& alpha, double& beta, double& gamma){
		// On a intersection si P = a*A + b*B + c*C et a + b + c = 1
		Vector AC = C - A;
		Vector AB = B - A;
		Vector AP = P - A;
		double m11 = AB.norm2();
		double m12 = dot(AB, AC);
		double m22 = AC.norm2();
		double denom = m11 * m22 - m12 * m12;

		double b11 = dot(AP, AB);
		double b21 = dot(AP, AC);
		double detb = b11 * m22 - b21 * m12;
		beta = detb / denom;

		double detg = m11 * b21 - m12 * b11;
		gamma = detg / denom;

		alpha = 1 - beta - gamma;
	}

	Vector getCenter() {
		return (A + B + C) / 3.;
	}

	Vector A;
	Vector B;
	Vector C;
};

std::ostream &operator<<(std::ostream &os, Triangle t) {
	std::stringstream ss;
	ss << t.A << " " << t.B << " " << t.C;
	return os << ss.str();
}

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class BBox {
public:
	BBox() {};
	BBox(const Vector& min, const Vector& max) : vmin(min), vmax(max) {};
	// vmin et vmax sont les vecteur contenant les coordonnées minimales et maximales suivant x, y et z

	bool intersect(const Ray& r) const {
		double tx1 = (vmin.x - r.C.x) / r.u.x;
		double tx2 = (vmax.x - r.C.x) / r.u.x;
		double txmin = std::min(tx1, tx2);
		double txmax = std::max(tx1, tx2);

		double ty1 = (vmin.y - r.C.y) / r.u.y;
		double ty2 = (vmax.y - r.C.y) / r.u.y;
		double tymin = std::min(ty1, ty2);
		double tymax = std::max(ty1, ty2);

		double tz1 = (vmin.z - r.C.z) / r.u.z;
		double tz2 = (vmax.z - r.C.z) / r.u.z;
		double tzmin = std::min(tz1, tz2);
		double tzmax = std::max(tz1, tz2);

		double tmin = std::max(std::max(txmin, tymin), tzmin);
		double tmax = std::min(std::min(txmax, tymax), tzmax);
		if (tmax - tmin > 0) {
			return true;
		}
		return false;
	}

	Vector vmin;
	Vector vmax;
};
std::ostream &operator<<(std::ostream &os, BBox bb) {
	std::stringstream ss;
	ss << "BBox(" << bb.vmin << ", " << bb.vmax << ")";
	return os << ss.str();
}

class BVHNode {
public:
	BVHNode() {};
	BVHNode(const int& debut, const int& fin, const BBox& box):debut(debut), fin(fin), bb(box) {};

	bool intersect(const Ray& r, std::vector<BVHNode*>& leaves) {
		if (bb.intersect(r)) {
			//std::cout << "intersect " << bb << ", " << debut << ", " << fin << endl;
			//if (fg) {
				//std::cout << fg->bb << endl;
			//}
			if (!fg || !fd) {
				
				leaves.push_back(this);
				return true;
			}
			bool found_fg = fg->intersect(r, leaves);
			bool found_fd = fd->intersect(r, leaves);
			return found_fd || found_fg;
		}
		return false;
	}

	BVHNode *fg, *fd;
	BBox bb;
	int debut, fin;
};

std::ostream &operator<<(std::ostream &os, BVHNode node) {
	std::stringstream ss;
	ss << "BVHNode(" << node.bb << ", " << node.debut << ", " << node.fin << ")";
	return os << ss.str();
}

class Geometry : public Object {
public:
	~Geometry() {}
	Geometry(const Vector& albedo = Vector(0, 0, 0),
		const bool& mirror = false,
		const bool& transparency = false,
		const bool& light = false) {
		this->albedo = albedo;
		this->mirror = mirror;
		this->transparency = transparency;
		this->light = light;
	};

	void add_textures(const char *filename) {
		int w, h, c;
		textures.push_back(stbi_load(filename, &w, &h, &c, 3));
		textures_width.push_back(w);
		textures_height.push_back(h);
	}

	void readOBJ(const char* obj, bool load_textures) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.x, &vec.y, &vec.z, &col.x, &col.z, &col.y) == 6) {
					col.x = std::min(1., std::max(0., col.x));
					col.y = std::min(1., std::max(0., col.y));
					col.z = std::min(1., std::max(0., col.z));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				}
				else {
					sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec.x, &vec.z, &vec.y); // as we rotate the girl by changing y and z we also need to change the normals
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec.x, &vec.y);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				}
				else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);
	}

	BBox initBox(int i0, int i1) {
		Vector min(1e99, 1e99, 1e99), max(-1e99, -1e99, -1e99);
		for (int i = i0; i < i1; i++) {
			min.x = std::min(min.x, vertices[indices[i].vtxi].x);
			min.x = std::min(min.x, vertices[indices[i].vtxj].x);
			min.x = std::min(min.x, vertices[indices[i].vtxk].x);

			min.y = std::min(min.y, vertices[indices[i].vtxi].y);
			min.y = std::min(min.y, vertices[indices[i].vtxj].y);
			min.y = std::min(min.y, vertices[indices[i].vtxk].y);

			min.z = std::min(min.z, vertices[indices[i].vtxi].z);
			min.z = std::min(min.z, vertices[indices[i].vtxj].z);
			min.z = std::min(min.z, vertices[indices[i].vtxk].z);

			max.x = std::max(max.x, vertices[indices[i].vtxi].x);
			max.x = std::max(max.x, vertices[indices[i].vtxj].x);
			max.x = std::max(max.x, vertices[indices[i].vtxk].x);

			max.y = std::max(max.y, vertices[indices[i].vtxi].y);
			max.y = std::max(max.y, vertices[indices[i].vtxj].y);
			max.y = std::max(max.y, vertices[indices[i].vtxk].y);

			max.z = std::max(max.z, vertices[indices[i].vtxi].z);
			max.z = std::max(max.z, vertices[indices[i].vtxj].z);
			max.z = std::max(max.z, vertices[indices[i].vtxk].z);
		}
		return BBox(min, max);
	}

	void buildTransform(const Vector& translation, const double& scale, const double& rz){
		transform_matrix.resize(3 * 4);
		transform_matrix[0] = sin(rz) * scale;
		transform_matrix[1] = cos(rz) * scale;
		transform_matrix[2] = 0;
		transform_matrix[3] = translation.x;
		transform_matrix[4] = -cos(rz) * scale;
		transform_matrix[5] = sin(rz) * scale;
		transform_matrix[6] = 0;
		transform_matrix[7] = translation.y;
		transform_matrix[8] = 0;
		transform_matrix[9] = 0;
		transform_matrix[10] = scale;
		transform_matrix[11] = translation.z;

		transform_inverse.resize(3 * 4);
		transform_inverse[0] =  sin(rz) / scale;
		transform_inverse[1] = -cos(rz) / scale;
		transform_inverse[2] = 0;
		transform_inverse[3] = (sin(rz) * translation.x + cos(rz) * translation.y) / scale;
		transform_inverse[4] = cos(rz) / scale;
		transform_inverse[5] = sin(rz) / scale;
		transform_inverse[6] = 0;
		transform_inverse[7] = (sin(rz) * translation.x - cos(rz) * translation.y) / scale;
		transform_inverse[8] = 0;
		transform_inverse[9] = 0;
		transform_inverse[10] = 1 / scale;
		transform_inverse[11] = (translation.z) / scale;
	}

	void transform(const Vector& translate = Vector(0, 0, 0), const double scale = 1., const double echange = true) {
		for (int i = 0; i < vertices.size(); i++) {
			if (echange) {
				vertices[i] = Vector(vertices[i].x*scale + translate.x, vertices[i].z*scale + translate.y, vertices[i].y*scale + translate.z);
			}
			else {
				vertices[i] = vertices[i] * scale + translate;
			}
		}
	}

	void initBVH() {
		BBox box = initBox(0, indices.size());
		root = BVHNode(0, indices.size(), box);
		buildBVH(root);
	}

	virtual bool intersect(const Ray& r, Vector& P, Vector& N, Vector& color) {
		int debut, fin;
		std::vector<BVHNode*> leaves;
		if (!root.intersect(r, leaves)) {
			return false;
		}
		else {
			bool found = false;
			double dist_min = 1E99;
			int debut, fin;
			for (int j = 0; j < leaves.size(); j++) {
				debut = leaves[j]->debut;
				fin = leaves[j]->fin;
				for (int i = debut; i < fin; i++) {
					Triangle t(vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk]);
					Vector locP, locN, locColor;
					if (t.intersect(r, locP, locN, locColor)) {
						double dist = (locP - r.C).norm2();
						found = true;
						if (dist < dist_min) {
							dist_min = dist;
							P = locP;
							// Get the normal to the triangle (according to the geometry)
							double alpha, beta, gamma;
							t.getBarycentric(P, alpha, beta, gamma);
							N = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
							N.normalize();
							//N = locN;
							if (dot(N, r.u) > 0) N = -N;
							// Get the color at this point
							Vector uv = alpha * uvs[indices[i].uvi] + beta * uvs[indices[i].uvj] + gamma * uvs[indices[i].uvk];
							int width  = textures_width[indices[i].group];
							int height = textures_height[indices[i].group];
							int x = fabs(fmod(uv.x, 1.)) * width;
							int y = height - fabs(fmod(uv.y, 1.)) * height - 1;
							color.x = textures[indices[i].group][(x + y * width) * 3    ] / 255.;
							color.y = textures[indices[i].group][(x + y * width) * 3 + 1] / 255.;
							color.z = textures[indices[i].group][(x + y * width) * 3 + 2] / 255.;
							//std::cout << color << endl;
						}

					}
				}
			}
			return found;
		}
	}

	void buildBVH(BVHNode& node) {
		BBox nodeBox = node.bb;
		Vector diagBox = nodeBox.vmax - nodeBox.vmin;
		Vector N;
		double criterion;
		// calculate where to divide in 2
		if (diagBox.x > diagBox.y && diagBox.x > diagBox.z) {
			N = Vector(1, 0, 0);
			criterion = (nodeBox.vmax.x + nodeBox.vmin.x) / 2;
		}
		else {
			if (diagBox.y > diagBox.x && diagBox.y > diagBox.z) {
				N = Vector(0, 1, 0);
				criterion = (nodeBox.vmax.y + nodeBox.vmin.y) / 2;
			}
			else {
				N = Vector(0, 0, 1);
				criterion = (nodeBox.vmax.z + nodeBox.vmin.z) / 2;
			}
		}
		// order the indices
		int p = node.debut;
		Vector tCenter;
		for (int i = node.debut; i < node.fin; i++) {
			Triangle t(vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk]);
			tCenter = t.getCenter();
			if (dot(N, tCenter) < criterion) { // à gauche
				std::swap(indices[p], indices[i]);
				p++;
			}
		}
		// recursive
		if (p != node.fin && p - node.debut > 3) {
			BBox box = initBox(node.debut, p);
			//std::cout << node.debut << p << endl;
			//BVHNode fg(node.debut, p, box);
			node.fg = new BVHNode(node.debut, p, box);
			//std::cout << node.fg.bb << endl;
			buildBVH(*node.fg);
		}
		if (p != node.debut && node.fin - p > 3) {
			BBox box = initBox(p, node.fin);
			//std::cout << node.fin << p << endl;
			//BVHNode fd(p, node.fin, box);
			node.fd = new BVHNode(p, node.fin, box);
			buildBVH(*node.fd);
		}
	}

	std::vector<TriangleIndices> indices; // liste des indices des points d'un triangle
	std::vector<Vector> vertices; // liste des points
	std::vector<Vector> normals; // normales aux points
	std::vector<Vector> uvs; // liste des uvs
	std::vector<Vector> vertexcolors;
	std::vector<unsigned char*> textures;
	std::vector<int> textures_width;
	std::vector<int> textures_height;
	std::vector<double> transform_matrix, transform_inverse;
	BVHNode root;
};

class Scene {
public:
	Scene(const Sphere& light) : light(light) {};
	void addObject(Object* s) {
		spheres.push_back(s);
	}
	bool intersect(const Ray& r, Vector& P, Vector& N, Object* &S, Vector& color) {
		double dist_min = 1E15;
		bool found = false;
		for (int i = 0; i < spheres.size(); i++) {
			Object* Sloc = spheres[i];
			Vector Ploc, Nloc, colorLoc;
			if (Sloc->intersect(r, Ploc, Nloc, colorLoc)) {
				double dist_cam = (Ploc - r.C).norm2();
				if (dist_cam < dist_min) {
					dist_min = dist_cam;
					P = Ploc;
					N = Nloc;
					S = Sloc;
					color = colorLoc;
					//std::cout << color << endl;
				}
				found = true;
			}
		}
		return found;
	}

	Vector getColor(Ray& ray, int bounce) {
		Vector rayColor;
		Vector P, N, color;
		Object* S;
		double epsilon = 0.01;
		if (bounce <= 0) {
			return Vector(0, 0, 0);
		}
		if (intersect(ray, P, N, S, color)) {
			// gestion de la réflection
			if (S->mirror) {
				Vector R = ray.u - 2.*dot(ray.u, N)*N;
				Ray rReflect(P + epsilon * N, R);
				rayColor = getColor(rReflect, bounce - 1);
			}// gestion de la transparence
			else if (S->transparency) {
				double n1 = 1.;
				double n2 = 1.4;
				double dotP = dot(ray.u, N);
				if (dotP > 0) {
					std::swap(n1, n2);
					N = -N;
					dotP = -dotP;
				}
				double refrac = n1 / n2;
				Vector T;
				Vector Tt = refrac * (ray.u - dotP * N);
				double sqrt_val = 1 - sqr(refrac) * (1 - sqr(dotP));
				if (sqrt_val < 0) {
					T = ray.u - 2.*dot(ray.u, N) * N;
					N = -N;
				}
				else {
					Vector Tn = -sqrt(sqrt_val) * N;
					T = Tt + Tn;
				}
				Ray rayT(P - epsilon * N, T);
				return getColor(rayT, bounce);
			} // Eclairage diffus
			else if (S->light) {
				rayColor = Vector(0, 0, 0);
				// rayColor = S->albedo * 20000.;
				// rayColor = S.I * (S.albedo) / (4 * std::pow(M_PI * S.R, 2.));
			}
			else {
				Vector dirAleatoireLight = dirAleatoire(P - light.O);
				Vector x_prime = light.O + dirAleatoireLight * light.R;
				Vector PL = x_prime - P;
				double dist_lum = PL.norm2();
				PL.normalize();

				// A-t'on un obstacle entre la sphere et la lumiere
				Ray rayon_lumière(P + epsilon * PL, PL);
				Vector Pprime, Nprime, colorPrime;
				Object* Sprime;
				Vector direct;
				bool obstacle = intersect(rayon_lumière, Pprime, Nprime, Sprime, colorPrime);
				// Contribution directe
				if (obstacle && (P - Pprime).norm2() < dist_lum) {
						direct = Vector(0, 0, 0);
				}
				else {
					Vector Ox_prime = x_prime - light.O;
					Ox_prime.normalize();
					Vector Ox = P - light.O;
					Ox.normalize();
					direct = (light.I / (4 * M_PI)) *
						(color / M_PI) *
						dot(PL, N) * dot(-PL, Ox_prime) / (dot(Ox, Ox_prime) * dist_lum);
				}
				//std::cout << color << endl;
				Vector directionAleatoire = dirAleatoire(N);
				Ray rReflectAlea(P + epsilon * N, directionAleatoire);
				rayColor = direct;
				rayColor += getColor(rReflectAlea, bounce - 1) * color;
			}
		}
		else {
			rayColor = Vector(0, 0, 0);
		}
		return rayColor;
	}

	Vector getPerpendicular(Vector N) {
		Vector T1;
		if (abs(N.x) <= abs(N.y) && abs(N.x) <= abs(N.z)){
			T1 = Vector(0., -N.z, N.y);
		}
		else{
			if (abs(N.y) <= abs(N.x) && abs(N.y) <= abs(N.z)){
				T1 = Vector(-N.z, 0., N.x);
			}
			else{
				T1 = Vector(-N.y, N.x, 0.);
			}
		}
		T1.normalize();
		return T1;
	}

	Vector dirAleatoire(Vector N) {
		N.normalize();
		double r1 = uniform(engine);
		double r2 = uniform(engine);
		Vector directionAleatoireLocale(
			cos(2 * M_PI * r1) * sqrt(1 - r2),
			sin(2 * M_PI * r1) * sqrt(1 - r2),
			sqrt(r2)
		);
		Vector t1 = getPerpendicular(N);
		Vector t2 = cross(N, t1);
		Vector directionAleatoireGlobale = directionAleatoireLocale.z * N +
			directionAleatoireLocale.x * t1 +
			directionAleatoireLocale.y * t2;
		// directionAleatoireGlobale.normalize();
		return  directionAleatoireGlobale;
	}

	std::vector<Object*> spheres;
	Sphere light;
};


int main() {
	high_resolution_clock::time_point time1 = high_resolution_clock::now();
	int W = 512;
	int H = 512;
	// Camera
	double fov = 60 * M_PI / 180.;
	Vector C(0, 0, 0);
	const int focal = 25;
	std::vector<unsigned char> image(W*H * 3, 0);

	// Light
	Sphere sLight(Vector(-10, 20, 40), 10, Vector(1., 1., 1.));
	sLight.light = true;
	//sLight.I = 500000000;
	sLight.I = 20000000000;
	/*Light L(
		sLight,
		10000000000
	);*/

	// Scene
	Sphere s1(Vector( 15, 5, -33), 5, Vector(1, 1, 1));
	Sphere s10(Vector(15, 5, -33), 4, Vector(1., 1., 1));
	s1.transparency = true;
	s10.transparency = true;
	Sphere s2(Vector(-15, 0, -30), 8, Vector(0.5, 1., 0.));
	s2.mirror = true;
	Sphere s9(Vector(10, -10, -25), 3, Vector(0., 1., 0.5));
	Sphere s3(Vector(0, -2000 - 15, 0), 2000, Vector(1, 1, 0)); // sol
	Sphere s4(Vector(0,  2100 + 100,0), 2000, Vector(1, 1, 0)); // plafond
	Sphere s5(Vector(-2100 - 50, 0, 0), 2000, Vector(0, 1, 0)); // gauche
	Sphere s6(Vector( 2100 + 50, 0, 0), 2000, Vector(0, 0, 1)); // droite
	Sphere s7(Vector(0, 0, -2100 - 200), 2000, Vector(0, 1, 1));  //fond
	Sphere s8(Vector(0, 0, 1500), 1000 - 60, Vector(0, 1., 1.));
	Triangle t1(Vector(15, 5, -30), Vector(15, 0, -30), Vector(0, 15, -20));
	t1.albedo = Vector(0, 0.8, 1);
	// Loading a geometry and transforming it
	Geometry femme;
	std::cout << "Loading the Girl" << endl;
	femme.readOBJ("Beautiful_Girl\\girl.obj", true);
	femme.add_textures("Beautiful_Girl\\visage.bmp");
	femme.add_textures("Beautiful_Girl\\cheveux.bmp");
	femme.add_textures("Beautiful_Girl\\corps.bmp");
	femme.add_textures("Beautiful_Girl\\pantalon.bmp");
	femme.add_textures("Beautiful_Girl\\accessoires.bmp");
	femme.add_textures("Beautiful_Girl\\mains.bmp");
	// scale ~ 15, translate -10, inverse y et z
	std::cout << "Transforming the Girl" << endl;
	femme.transform(Vector(0, -10, -25), 10, true);
	femme.initBVH();
	// femme.albedo = Vector(1, 1, 1);
	std::cout << femme.root.bb.vmin << endl;
	std::cout << femme.root.bb.vmax << endl;

	// add objects to the scene
	Scene scene(sLight);
	scene.light = sLight;
	scene.addObject(&sLight);
	//scene.addObject(&t1);
	scene.addObject(&femme);
	//scene.addObject(&ironMan);
	scene.addObject(&s1);
	scene.addObject(&s10);
	scene.addObject(&s2);
	scene.addObject(&s3);
	scene.addObject(&s4);
	scene.addObject(&s5);
	scene.addObject(&s6);
	scene.addObject(&s7);
	scene.addObject(&s8);
	scene.addObject(&s9);

	double epsilon = 0.001;
	int max_bounce = 5;
	const int nrays = 30;
	std::cout << "Creating the pixels" << endl;
#pragma omp parallel for
	for (int i = 0; i < W; i++) {
		std::cout << i << endl;
		for (int j = 0; j < H; j++) {
			Vector pixelColor;
			Vector u = Vector(j - W / 2, -i + H / 2, -W / (2 * tan(fov / 2)));
			for (int k = 0; k < nrays; k++) {
				// Anti aliasing
				double r1 = uniform(engine);
				double r2 = uniform(engine);
				Vector delta(
					cos(2 * M_PI * r1) * sqrt(-2 * log(r2)) * 0.3,
					sin(2 * M_PI * r1) * sqrt(-2 * log(r2)) * 0.3,
					0
				);
				Vector direction = u + delta;
				direction.normalize();

				// focal blur
				double r3 = uniform(engine);
				double r4 = uniform(engine);
				Vector delta_blur(
					cos(2 * M_PI * r3) * sqrt(-2 * log(r4))*0.5,
					sin(2 * M_PI * r3) * sqrt(-2 * log(r4))*0.5,
					0
				);
				Vector focalPoint = C + direction * focal;
				Vector C_prime = C + delta_blur;
				Vector u_prime = focalPoint - C_prime;
				u_prime.normalize();
				Ray r(C_prime, u_prime);
				// pixel color
				pixelColor += scene.getColor(r, max_bounce) /nrays;
			}
			image[(i*W + j) * 3 + 0] = std::min(255., std::pow(pixelColor.x, 0.45));
			image[(i*W + j) * 3 + 1] = std::min(255., std::pow(pixelColor.y, 0.45));
			image[(i*W + j) * 3 + 2] = std::min(255., std::pow(pixelColor.z, 0.45));
		}
	}
	std::cout << "Writing image" << endl;
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
	std::cout << "Finished writing image" << endl;
	high_resolution_clock::time_point time2 = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(time2 - time1).count();
	std::cout << "Le script s'est execute pendant " << duration << " secondes" << endl;
	return 0;
}
